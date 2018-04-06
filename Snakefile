import pysam
picard_tmpdir_switch=''
illumina_files = ['fAnaTes1_S1_L008', 'fAnaTes1_S2_L008', 'fAnaTes1_S3_L008', 'fAnaTes1_S4_L008']
refs = ['fAnaTes1.1']
pacbio_file = ['m54097_170223_195514', 'm54097_170225_002438', 'm54097_170225_103404', 'm54097_170225_204342', 'm54097_170226_065311', 'm54097_170226_170255', 'm54097_170227_031236', 'm54097_170227_112229']

# Tools assumed to be installed somewhere on the PATH.
samtools = 'samtools'
bwa = 'bwa'
picard = 'picard'
seqtk = 'seqtk'
splitVCFbyblocks = 'scripts/splitVCFbyblocks.py'
whatshap = 'scripts/whatshap/bin/whatshap'
splitVCFbyblocks = 'splitVCFbyblocks.py'

rule master:
	input:
		expand('phased/{refs}/whatshap.vcf.gz', refs=refs),
	message: 'MASTER rule'

# make sure illumina data, pacbio and reference in their respective directories.

rule index_reference:
	output:
		'ref/{refs}.fa.gz.bwt'
	input:
		'ref/{refs}.fa.gz'
	shell:
		"bwa index {input}"

rule tenX_to_illumina:
	input: '10X/{file}_R1_001.fastq.gz', '10X/{file}_R2_001.fastq.gz', 
	output: 'illumina/{file}_R1_001.fastq.gz', 'illumina/{file}_R2_001.fastq.gz'
	shell: 'seqtk trimfq -b 23 t12 {input[0]} | gzip > {output[0]} && seqtk trimfq -e 23 t12 {input[1]} | gzip > {output[1]}'

rule align_illumina:
	input: 'illumina/{file}_R1_001.fastq.gz', 'illumina/{file}_R2_001.fastq.gz', 'ref/{refs}.fa.gz', expand('ref/{refs}.fa.gz.bwt', refs=refs),
	output: 'illumina/{refs}/{file}.bam'
	shell: 'bwa mem -t 12 {input[2]} {input[0]} {input[1]} | samtools view -Sb - > {output}'

rule merge_illumina:
	input: expand('illumina/{refs}/{file}.bam', file=illumina_files, refs=refs),
	output: 'illumina/{refs}/illumina.bam'
	shell: 'samtools merge {output} {input[0]}'

rule sort_illumina_bam:
	input: 'illumina/{refs}/illumina.bam'
	output: 'illumina/{refs}/illumina.sorted.bam'
	log:'illumina/{refs}/illumina.sorted.bam.log'
	shell: 'picard -Xmx16g SortSam VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=50000 SORT_ORDER=coordinate CREATE_INDEX=true CREATE_MD5_FILE=true I={input} O={output} > {log} 2>&1'

rule call_variants:
	input: 'illumina/{refs}/illumina.sorted.bam', 'ref/{refs}.fa.gz'
	output: 'illumina/{refs}/illumina.vcf'
	shell: 'freebayes -f {input[1]} {input[0]} | vcf-sort  > {output}'

rule align_pacbio:
	input:'PacBio/{file}.subreads.fasta.gz','ref/{refs}.fa.gz', expand('ref/{refs}.fa.gz.bwt', refs=refs)
	output: 'pacbio/{refs}/{file}.bam'
	shell: 'bwa mem -t 12 {input[1]} {input[0]} | samtools view -Sb - > {output}'

rule merge_pacbio:
	input: expand('pacbio/{refs}/{file}.bam', file=pacbio_file, refs=refs),
	output: 'pacbio/{refs}/pacbio.bam'
	shell: 'samtools merge {output} {input}'

rule sort_pacbio_bam:
	input: 'pacbio/{refs}/pacbio.bam'
	output: 'pacbio/{refs}/pacbio.sorted.bam'
	log:'pacbio/{refs}/pacbio.sorted.bam.log'
	shell: 'picard -Xmx16g SortSam VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=50000 SORT_ORDER=coordinate CREATE_INDEX=true CREATE_MD5_FILE=true I={input} O={output} > {log} 2>&1'

rule run_whatshap:
	input: 'illumina/{refs}/illumina.vcf', 'pacbio/{refs}/pacbio.sorted.bam', 'ref/{refs}.fa.gz'
	output: 'phased/{refs}/whatshap.vcf'
	log: 'phased/{refs}/whatshap.vcf.log'
	shell:'{whatshap} phase --max-coverage 15 --reference {input[2]} --ignore-read-groups {input[0]} {input[1]} --output {output} > {log} 2>&1'

rule zip_vcf:
	input: 'phased/{refs}/whatshap.vcf'
	output: 'phased/{refs}/whatshap.vcf.gz'
	shell: 'bgzip -c {input[1]} > {output} && tabix -p vcf {output}'

chrs= []
import gzip
import re
for ref in refs:
	with gzip.open("phased/%s/whatshap.vcf.gz" % ref, "rb") as ifile:
		for line in ifile:
			line=line.decode()
			# look for ##contig=<ID=Contig99arrow,length=72881>
			if not re.match("^##", line):
				break
			if re.match("^##contig", line):
				chrid=line.split("=")[2].split(",")[0]
				chrs.append("phased/%s/split/%s.vcf" %(ref,chrid))

rule separate_contigs:
	input:'ref/{refs}.fa.gz'
	output:'phased/{refs}/split/{chrid}.fa'
	shell: '(mkdir -p phased/{refs} && cd phased/{refs} && faidx -x ../../{input[0]})'

rule master2:
	input: chrs
	message: "Split"

rule split_vcf:
	input: 'phased/{refs}/whatshap.vcf.gz'
	output: 'phased/{refs}/split/{chrid}.vcf'
	shell: 'tabix {input} -h {wildcards.chrid} > {output}'

rule assemblies:
	input:'phased/{refs}/split/{chrid}.vcf', 'phased/{refs}/split/{chrid}.fa'
	output:'phased/{refs}/split/{chrid}.vcf__linear.fasta'
	shell: 'python3 {splitVCFbyblocks}' {input[1]} {input[0]}'

