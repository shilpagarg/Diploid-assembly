import vcf
import sys
from collections import defaultdict
# get haplotigs based on vcf and canu contigs
input_file = sys.argv[1]
canu_contig_file = sys.argv[2]

vcf_reader = vcf.Reader(open(input_file, 'r', encoding='utf-8'))
block_ids = set()

for record in vcf_reader:
    try:
        block_ids.add(int(record.samples[0]['PS'])-1)
    except:
        continue

canu_contig_ref = list()
sorted_block_ids = sorted(block_ids)
del sorted_block_ids[0]
with open(canu_contig_file) as fp:
	for line in fp:
		var=line.rstrip()[0]
		if var!='>':
                    for x in list(line.rstrip()):
                        canu_contig_ref.append(x)

canu_contig_alt = list()
for a in canu_contig_ref:
	canu_contig_alt.append(a)

#print(len(block_ids))
vcf_writer = open(input_file+ '_linear' + '.fasta', 'w')
hap1_seq = ''
hap2_seq = ''
tmp =0
#print(canu_contig_alt)

hap1_dist = defaultdict()
hap2_dist = defaultdict()

vcf_reader = vcf.Reader(open(input_file, 'r', encoding='utf-8'))
for record in vcf_reader:
	try:

		pos=int(record.POS)
		ref=record.REF
		alt=record.ALT[0]
		#print('hello')
		#print(pos)
		hap1=record.samples[0]['GT']
		allele1=hap1.split('|')
		if int(allele1[0])==0:
			hap1_dist[pos-1]=str(ref)
		else:
			hap1_dist[pos-1]=str(alt)

		if int(allele1[1])==0:
			hap2_dist[pos-1]=str(ref)
		else:
			hap2_dist[pos-1]=str(alt)				
	except:
		continue 

count =0
for x in range(0,len(canu_contig_ref)):
	if x in hap1_dist:
		hap1_seq = hap1_seq+ hap1_dist[x]
		hap2_seq = hap2_seq+ hap2_dist[x]
	else:
		hap1_seq = hap1_seq+ canu_contig_ref[x]
		hap2_seq = hap2_seq+ canu_contig_ref[x]
	if x in sorted_block_ids:
		vcf_writer.write(">" + input_file+ str(count) + "_1" + "\n")
		vcf_writer.write(hap1_seq + "\n")
		vcf_writer.write(">" + input_file+ str(count) + "_2" + "\n")
		vcf_writer.write(hap2_seq + "\n")
		count=count+1
		hap1_seq = ''
		hap2_seq = ''

vcf_writer.write(">" + input_file+ str(count) + "_1" + "\n")
vcf_writer.write(hap1_seq + "\n")
vcf_writer.write(">" + input_file+ str(count) + "_2" + "\n")
vcf_writer.write(hap2_seq + "\n")





