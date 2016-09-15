#!/usr/bin/env python
from __future__ import print_function
from Bio import SeqIO
import vcf
import os
from itertools import tee

'''
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mass_auto_kar3_sorted.nodups.bam
'''
def get_vcf_readers(path):
    vcf_readers = []
    vcf_files = [i for i in os.listdir(path) if i.endswith('.vcf')]
    for vcf_file in vcf_files:
        vcf_path = os.path.join(path, vcf_file)
        vcf_reader = vcf.Reader(open(vcf_path), 'r')
        vcf_readers.append((os.path.basename(vcf_path), vcf_reader))
    return vcf_readers


def neighborhood(iterable):
    iterator = iter(iterable)
    prev = (None, None)
    item = iterator.next()  # throws StopIteration if empty.
    for next in iterator:
        yield (prev,item,next)
        prev = item
        item = next
    yield (prev,item,(None, None))


def walk_together(contig_name_len_gen, vcf_readers, shared_pos=4):
    contig_name_len_gen, contig_name_len_gen_2 = tee(contig_name_len_gen)
    #vcf_writer = vcf.Writer(open('merged_vcf_4.vcf', 'w+'), vcf_readers[0][1])
    vcf_writer = open('merged_vcf_4_rewrite_2.vcf', 'w+')
    for contig, seq_len in contig_name_len_gen_2:
        vcf_writer.write("##contig=<ID={},length={}>\n".format(contig, seq_len))

    line = '#CHROM\tPOS\tID\tREF\t{}\n'.format("\t".join([i for i,j in vcf_readers]))
    vcf_writer.write(line)

    for prev, curr, nex in neighborhood(contig_name_len_gen):

        p_c, p_l = prev
        contig, seq_len = curr
        n_c, n_l = nex
	
        try:
            int_ctg = int(contig[3:])
	except ValueError:
            int_ctg = int(1)
        contig_variants = {"POSITIONS":{},
                           "RECORDS":{fn:{} for fn,rd in vcf_readers}}


        for file_name, vcf_reader in vcf_readers:
            in_contig = False
            for item in vcf_reader:
                try:
                    if int(item.CHROM[3:]) > int(int_ctg):
                        break
                except ValueError:
                    pass
                if in_contig and item.CHROM != contig:
                    break
                if item.CHROM == contig:
                    in_contig = True
                if in_contig:
                    contig_variants["RECORDS"][file_name][item.POS] = item
                    if item.POS not in contig_variants["POSITIONS"]:
                        contig_variants["POSITIONS"][item.POS] = 0
                    contig_variants["POSITIONS"][item.POS] += 1

        for key in contig_variants['POSITIONS']:
            if contig_variants["POSITIONS"][key] >= shared_pos:
                #print (contig, key, contig_variants["POSITIONS"][key])
                alts = []
                ref = None
                for file_name, vcf_reader in vcf_readers:
                    if key in contig_variants["RECORDS"][file_name]:
                        ref = contig_variants["RECORDS"][file_name][key].REF
                if ref:
                    for file_name, vcf_reader in vcf_readers:
                        if key in contig_variants["RECORDS"][file_name]:
                            record = contig_variants["RECORDS"][file_name][key]
                            alts.append(",".join([str(i) for i in record.ALT]))
                        else:
                            alts.append(ref)

                line = "{}\t{}\t{}\t{}\n".format(contig, key, str(record.REF), "\t".join([str(i) for i in alts]))
                vcf_writer.write(line)


def contig_name_len_gen(fasta_path):
    with open(fasta_path) as input_handle:
        for record in SeqIO.parse(input_handle, 'fasta'):
            yield (record.id, len(record.seq))


def main():
    path = "/Users/genetics/Desktop/vcf_intersect"
    fasta_path = "/Users/genetics/Desktop/vcf_intersect/masurca_mito_y_x_removed.final.contigs.fasta"
    contig_gen = contig_name_len_gen(fasta_path)
    readers = get_vcf_readers(path)
    walk_together(contig_gen, readers)


main()
