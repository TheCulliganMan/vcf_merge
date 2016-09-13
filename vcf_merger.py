#!/usr/bin/env python
from __future__ import print_function
from Bio import SeqIO
import vcf
import os

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
    vcf_writer = vcf.Writer(open('merged_vcf_4.vcf', 'w+'), vcf_readers[0][1])
    for prev, curr, nex in neighborhood(contig_name_len_gen):

        p_c, p_l = prev
        contig, seq_len = curr
        n_c, n_l = nex

        int_ctg = int(contig[3:])

        contig_variants = {"POSITIONS":{},
                           "RECORDS":{fn:{} for fn,rd in vcf_readers}}


        for file_name, vcf_reader in vcf_readers:
            in_contig = False
            for item in vcf_reader:
                if int(item.CHROM[3:]) > int(int_ctg):
                    break
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
                alts = set()
                for file_name, vcf_reader in vcf_readers:
                    if key in contig_variants["RECORDS"][file_name]:
                        record = contig_variants["RECORDS"][file_name][key]
                        for i in record.ALT:
                            alts.add(str(i))
                record.ALT = [vcf.model._Substitution(i) for i in alts]
                vcf_writer.write_record(record)


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
