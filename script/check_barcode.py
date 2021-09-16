#!/usr/bin/python
import os
import sys
import glob
import itertools
from Bio import SeqIO

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "NNK":"X"}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def extract_3rd_nucleotide(seq):
  nucleotide_3rd = ''
  for n in range(len(seq)):
    if (n + 1) % 3 == 0:
      nucleotide_3rd += seq[n]
  return nucleotide_3rd

def barcode_permutation(nucleotide_3rd):
  barcode_perms = []
  for i in nucleotide_3rd:
    if i == 'K':
      barcode_perms.append(['G','T'])
    else: 
      assert i in ['A','C','T','G']
      barcode_perms.append([i])
  barcode_list = []
  for barcode_perm in itertools.product(*barcode_perms):
    barcode_list.append(''.join(barcode_perm))
  return barcode_list
  
def extract_barcodes(infile, outfile):
  records = SeqIO.parse(infile,"fasta")
  outfile = open(outfile,'w')
  barcode_dict = {}
  outfile.write("barcode"+"\t"+"position"+"\n")
  for record in records:
    ID  = record.id
    seq = record.seq
    seq = seq.upper()
    nucleotide_3rd = extract_3rd_nucleotide(seq)
    barcode_list = barcode_permutation(nucleotide_3rd)
    barcode_dict[ID] = barcode_list
    print ("Barcode:"+"\t"+ID+"\t"+','.join(barcode_list))
    for barcode in barcode_list:
      outfile.write(barcode+"\t"+ID+"\n")
  outfile.close()
  return (barcode_dict)

def compare_barcodes(barcode_dict):
  ID_list = list(barcode_dict.keys())
  overlap_count = 0
  for n in range(len(ID_list)):
    for m in range(len(ID_list)):
      if n < m:
        ID1 = ID_list[n]
        ID2 = ID_list[m]
        ID1_barcodes = set(barcode_dict[ID1])
        ID2_barcodes = set(barcode_dict[ID2])
        overlap_barcodes = ID1_barcodes.intersection(ID2_barcodes)
        if len(overlap_barcodes) != 0:
          overlap_count += 1
          print ("Overlapping barcodes between %s and %s:" % (ID1, ID2))
          print (overlap_barcodes)
  if overlap_count == 0: 
    print ("Good job! No overlapping barcodes!")

def main():
  infile  = "primer/SARS2-NTD_lib_Fprimer_bc.fa"
  outfile = "data/barcodes.tsv"
  barcode_dict = extract_barcodes(infile, outfile)
  compare_barcodes(barcode_dict)

if __name__ == "__main__":
  main()
