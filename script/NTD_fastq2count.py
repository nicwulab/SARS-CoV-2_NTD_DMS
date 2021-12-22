#!/usr/bin/python
import os
import sys
import glob
from Bio import SeqIO
from collections import Counter

def ProcessMultilib(infile, sampleID, fragID, primer_F, primer_R, var_dict):
    print ("Reading %s" % infile)
    records = SeqIO.parse(infile,"fastq")
    variants = [] 
    record_count = 0
    for record in records:
        record_count += 1
        Rseq  = record.seq
        Rroi = Rseq
        # Only include reads with correct forward/reverse primers and the correct length
        if Rroi[0:len(primer_F)] == primer_F and Rroi[-len(primer_R):] == primer_R:
            Rroi = Rroi[len(primer_F):-len(primer_R)] # Trim forward and reverse primers
            if len(Rroi) == 432:
                variants.append(Rroi)
        #if record_count == 1000: break
    var_dict[sampleID] = Counter(variants)
    return var_dict 

def writing_file(outfile, var_dict, sampleIDs, muts, fragID):
    print ("Compiling results into %s" % outfile)
    outfile = open(outfile,'w')
    sampleIDs_out = [sampleID[0:-2] for sampleID in sampleIDs]
    outfile.write('variant'+"\t"+'fragID'+"\t"+"\t".join(sampleIDs_out)+"\n")
    for mut in muts:
        out = [mut, fragID]
        for sampleID in sampleIDs:
            count = 0 if mut not in var_dict[sampleID].keys() else var_dict[sampleID][mut]
            out.append(count)
        outfile.write("\t".join(map(str,out))+"\n")
    outfile.close()

def Output(var_dict, outfile_A, outfile_B):
    sampleIDs_A = sorted([sampleID for sampleID in list(var_dict.keys()) if sampleID[-1]=='A'])
    sampleIDs_B = sorted([sampleID for sampleID in list(var_dict.keys()) if sampleID[-1]=='B'])
    muts_A    = sorted(list(set([mut for sampleID in sampleIDs_A for mut in var_dict[sampleID].keys()])))
    muts_B    = sorted(list(set([mut for sampleID in sampleIDs_B for mut in var_dict[sampleID].keys()])))
    writing_file(outfile_A, var_dict, sampleIDs_A, muts_A, 'A')
    writing_file(outfile_B, var_dict, sampleIDs_B, muts_B, 'B')

def main():
    outfile_A = 'result/NTD_DMS_count_nuc_A.tsv'
    outfile_B = 'result/NTD_DMS_count_nuc_B.tsv'
    infiles = glob.glob('fastq_merged/*.assembled.fastq')
    var_dict = {}
    for infile in infiles:
        sampleID = '_'.join(infile.rsplit('/')[1].rsplit('-')[0].rsplit('_')[0:-1])
        fragID   = sampleID[-1]
        if fragID == 'A':
            primer_F = 'CTGCTGCCTCTGGTGTCCAGC'
            primer_R = 'CGGGTGTACAGCAGCGCCAAC'
        if fragID == 'B':
            primer_F = 'AGCTGGATGGAAAGCGAGTTC'
            primer_R = 'ACCCTGAAGTCCTTCACCGTG'
        var_dict = ProcessMultilib(infile, sampleID, fragID, primer_F, primer_R, var_dict)
    Output(var_dict, outfile_A, outfile_B)

if __name__ == "__main__":
  main()
