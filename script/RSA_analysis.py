#!/usr/bin/python
#Required dssp_parser: http://openwetware.org/wiki/Wilke:ParseDSSP
import os
import sys
import glob
from collections import defaultdict
from Bio.PDB.DSSP import dssp_dict_from_pdb_file

def reading_ASA(file_asa):
  infile = open(file_asa,'r')
  dict_asa = {}
  for line in infile.readlines():
    if 'aa' in line: continue
    line = line.rstrip().rsplit("\t")
    aa  = line[0]
    asa = line[3]
    dict_asa[aa] = asa
  return dict_asa

def parse_dssp(dssp, dict_asa, chainID, positions):
  RSA_dict = {}
  for pos in positions:
    ID = (chainID,(' ', pos, ' '))
    if ID not in dssp.keys(): continue
    dssp_info = list(dssp[ID])
    aa = dssp_info[0]
    if aa.islower(): aa = 'C'
    SA = dssp_info[2]
    RSA = float(SA)/float(dict_asa[aa])
    RSA_dict[pos] = RSA
  return RSA_dict

def classify_pos_type(pos, RSA_monomer, delta_RSA, site_dict):
  pos_type = []
  for site in site_dict.keys():
    if int(pos) in site_dict[site]:
      pos_type.append(site)
  if len(pos_type) == 0:
    if RSA_monomer < 0.05:
      pos_type.append('buried')
    elif delta_RSA > 0.5*RSA_monomer:
      pos_type.append('interface')
  if len(pos_type) > 1:
    print ('position %s belongs to more than 1 classification (%s)' % (pos, ', '.join(pos_type)))
  if len(pos_type) == 0:
    pos_type.append('exposed')
  pos_type = '-'.join(pos_type)
  return pos_type

def compile_out(outfile, RSA_dict_monomer, RSA_dict_trimer, positions):
  print ('writing: %s' % outfile)
  outfile = open(outfile,'w')
  outfile.write("\t".join(['pos', 'RSA_monomer', 'RSA_trimer', 'delta_RSA'])+"\n")
  for pos in positions:
    RSA_monomer  = RSA_dict_monomer[pos]
    RSA_trimer   = RSA_dict_trimer[pos] if pos in RSA_dict_trimer.keys() else RSA_monomer
    delta_RSA    = RSA_monomer-RSA_trimer
    if delta_RSA < 0:
      delta_RSA = 0
    outfile.write("\t".join(map(str, [pos, RSA_monomer, RSA_trimer, delta_RSA]))+"\n")
  outfile.close()

def main():
  outfile   = 'result/NTD_RSA.tsv'
  dssp_monomer  = dssp_dict_from_pdb_file('PDB/7b62.pdb', DSSP='mkdssp')[0]
  dssp_trimer   = dssp_dict_from_pdb_file('PDB/6zge.pdb', DSSP='mkdssp')[0]
  file_asa  = 'data/ASA.table'
  dict_asa  = reading_ASA(file_asa)
  chainID   = 'A'
  positions = list(range(14, 302))
  RSA_dict_monomer  = parse_dssp(dssp_monomer, dict_asa, chainID, positions)
  RSA_dict_trimer   = parse_dssp(dssp_trimer, dict_asa, chainID, positions)
  compile_out(outfile, RSA_dict_monomer, RSA_dict_trimer, positions)

if __name__ == "__main__":
  main()
