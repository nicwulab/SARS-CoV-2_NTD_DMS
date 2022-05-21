#!/usr/bin/python
import os
import sys
import glob
import pandas as pd
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

def main():
  outfile  = 'result/RBD_exp_RSA.tsv' 
  RBD_data = pd.read_csv('data/RBD_DMS.tsv', sep="\t")
  RBD_data = RBD_data.groupby('site_SARS2')['expr_avg'].mean().reset_index()
  file_asa  = 'data/ASA.table'
  dict_asa  = reading_ASA(file_asa)
  dssp_RBD  = dssp_dict_from_pdb_file('PDB/6zge_RBD.pdb', DSSP='mkdssp')[0]
  chainID   = 'A'
  positions = list(range(331, 532))
  RSA_dict_RBD = parse_dssp(dssp_RBD, dict_asa, chainID, positions)
  RBD_data['RSA'] = RBD_data['site_SARS2'].map(RSA_dict_RBD)
  RBD_data.rename(columns = {'site_SARS2':'pos'}, inplace = True)
  print ('writing: %s' % outfile)
  RBD_data.to_csv(outfile, sep="\t", index=False)
   
if __name__ == "__main__":
  main()
