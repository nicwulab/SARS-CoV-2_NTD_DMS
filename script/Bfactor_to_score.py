#!/usr/bin/python
import os
import sys
from collections import defaultdict

def read_expfile(file_file):
  infile = open(file_file, 'r')
  exp_dict = {}
  low_count_pos = []
  for line in infile.readlines():
    if 'pos' in line: continue
    info = line.rstrip().rsplit("\t")
    if len(info) < 4:
      resi, pos = info
      low_count_pos.append(pos)
    else:
      resi, pos, count, exp = info
      if float(count) >= 6: 
        exp_dict[pos] = float(exp)
      else:
        low_count_pos.append(pos)
  infile.close()
  print ('Positions with low counts:\n'+'+'.join(low_count_pos))
  return (exp_dict)

def normalizing_exp(exp_dict):
  max_exp  = max(exp_dict.values())
  min_exp  = min(exp_dict.values())
  print ('exp range (pre-normalized): %f to %f' % (min_exp, max_exp))
  for pos in exp_dict.keys():
    exp = exp_dict[pos]
    exp = (exp-min_exp)/(max_exp-min_exp)*100
    exp_dict[pos] = exp
  print ('exp range (post-normalized): %f to %f' % (min(exp_dict.values()), max(exp_dict.values())))
  return (exp_dict)

def add_exp_to_pdb(exp_dict, pdb_file):
  assert('.pdb' in pdb_file)
  new_pdb_file = pdb_file.replace('.pdb', '_exp.pdb')
  print ("writing: %s" % new_pdb_file)
  infile  = open(pdb_file, 'r')
  outfile = open(new_pdb_file,'w')
  for line in infile.readlines():
    if "ATOM" == line[0:4]:
      pos      = int(line[23:26])
      chain    = line[21:23].rstrip()
      b_factor = line[60:66]
      if chain == 'A' and str(pos) in exp_dict.keys():
        exp = str(round(exp_dict[str(pos)],2))
        exp = exp+'0'*(2-len(exp.rsplit('.')[-1]))
        exp = ' '*(6-len(exp))+exp
        new_line = line[0:60]+exp+line[66::]
        outfile.write(new_line)
      else:
        outfile.write(line)
    else:
      outfile.write(line)
  outfile.close()
  

def main():
  pdb_file = "PDB/7b62.pdb"
  exp_file = 'result/NTD_DMS_scores_by_resi.tsv'
  exp_dict = read_expfile(exp_file)
  exp_dict = normalizing_exp(exp_dict)
  add_exp_to_pdb(exp_dict, pdb_file)

if __name__ == "__main__":
  main() 
