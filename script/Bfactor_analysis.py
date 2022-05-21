#!/usr/bin/python
import sys
from prody import *

def extract_resi_Bfactor(p):
  Bfactor_dict = {}
  p = p.select('name CA')
  for chain in sorted(list(set(p.getChids()))):
    p_chain = p.select('chain '+chain)
    for resnum in sorted(list(set(p_chain.getResnums()))):
      Bfactor = p_chain.select('resnum '+str(resnum)).getBetas()
      Bfactor_dict[chain+'-'+str(resnum)] = Bfactor[0]
  return (Bfactor_dict)

def writing_Bfactor(Bfactor_dict, chain, residues, outfile):
  print ("writing: %s" % outfile)
  outfile = open(outfile, 'w')
  outfile.write("pos"+"\t"+"Bfactor"+"\n")
  for residue in residues:
    resi_ID = chain+'-'+str(residue)
    if resi_ID in Bfactor_dict.keys():
      B = Bfactor_dict[chain+'-'+str(residue)]
      outfile.write("\t".join(map(str,[residue, B]))+"\n")
  outfile.close()

def main():
  outfile = 'result/NTD_Bfactor.tsv'
  p = parsePDB('PDB/6zge.pdb')
  chain = 'A'
  residues = list(range(14,302))
  Bfactor_dict = extract_resi_Bfactor(p)
  writing_Bfactor(Bfactor_dict, chain, residues, outfile)
  

if __name__ == "__main__":
  main()
