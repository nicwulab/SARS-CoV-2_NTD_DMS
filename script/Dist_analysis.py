#!/usr/bin/python
import sys
from prody import *
from scipy.spatial import distance 

def extract_resi_centroid(p):
  centroid_dict = {}
  for chain in sorted(list(set(p.getChids()))):
    p_chain = p.select('chain '+chain)
    for resnum in sorted(list(set(p_chain.getResnums()))):
      res = p_chain.select('resnum '+str(resnum))
      centroid = (calcCenter(res))
      centroid_dict[chain+'-'+str(resnum)] = centroid
  return (centroid_dict)

def min_distance(A_centroid_dict, B_centroid_dict):
  A_min_dist_to_B = {}
  for A_res in A_centroid_dict.keys():
    A_res_centroid = A_centroid_dict[A_res]
    min_dist = 9999
    for B_res in B_centroid_dict.keys():
      B_res_centroid = B_centroid_dict[B_res]
      dist = distance.euclidean(A_res_centroid, B_res_centroid)
      if dist < min_dist: min_dist = dist
    A_min_dist_to_B[A_res] = min_dist
  return (A_min_dist_to_B)

def writing_dist_file(outfile, NTD_min_dist_to_RBD_dict, NTD_min_dist_to_S2_dict):
  print ("writing: %s" % outfile)
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['pos', 'dist_to_RBD', 'dist_to_S2', 'dist_to_RBD/S2'])+"\n")
  NTD_residues = list(set(list(NTD_min_dist_to_RBD_dict.keys())+list(NTD_min_dist_to_S2_dict.keys())))
  NTD_residues = sorted(list(map(lambda x:int(x.rsplit('-')[1]), NTD_residues)))
  for NTD_resi in NTD_residues:
    dist_to_RBD = NTD_min_dist_to_RBD_dict['X-'+str(NTD_resi)]
    dist_to_S2  = NTD_min_dist_to_S2_dict['X-'+str(NTD_resi)]
    outfile.write("\t".join(map(str,[NTD_resi, dist_to_RBD, dist_to_S2, min([dist_to_RBD, dist_to_S2])]))+"\n")
  outfile.close()

def main():
  outfile = 'result/Dist_NTD_to_RBD_S2.tsv'
  p = parsePDB('PDB/spike_with_complete_NTD.pdb')
  NTD = p.select('chain X resnum 14 to 301')
  RBD = p.select('chain A B C resnum 331 to 528')
  S2  = p.select('chain A B C resnum 685 to 1146')
  NTD_centroid_dict = extract_resi_centroid(NTD)
  RBD_centroid_dict = extract_resi_centroid(RBD)
  S2_centroid_dict = extract_resi_centroid(S2)
  NTD_min_dist_to_RBD_dict = min_distance(NTD_centroid_dict, RBD_centroid_dict)
  NTD_min_dist_to_S2_dict  = min_distance(NTD_centroid_dict, S2_centroid_dict)
  writing_dist_file(outfile, NTD_min_dist_to_RBD_dict, NTD_min_dist_to_S2_dict)
  

if __name__ == "__main__":
  main()
