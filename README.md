## Studying the impact of NTD mutations on SARS-CoV-2 spike expression using deep mutational scanning

### Dependencies
* [Python](https://www.python.org/) (version 3.9)
* [Biopython](https://github.com/biopython/biopython)
* [R](https://www.r-project.org/) (version 4.1)
* [PEAR](https://github.com/tseemann/PEAR)
* [ProDy](http://prody.csb.pitt.edu/)

### Input files
* [./Fasta/SARS2-NTD.fa](./Fasta/SARS2-NTD.fa): NTD sequences with 21 nt upstream (5' flank)
* [./Fasta/Amplicon.fa](./Fasta/Amplicon.fa): Amplicon sequences for NGS (including recovery primer regions)
* [./Fasta/NTD_ref.fa](./Fasta/NTD_ref.fa): Reference (i.e. wild type) amino acid sequences (primer regions not included)
* [./PDB/6zge.pdb](./PDB/6zge.pdb): Spike cryoEM structure from [Wrobel et al. (2020)](https://www.nature.com/articles/s41594-020-0468-7)
* [./PDB/6zge_RBD.pdb](./PDB/6zge_RBD.pdb): RBD only from the Spike cryoEM structure from [Wrobel et al. (2020)](https://www.nature.com/articles/s41594-020-0468-7)
* [./PDB/7b62.pdb](./PDB/7b62.pdb): NTD crystal structure from [Rosa et al. (2021)](https://www.science.org/doi/10.1126/sciadv.abg7607)
* [./PDB/spike_with_complete_NTD.pdb](./PDB/spike_with_complete_NTD.pdb): 6zge with NTD in chain A being replaced by 7b62 (chain X)
* [./data/ASA.table](./data/ASA.table): Maximum accessible surface area for individual amino acids from [Tien et al. (2013)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0080635)
* [./data/RBD_DMS_data.csv](./data/RBD_DMS_data.csv): DMS data of RBD from [Chan et al. (2021)](https://www.science.org/doi/10.1126/sciadv.abf1738)
* Raw read files in fastq format from NIH SRA database [BioProject PRJNA792013](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA792013)

### Primer design for DMS library construction
1. Generating foward (NNK + internal barcode) and reverse primers (constant)   
``python3 script/lib_primer_design.py``
    - Input file:
      - [./Fasta/SARS2-NTD.fa](./Fasta/SARS2-NTD.fa)
    - Output files:
      - [./primer/SARS2-NTD_lib_Fprimer_bc.fa](./primer/SARS2-NTD_lib_Fprimer_bc.fa)
      - [./primer/SARS2-NTD_lib_Rprimer_bc.fa](./primer/SARS2-NTD_lib_Rprimer_bc.fa)

2. Generating barcode file   
``python3 script/check_barcode.py``
    - Input files:
      - [./primer/SARS2-NTD_lib_Fprimer_bc.fa](./primer/SARS2-NTD_lib_Fprimer_bc.fa)
      - [./Fasta/Amplicon.fa](./Fasta/Amplicon.fa)
    - Output file:
      - [./data/barcodes.tsv](./data/barcodes.tsv)

### Calculating experssion score from DMS data
1. Merge overlapping paired-end reads using [PEAR](https://github.com/tseemann/PEAR)   
``pear -f [FASTQ FILE FOR FORWARD READ] -r [FASTQ FILE FOR FORWARD READ] -o [OUTPUT FASTQ FILE]``   
    - Output files should be placed in the fastq_merged/ folder and named as described in [./doc/filename_merged_fastq.txt](./doc/filename_merged_fastq.txt)

2. Counting variants based on nucleotide sequences   
``python3 script/NTD_fastq2count.py``   
    - Input files:
      - Merged read files in fastq_merged/ folder
    - Output files:
      - result/NTD_DMS_count_nuc_A.tsv
      - result/NTD_DMS_count_nuc_B.tsv

3. Convert nucleotide sequences to amino acid mutations   
``python3 script/NTD_count_nuc2aa.py``   
    - Input files:
      - [./data/barcodes.tsv](./data/barcodes.tsv)
      - [./Fasta/NTD_ref.fa](./Fasta/NTD_ref.fa)
      - result/NTD_DMS_count_nuc_A.tsv
      - result/NTD_DMS_count_nuc_B.tsv
    - Output files:
      - result/NTD_DMS_count_aa_A.tsv
      - result/NTD_DMS_count_aa_B.tsv

4. Compute expression score   
``python3 script/NTD_count2score.py``   
    - Input files:
      - result/NTD_DMS_count_aa_A.tsv
      - result/NTD_DMS_count_aa_B.tsv
    - Output file:
      - [./result/NTD_DMS_scores.tsv](./result/NTD_DMS_scores.tsv)
      - [./result/NTD_DMS_scores_by_resi.tsv](./result/NTD_DMS_scores_by_resi.tsv)

5. Plot correlation between replicates as quality control  
``Rscript script/plot_QC.R``   
    - Input file:
      - [./result/NTD_DMS_scores.tsv](./result/NTD_DMS_scores.tsv)
    - Output files:
      - [./graph/QC_replicate_exp.png](./graph/QC_replicate_exp.png)
      - [./graph/Exp_by_class.png](./graph/Exp_by_class.png)
      - [./result/NTD_DMS_expression_score.tsv](./result/NTD_DMS_expression_score.tsv)
``Rscript script/plot_QC_mean_exp_by_replicates.R``   
    - Input file:
      - [./result/NTD_DMS_scores_by_resi.tsv](./result/NTD_DMS_scores_by_resi.tsv)
    - Output files:
      - [./graph/QC_replicate_mean_exp.png](./graph/QC_replicate_mean_exp.png)
``Rscript script/QC_bin0_vs_bin3_counts.R``
    - Input file:
      - [./result/NTD_DMS_scores.tsv](./result/NTD_DMS_scores.tsv)
    - Output files:
      - [./graph/bin0_bin_3_ratio_rep_1_all.png](./graph/bin0_bin_3_ratio_rep_1_all.png)
      - [./graph/bin0_bin_3_ratio_rep_2_all.png](./graph/bin0_bin_3_ratio_rep_2_all.png)
``Rscript script/QC_bin0_vs_bin3_biomodality.R``
    - Input file:
      - [./result/NTD_DMS_scores.tsv](./result/NTD_DMS_scores.tsv)
    - Output files:
      - [./graph/bin0_by_replicates.png](./graph/bin0_by_replicates.png)
      - [./graph/bin3_by_replicates.png](./graph/bin3_by_replicates.png)

6. Plot heatmap for input frequencies of individual mutations
``Rscript script/plot_input_freq_heatmap.R``   
    - Input file:
      - [./result/NTD_DMS_scores.tsv](./result/NTD_DMS_scores.tsv)
    - Ouput file:
      - [./graph/NTD_input_freq_heatmap.png](./graph/NTD_input_freq_heatmap.png)

7. Plot heatmap for the expression scores of individual mutations   
``Rscript script/plot_score_heatmap.R``   
    - Input file:
      - [./result/NTD_DMS_scores.tsv](./result/NTD_DMS_scores.tsv)
    - Ouput file:
      - [./graph/NTD_exp_heatmap.png](./graph/NTD_exp_heatmap.png)

8. Plot mutational tolerability in loops vs others   
``Rscript script/NTD_loop_vs_other_residues.R``   
    - Input file:
      - [./result/NTD_DMS_scores_by_resi.tsv](./result/NTD_DMS_scores_by_resi.tsv)
    - Output file:
      - [./graph/ntd_loop_vs_other_residues.png](./graph/ntd_loop_vs_other_residues.png)
    
9. Plot the mutational tolerability in selected hotspot regions
``Rscript script/Hot_spots_vs_other_residues.R``   
    - Input file:
      - [./result/NTD_DMS_scores_by_resi.tsv](./result/NTD_DMS_scores_by_resi.tsv)
    - Output file:
      - [./graph/hot_spots_vs_other_residues.png](./graph/hot_spots_vs_other_residues.png)
      
### Structural analysis
1. Computing relative solvent accessibility (RSA) for individual residues   
``python3 script/RSA_analysis.py``   
    - Input files:
      - [./PDB/7b62.pdb](./PDB/7b62.pdb)
      - [./PDB/6zge.pdb](./PDB/6zge.pdb)
      - [./data/ASA.table](./data/ASA.table)
    - Ouput file:
      - [./result/NTD_RSA.tsv](./result/NTD_RSA.tsv)

2. Compute expression score for RBD DMS data
``Rscript script/compute_exp_score_RBD.R``
      - Input file:
        - [./data/RBD_DMS.csv](./data/RBD_DMS.csv)
      - Output file:
        - [./data/RBD_DMS_exp.tsv](./data/RBD_DMS_exp.tsv)

3. Compute RSA for RBD residues
``python3 script/RBD_analysis.py``
    - Input files:
      - [./data/RBD_DMS_exp.tsv](./data/RBD_DMS_exp.tsv)
      - [./PDB/6zge_RBD.pdb](./PDB/6zge_RBD.pdb)
      - [./data/ASA.table](./data/ASA.table)
    - Output file:
      - [./result/RBD_exp_RSA.tsv](./result/RBD_exp_RSA.tsv)

4. Computing the distance of individual NTD residues to RBD or S2   
``python3 script/Dist_analysis.py``   
    - Input file:
      - [./PDB/spike_with_complete_NTD.pdb](./PDB/spike_with_complete_NTD.pdb)
    - Output file:
      - [./result/Dist_NTD_to_RBD_S2.tsv](./result/Dist_NTD_to_RBD_S2.tsv)

5. Replace the B-factor by expression score in the PDB file   
``python3 script/Bfactor_to_score.py``   
    - Input file:
      - [./PDB/7b62.pdb](./PDB/7b62.pdb)
      - [./result/NTD_DMS_scores_by_resi.tsv](./result/NTD_DMS_scores_by_resi.tsv)
    - Output file:
      - [./PDB/7b62_exp.pdb](./PDB/7b62_exp.pdb)

6. Plot mean expression score vs RSA for individual residues   
``Rscript script/plot_RSA.R``   
    - Input file:
      - [./result/NTD_DMS_scores_by_resi.tsv](./result/NTD_DMS_scores_by_resi.tsv)
      - [./result/NTD_RSA.tsv](./result/NTD_RSA.tsv)
    - Output file:
      - [./graph/Exp_vs_RSA.png](./graph/Exp_vs_RSA.png)

7. Plot mutational tolerability vs RSA for RBD DMS data   
``Rscript script/Mean_expression_score_in_mammalian_system.R``   
    - Input file:
      - [./result/RBD_DMS_data.csv](./result/RBD_DMS_data.csv)
    - Output file:
      - [./result/RBD_exp_RSA.tsv](./result/RBD_exp_RSA.tsv)
``Rscript script/plot_RSA_RBD.R``   
    - Input file:
      - [./result/RBD_exp_RSA.tsv](./result/RBD_exp_RSA.tsv)
    - Output file:
      - [./graph/Exp_vs_RSA_RBD.png](./graph/Exp_vs_RSA_RBD.png)

8. Plot mutational tolerability vs distance to RBD/S2 for individual residues and categorized by antibody epitopes   
``Rscript script/Dist_to_RBD_S2_exp_by_Ab.R``   
    - Input file:
      - [./result/NTD_DMS_scores_by_resi.tsv](./result/NTD_DMS_scores_by_resi.tsv)
      - [./result/Dist_NTD_to_RBD_S2.tsv](./result/Dist_NTD_to_RBD_S2.tsv)
    - Output file:
      - [./graph/Exp_vs_dist.png](./graph/Exp_vs_dist.png)
      - [./graph/antibody_epi_vs_mean_exp.png](./graph/antibody_epi_vs_mean_exp.png)
``Rscript script/antibody_epitopes_vs_distances.R``   
    - Input file:
      - [./result/NTD_DMS_scores_by_resi.tsv](./result/NTD_DMS_scores_by_resi.tsv)
      - [./result/Dist_NTD_to_RBD_S2.tsv](./result/Dist_NTD_to_RBD_S2.tsv)
    - Output file:
      - [./graph/Exp_vs_dist.png](./graph/Exp_vs_dist.png)
      - [./graph/antibody_epitopes_vs_distance.png](./graph/antibody_epitopes_vs_distance.png)

9. Plot mutational tolerability vs sequence conservation for individual residues   
``Rscript script/align_freq_vs_score.R``   
    - Input files:
      - [./result/NTD_DMS_scores_by_resi.tsv](./result/NTD_DMS_scores_by_resi.tsv)
      - [./result/residue_freq.csv](./result/residue_freq.csv)
    - Output file:
      - [./graph/dms_expression_vs_seq_conservation.png](./graph/dms_expression_vs_seq_conservation.png)
 
10. Plot mutational tolerability vs sequence conservation for individual residues   
``Rscript script/RSA_vs_alignment_frequency.R``   
    - Input files:
      - [./result/residue_freq.csv](./result/residue_freq.csv)
      - [./result/NTD_RSA.tsv](./result/NTD_RSA.tsv)
    - Output files:
      - [./graph/RSA_vs_seq_con.png](./graph/RSA_vs_seq_con.png)

11. Visualizing the mutational tolerability on the S protein structure
``Pymol script/plot_Bfactor_as_exp.pml`` 
    - Input files:
      - [./PDB/7b62_exp.pdb](./PDB/7b62_exp.pdb)
      - [./PDB/6zge.pdb](./PDB/6zge.pdb)
    - Output files:
      - [./PDB/overlay.png](./PDB/overlay.png)
      - [./PDB/Overlay_closeup_0.png](./PDB/Overlay_closeup_0.png)
      - [./PDB/Overlay_closeup_180.png](./PDB/Overlay_closeup_180.png)
      - [./PDB/Overlay_closeup_topdown.png](./PDB/Overlay_closeup_topdown.png)
      - [./PDB/circulating_mutation_1.png](./PDB/circulating_mutation_1.png)
      - [./PDB/circulating_mutation_2.png](./PDB/circulating_mutation_2.png)
      - [./PDB/circulating_mutation_3.png](./PDB/circulating_mutation_3.png)
      - [./PDB/circulating_indel_1.png](./PDB/circulating_indel_1.png)
      - [./PDB/circulating_indel_2.png](./PDB/circulating_indel_1.png)
      - [./PDB/circulating_indel_3.png](./PDB/circulating_indel_1.png)

12. Analysis of the ciculating NTD mutations/indels among 17 major variants
``Rscript script/NTD_circulating_mutation_vs_other_residues.R``  
    - Input files:
      - [./result/NTD_DMS_scores.tsv](./result/NTD_DMS_scores.tsv)
      - [./result/NTD_DMS_scores_by_resi.tsv](./result/NTD_DMS_scores_by_resi.tsv)
    - Output files:
      - [./graph/ntd_circulating_mutations_vs_other_mutations.png](./graph/ntd_circulating_mutations_vs_other_mutations.png)
      - [./graph/ntd_indel_sites_vs_other_sites.png](./graph/ntd_indel_sites_vs_other_sites.png)

13. Mutational tolerability of selected regions compared to other residues
``Rscript script/NTD_loop_vs_other_residues.R``  
    - Input files:
      - [./result/NTD_DMS_scores_by_resi.tsv](./result/NTD_DMS_scores_by_resi.tsv)
    - Output files:
      - [./graph/NTD_loop_vs_other_residues.png](./graph/NTD_loop_vs_other_residues.png)
``Rscript script/Hot_spots_vs_other_residues.R``  
     - Input files:
      - [./result/NTD_DMS_scores_by_resi.tsv](./result/NTD_DMS_scores_by_resi.tsv)
      - Output files:
      - [./graph/hot_spots_vs_other_residues.png](./graph/hot_spots_vs_other_residues.png)
