## Studying the impact of NTD mutations on SARS-CoV-2 spike expression using deep mutational scanning

### Dependencies
* [Python](https://www.python.org/) (version 3.9)
* [Biopython](https://github.com/biopython/biopython)
* [R](https://www.r-project.org/) (version 4.1)
* [PEAR](https://github.com/tseemann/PEAR)

### Input files
* [./Fasta/SARS2-NTD.fa](./Fasta/SARS2-NTD.fa): NTD sequences with 21 nt upstream (5' flank)
* [./Fasta/Amplicon.fa](./Fasta/Amplicon.fa): Amplicon sequences for NGS (including recovery primer regions)
* [./Fasta/NTD_ref.fa](./Fasta/NTD_ref.fa): Reference (i.e. wild type) amino acid sequences (primer regions not included)

### Primer design for DMS library construction
1. Generating foward (NNK + internal barcode) and reverse primers (constant)   
```python3 script/lib_primer_design.py```
    - Input file:
      - [./Fasta/SARS2-NTD.fa](./Fasta/SARS2-NTD.fa)
    - Output files:
      - [./primer/SARS2-NTD_lib_Fprimer_bc.fa](./primer/SARS2-NTD_lib_Fprimer_bc.fa)
      - [./primer/SARS2-NTD_lib_Rprimer_bc.fa](./primer/SARS2-NTD_lib_Rprimer_bc.fa)

2. Generating barcode file   
```python3 script/check_barcode.py```
    - Input files:
      - [./primer/SARS2-NTD_lib_Fprimer_bc.fa](./primer/SARS2-NTD_lib_Fprimer_bc.fa)
      - [./Fasta/Amplicon.fa](./Fasta/Amplicon.fa)
    - Output file:
      - [./data/barcodes.tsv](./data/barcodes.tsv)

### Calculating experssion score from DMS data
1. Merge overlapping paired-end reads using [PEAR](https://github.com/tseemann/PEAR)   
```pear -f [FASTQ FILE FOR FORWARD READ] -r [FASTQ FILE FOR FORWARD READ] -o [OUTPUT FASTQ FILE]```   
    - Output files should be placed in the fastq_merged/ folder

2. Counting variants based on nucleotide sequences   
```python3 script/NTD_fastq2count.py```   
    - Input files:
      - Merged read files in fastq_merged/ folder
    - Output files:
      - result/NTD_DMS_count_nuc_A.tsv
      - result/NTD_DMS_count_nuc_B.tsv

3. Convert nucleotide sequences to amino acid mutations   
```python3 script/NTD_count_nuc2aa.py```   
    - Input files:
      - [./data/barcodes.tsv](./data/barcodes.tsv)
      - result/NTD_DMS_count_nuc_A.tsv
      - result/NTD_DMS_count_nuc_B.tsv
    - Output files:
      - result/NTD_DMS_count_aa_A.tsv
      - result/NTD_DMS_count_aa_B.tsv

4. Compute expression score   
```python3 script/NTD_count2score.py```   
    - Input files:
      - result/NTD_DMS_count_aa_A.tsv
      - result/NTD_DMS_count_aa_B.tsv
    - Output file:
      - [./result/NTD_DMS_scores.tsv](./result/NTD_DMS_scores.tsv)
