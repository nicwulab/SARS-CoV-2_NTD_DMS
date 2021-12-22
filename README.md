## Studying the impact of NTD mutations on SARS-CoV-2 spike expression using deep mutational scanning

### Dependencies
* [Python](https://www.python.org/) (version 3.9)
* [Biopython](https://github.com/biopython/biopython)
* [R](https://www.r-project.org/) (version 4.1)
* [PEAR](https://github.com/tseemann/PEAR)

### Input files
* [./Fasta/SARS2-NTD.fa](./Fasta/SARS2-NTD.fa): NTD sequences with 21 nt upstream (5' flank)
* [./Fasta/Amplicon.fa](./Fasta/Amplicon.fa): Amplicons for sequencing (including recovery primer regions)

### Primer design for DMS library construction
1. Generating foward (NNK + internal barcode) and reverse primers (constant)   
~~~python3 script/lib_primer_design.py~~~
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
