## SARS2_NTD_DMS

### Dependencies
* python=3.9
* [biopython](https://github.com/biopython/biopython)

### Input files
* [./Fasta/SARS2-NTD.fa](./Fasta/SARS2-NTD.fa): NTD sequences with 21 nt upstream (5' flank)

### Primer design for DMS library construction
1. Generating foward (NNK + internal barcode) and reverse primers (constant)   
```python3 script/lib_primer_design.py```

2. Checking internal barcode   
```python3 script/check_barcode.py```
