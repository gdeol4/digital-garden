---
title: "Genomics ETL pipeline"
date: 2022-07-04
draft: false
project_tags: ["Bioinformatics", "Python", "Pandas", ETL]
status: "evergreen"
weight: 3
summary: "Building an ETL pipeline with Pandas to annotate mutations in the human genome"
---



# Introduction

This project was a take home technical assessment from an AI based Biotech company. The goal of the project is to build an ETL pipeline with python that takes in a VCF file annotates the mutations, creates an SQL database, and loads the transformed data into a table.

### Sources of data

Human genome assembly GRCh38 (hg38) was used as the reference genome. A VCF file from dbSNP and a GFF3 file for annotation was obtained from the Human Genome Resources at NCBI[(https://www.ncbi.nlm.nih.gov/genome/guide/human/)]

![plot pca](/deseq2pca.png)

# References
1. Bornstein P, Sage EH. Matricellular proteins: extracellular modulators of cell function. Curr Opin Cell Biol 2002; 14:608–616. doi: 10.1016/S0955-0674(02)00361-7
2. Vannahme C, Gosling S, Paulsson M, Maurer P, Hartmann U. Characterization of SMOC-2, a modular extracellular calcium-binding protein. Biochem J 2003; 373:805–814. doi: 10.1042/bj20030532.