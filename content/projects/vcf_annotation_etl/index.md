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

Human genome assembly GRCh38 (hg38) was used as the reference genome. A VCF file from dbSNP and a GFF3 file for annotation was obtained from the Human Genome Resources at NCBI[(https://www.ncbi.nlm.nih.gov/genome/guide/human/)].

### Tools used

The pipeline is built using Pandas, Pysam's tabix interface, and SQLite3. Tabix is used to index and access the VCF file as it's ~26gb compressed and extracting the file would be resource consuming.

I decided to use Streamlit to build a web app for this pipeline so that it can be shared and modifed to print logs, stream data from AWS S3, or build it out into a multipage app.

Code refactoring was done using Sourcery and Autopep8 was used to ensure PEP8 formatting

### Project structure

The pipeline is built using three modules that are  named: extract, transform, load
##### Extract:
The extraction phase calls the extract functions to parse the data of the VCF and GFF3 files into pandas dataframes.

##### Transform:

The transformation phase formats the VCF and then takes the gff dataframe as input into the pandas pipe, which is a group of functions that take the output of the preceding function as their input. The annotation step merges the two processed dataframes on the Gene column

##### Load:
The load phase simply creates an SQLite3 database and the pandas 'to_sql' function handles the upload to the SQL table. Lastly, a list of genes is parameterized and sent as a uniprot query to return protein sequences.


![etl](/etl.png)

# References
1. Bornstein P, Sage EH. Matricellular proteins: extracellular modulators of cell function. Curr Opin Cell Biol 2002; 14:608–616. doi: 10.1016/S0955-0674(02)00361-7
2. Vannahme C, Gosling S, Paulsson M, Maurer P, Hartmann U. Characterization of SMOC-2, a modular extracellular calcium-binding protein. Biochem J 2003; 373:805–814. doi: 10.1042/bj20030532.