---
title: "Accessing the NCBI Entrez database using Biopython"
date: 2022-05-14
lastmod: 2022-05-14
draft: false
garden_tags: ["python", "bioinformatics", "datasets"]
summary: "This tutorial uses the NCBI API to interface with Entrez."
status: "seeding"
---

#### Background

This tutorial uses the NCBI API to interface with Entrez. To get started, the necessary libraries need to be imported and an email needs to be provided (so NCBI can contact you about your query if needed).

```python
email_address = "gurkamal.dev@gmail.com" 
```

The goal will be to find the chloroquine resistance transporter (CRT) gene in the parasite Plasmodium flaciparum in the nucleotide database.

#### Loading Libraries

```python
from Bio import Entrez, Medline, SeqIO
```

#### What is a handle?

A handle is essentially a “wrapper” around text information.

Handles provide two benefits over plain text information:

They provide a standard way to deal with information stored in different ways. The text information can be in a file, or in a string stored in memory, or the output from a command line program, or at some remote website, but the handle provides a common way of dealing with information in all of these formats.

They allow text information to be read incrementally, instead of all at once. This is really important when dealing with huge text files which would use up all of the memory if you had to load them all.

Handles can deal with text information that is being read (reading from a file) or written (writing information to a file). In the case of a “read” handle, commonly used functions are ```read()```, which reads the entire text information from the handle, and ```readline()```, which reads information one line at a time. For “write” handles, the function ```write()``` is regularly used.

#### Retrieving information

To see the available databases:

```python
#This gives you the list of available databases
handle = Entrez.einfo()

#Read and store the Entrez query record returned
rec = Entrez.read(handle)
print(rec)


{'DbList': ['pubmed', 'protein', 'nuccore', 'ipg', 'nucleotide', 'structure', 'genome', 'annotinfo', 'assembly', 'bioproject', 'biosample', 'blastdbinfo', 'books', 'cdd', 'clinvar', 'gap', 'gapplus', 'grasp', 'dbvar', 'gene', 'gds', 'geoprofiles', 'homologene', 'medgen', 'mesh', 'ncbisearch', 'nlmcatalog', 'omim', 'orgtrack', 'pmc', 'popset', 'proteinclusters', 'pcassay', 'protfam', 'pccompound', 'pcsubstance', 'seqannot', 'snp', 'sra', 'taxonomy', 'biocollections', 'gtr']}
```

#### Searching for a specific gene

```python
#search the nucleotide database for our gene and organism
handle = Entrez.esearch(db="nucleotide", term='CRT[Gene Name] AND "Plasmodium falciparum"[Organism]')

#read the result that is returned
rec_list = Entrez.read(handle)


{'Count': '2022', 'RetMax': '20', 'RetStart': '0', 'IdList': ['2196471109', '2196471107', '2196471105', '2196471103', '2196471101', '2196471099', '2196471097', '2196471095', '2196471093', '2196471091', '2196471089', '2196471087', '2196471085', '2196471083', '2196471081', '2196471079', '2196471077', '2196471075', '2196471073', '2196471071'], 'TranslationSet': [{'From': '"Plasmodium falciparum"[Organism]', 'To': '"Plasmodium falciparum"[Organism]'}], 'TranslationStack': [{'Term': 'CRT[Gene Name]', 'Field': 'Gene Name', 'Count': '4778', 'Explode': 'N'}, {'Term': '"Plasmodium falciparum"[Organism]', 'Field': 'Organism', 'Count': '258609', 'Explode': 'Y'}, 'AND'], 'QueryTranslation': 'CRT[Gene Name] AND "Plasmodium falciparum"[Organism]'}
```

#### Returning all records

The standard search will limit the number of record references to 20, so if you have more, you may want to repeat the query with an increased maximum limit. In this case, we will actually override the default limit with retmax. The Entrez system provides quite a few sophisticated ways to retrieve large number of results.

Be careful with this technique, because you will retrieve a large amount of complete records, and some of them will have fairly large sequences inside.

```python
if rec_list['RetMax'] < rec_list['Count']:
    handle = Entrez.esearch(db="nucleotide", term='CRT[Gene Name] AND "Plasmodium falciparum"[Organism]',
                            retmax=rec_list['Count'])
    rec_list = Entrez.read(handle)


```
#### Downloading nucleotide sequences

Now that we have the IDs of all of the records, you still need to retrieve the records properly.

This will retrieve a list of records in the GenBank format (including sequences and metadata)

```python
# query will download all matching nucleotide sequences from GenBank
id_list = rec_list['IdList']

handle_2 = Entrez.efetch(db='nucleotide', id=id_list, rettype='gb', retmax=rec_list['Count'])

```

The result of ```SeqIO.parse``` is an iterator and is converted to a list. The advantage of doing this is that we can use the result as many times as we want(for example, iterate many times over), without repeating the query on the server.This saves time, bandwidth, and server usage if you plan to iterate many timesover. 

The disadvantage is that it will allocate memory for all records. This will not work for very large datasets

```python
recs = list(SeqIO.parse(handle_2, 'gb'))
```

#### Reading a record

The ```rec```  variable now contains the record of interest. The ```rec.description``` will contain its human-readable description.

```python
for rec in recs:
    if rec.name == 'KM288867':
        break
print(rec.name)
print(rec.description)


KM288867
Plasmodium falciparum clone PF3D7_0709000 chloroquine resistance transporter (CRT) gene, complete cds
```
#### Extracting sequences features 
Extract sequence features which contain information such as gene products and exon positions on the sequence:

```python
#If the feature.type is gene, print its name
for feature in rec.features:
    if feature.type == 'gene':
        # gene name will be in the feature.qualifiers dictionary.
        print(feature.qualifiers['gene'])
        # print the start, end, and originating strand of the exon
    elif feature.type == 'exon':
        loc = feature.location
        print('Exon', loc.start, loc.end, loc.strand)
    else:
        print('not processed:\n%s' % feature)


not processed:
type: source
location: [0:10000](+)
qualifiers:
    Key: clone, Value: ['PF3D7_0709000']
    Key: db_xref, Value: ['taxon:5833']
    Key: mol_type, Value: ['genomic DNA']
    Key: organism, Value: ['Plasmodium falciparum']

['CRT']
not processed:
type: mRNA
location: join{[2751:3543](+), [3720:3989](+), [4168:4341](+), [4513:4646](+), [4799:4871](+), [4994:5070](+), [5166:5249](+), [5376:5427](+), [5564:5621](+), [5769:5862](+), [6055:6100](+), [6247:6302](+), [6471:7598](+)}
qualifiers:
    Key: gene, Value: ['CRT']
    Key: product, Value: ['chloroquine resistance transporter']

not processed:
type: 5'UTR
location: [2751:3452](+)
qualifiers:
    Key: gene, Value: ['CRT']

not processed:
type: primer_bind
...
type: primer_bind
location: [7833:7856](-)
qualifiers:
```
#### Record annotations

We will now look at the annotations on the record, which are mostly metadata that is not related to the sequence position. Note that some values are not strings; they can be numbers or even lists (for example, the taxonomy annotation is a list)

```python
for name, value in rec.annotations.items():
    print('%s=%s' % (name, value))



molecule_type=DNA
topology=linear
data_file_division=INV
date=12-NOV-2014
accessions=['KM288867']
sequence_version=1
keywords=['']
source=Plasmodium falciparum (malaria parasite P. falciparum)
organism=Plasmodium falciparum
taxonomy=['Eukaryota', 'Sar', 'Alveolata', 'Apicomplexa', 'Aconoidasida', 'Haemosporida', 'Plasmodiidae', 'Plasmodium', 'Plasmodium (Laverania)']
references=[Reference(title='Versatile control of Plasmodium falciparum gene expression with an inducible protein-RNA interaction', ...), Reference(title='Direct Submission', ...)]
```

#### Accessing the sequence data

Last but not least, you can access the fundamental piece of information, the sequence

```python
print(rec.seq[0:100])

ATATGTAAAACCAAAATAAATTAAACAGAATTTATTTTT
AAAAGATTTATTTGTAACAATATTACCATGATGATTTAT
TAAAGTAAAATCACCACCTATT
```
