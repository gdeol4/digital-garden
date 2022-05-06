---
title: "A python function to perform restriction enzyme digest"
date: 2021-10-18
lastmod: 2021-10-18
draft: false
garden_tags: ["python", "bioinformatics"]
summary: "Implementing restriction digestion using regular expressions, lists, and dictionaries."
status: "evergreen"
---

Restriction enzymes are proteins produced by bacteria that cleave DNA at specific sites along the molecule. The enzyme functions on a specific, short nucleotide sequence and cuts the DNA only at that specific site, which is known as restriction site or target sequence. In the bacterial cell, restriction enzymes cleave foreign DNA, thus eliminating infecting organisms. The activity of a restriction enzyme can be defined by its recognition site on the DNA sequence and the position relative to the recognition site, at which it cuts the DNA.

```python
# create enzyme dictionary
restrictionEnzymes = {}

# add "bamH1" and "sma1" enzymes, their target sequence and
# position releative to the recognition site
restrictionEnzymes['bamH1'] = ['ggatcc',0]
restrictionEnzymes['sma1'] = ['cccggg',2]

# a function to calculate the molecular weight of dna sequences
def oligoMolecularWeight(sequence):

    # create a dictionairy of DNA basepair molecular weights
    dnaMolecularWeight = {'a':313.2,'c':289.2,'t':304.2,'g':329.2}

    # initialize molecular weight
    molecularWeight = 0.0

    # iterate through DNA sequnce and update weight of sequence
    for base in sequence:
        molecularWeight += dnaMolecularWeight[base]
    return molecularWeight

# the primary function for restriction digest
def digest(sequence, enzyme):
    # set target sequence
    target = restrictionEnzymes[enzyme][0]

    # enzyme cut position relative to recognition site
    cutPosition = restrictionEnzymes[enzyme][1]

    # a list to collect DNA fragments
    fragments = []

    # counter for the position of the last restriction site; 
    # beginning of sequence
    found = 0

    # a variable to store the position of the last cut;
    # end of sequence
    lastCut = found

    # variable to set where to search for the next site from
    searchFrom = lastCut

    while found != -1:
        found = sequence.find(target, searchFrom)
        if found != -1:
            fragment = sequence[lastCut:found+cutPosition]
            mwt = oligoMolecularWeight(fragment)
            fragments.append((fragment,mwt))
        else:
            fragment = sequence[lastCut:]
            mwt = oligoMolecularWeight(fragment)
            fragments.append((fragment,mwt))
        lastCut = found + cutPosition
        searchFrom = lastCut + 1
    
    return fragments
```
Running the function on a test sequence results in the following:

```python
digestSequence = "gcgatgctaggatccgcgatcgcgtacgatcgtacgcggtacggacggatccttctc"
```

```python
digested_dna = digest(digestSequence,'bamH1')
```
```python
print(digested_dna)
```
```python
[('gcgatgcta', 2800.7999999999997), 
('ggatccgcgatcgcgtacgatcgtacgcggtacggac', 11478.400000000005), 
('ggatccttctc', 3345.1999999999994)]
```