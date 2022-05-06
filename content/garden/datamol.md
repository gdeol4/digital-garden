---
title: "Preprocessing molecular data with Datamol"
date: 2021-10-22
lastmod: 2021-10-22
draft: false
garden_tags: ["python", "cheminformatics", "drug discovery", "data"]
summary: "Sanitizing and manipulating labeled blood-brain barrier permeable molecules."
status: "growing"
---

SMILES (Simplified Molecular Input Line Entry System) is a standard notation representing the molecular structure of a compound as a string representation that can be understood by a computer. The SMILES notation consists of a handful of rules which allow for converting the string to an image or graph. SMILES can then be easily used for generating further representations to train machine learning models with.

```python
import datamol as dm
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
```
```python
BBBP_df = pd.read_csv("data/BBBP.csv")
BBBP_df.head()
```

---

| num | name | p_np | smiles |
|:---------|:---------|:---------|:---------|
| 1 | Propanolol | 1 | [Cl].CC(C)NCC(O)COc1cccc2ccccc12 | 
| 2 | Terbutylchlorambucil | 1 | C(=O)(OC(C)(C)C)CCCc1ccc(cc1)N(CCCl)CCCl |
| 3 | 40730 | 1 | c12c3c(N4CCN(C)CC4)c(F)cc1c(c(C(O)=O)cn2C(C)CO... |
| 4 | 24 | 1 | C1CCN(CC1)Cc1cccc(c1)OCCCNC(=O)C | 
| 5 | cloxacillin | 1 | Cc1onc(c2ccccc2Cl)c1C(=O)N[C@H]3[C@H]4SC(C)(C)... | 

The dataframe shows 4 named columns, including the "num" of the molecule, the name, a binary label for blood brain barrier permeability status "p_np", and the SMILES string.

```python
# The name and number can be dropped
BBBP_df = BBBP_df.drop(["num", "name"], axis=1)

# Checking the data for null values
BBBP_df["smiles"].isnull().values.any()

# Renaming the binary label to "BBB+/BBB-" for clarity
BBBP_df.columns = ['BBB+/BBB-', 'SMILES']
```
```python
BBBP_df
```

---

| BBB+/BBB- | SMILES |
|:---------|:---------|
| 0 | 1 | [Cl].CC(C)NCC(O)COc1cccc2ccccc12 |
| 1 | 1 | C(=O)(OC(C)(C)C)CCCc1ccc(cc1)N(CCCl)CCCl |
| 2048 | 1 | C1=C(OC)C(=CC2=C1C(=[N+](C(=C2CC)C)[NH-])C3=CC... |
| 2049 | 1 | [N+](=NCC(=O)N[C@@H]([C@H](O)C1=CC=C([N+]([O-]... |


Mols and smiles need to be sanitized as it will leave us with SMILES that are complete nonesense, for example, errors resulting from kekulization.

![kekul](/kekul.jpg)

RDkit generates the alternate position of double bonds, and then (in a second step they call "aromatization") labels the ring as aromatic. In panel (2), there are three possible Lewis structures contributing to the actual structure (i.e. there is resonance), so the software would have to generate all three to be able to search for identical structures.[^1]

Below is a function using datamol to preprocess the dataset, including steps to generate mol objects, SELFIES, inchi, and inchikeys for each molecule. The function also standardizes mols and SMILES, drops NA values, and returns a dataframe.

```python
def preprocess_smiles(df):
    df["mol"] = [dm.to_mol(x) for x in df['SMILES']] # generating mols from SMILES
    df["mol"] = [dm.fix_mol(x) for x in df['mol']] # Fixing mols

    df = df.dropna() # dropping NA values

    df["mol"] = [dm.sanitize_mol(x, sanifix=True, charge_neutral=False) for x in df['mol']] # sanitize mol objects
    df["mol"] = [dm.standardize_mol(x, disconnect_metals=False, normalize=True, reionize=True, uncharge=False, stereo=True) for x in df['mol']] # standardize mol objects

    df["standard_smiles"] = [dm.standardize_smiles(x) for x in df['SMILES']] # standardize SMILES
    df["selfies"] = [dm.to_selfies(x) for x in df['mol']] # generate SELFIES
    df["inchi"] = [dm.to_inchi(x) for x in df['mol']] # Generating InChi
    df["inchikey"] = [dm.to_inchikey(x) for x in df['mol']] # Generating InChIKey

    return df
```

Running the function and taking a look at the outputs

```python
data_clean = preprocess_smiles(BBBP_df)
```

```python
data_clean.shape
```

```python
(2039, 7)
```

The data contains a 3:1 ratio of positive to negative labels, which creates a bias towards molecules with blood brain permeability properties. This may need to be addressed when training models. The next steps are to save the cleaned data for further analysis.

```python
counts = data_clean['BBB+/BBB-'].value_counts().to_dict()
print(counts)
```

```python
{1: 1560, 0: 479}
```

```python
data_clean.to_csv('./data/MoleculeNet.csv', index=False)
```

# References
[^1]: Urbaczek, Sascha. A consistent cheminformatics framework for automated virtual screening. Ph.D. Thesis, Universität Hamburg, August 2014. URL: http://ediss.sub.uni-hamburg.de/volltexte/2015/7349/; URN: urn:nbn:de:gbv:18-73491; PDF via Semantic Scholar