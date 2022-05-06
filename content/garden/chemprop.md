---
title: "Scaffold splitting and model training using chemprop"
date: 2021-10-22
lastmod: 2021-10-22
draft: false
garden_tags: ["python", "cheminformatics", "Drug Discovery", "Deep Learning"]
summary: "Chemically aware data splitting  using scaffold splitting. Data is then trained using a directed message passing neural network (D-MPNN)."
status: "seeding"
---

# Chemprop

I found the chemprop library in the paper A Deep Learning Approach to Antibiotic Discovery30102-1"), it uses the directed message passing neural network (D-MPNN) described in Analyzing Learned Molecular Representations for Property Prediction.

The drug discovery workflow described in the paper uses a training dataset of 2335 molecules with binary e.coli growth inhibition labels to train a classification model.

The D-MPNN architecture as described in the paper:

>The D-MPNN architecture translates the graph representation of a molecule into a continuous vector via a directed bond-based message-passing approach. This builds a molecular representation by iteratively aggregating the features of individual atoms and bonds. The model operates by passing “messages” along bonds that encode information about neighboring atoms and bonds. By applying this message passing operation multiple times, the model constructs higher-level bond messages that contain information about larger chemical substructures. The highest-level bond messages are then combined into a single continuous vector representing the entire molecule.

The researchers generated molecular features using RDKit to tackle overfitting, increased algorithm robustness by using an ensemble of classifiers, and optimized hyperparameters with bayesian optimization. The model achieved a ROC-AUC of 0.896 on test data.

This notebook will apply the neural network to the dataset of blood-brain barrier molecules gathered in the previous notebooks.


# Scaffold splitting

A distributional shift is a change in the data distribution between an algorithm's training dataset, and a dataset it encounters when deployed. Distributional shifts are common in the field of molecular property prediction where chemical space is huge and different areas of it exhibit high structural heterogeneity. The shifts tend to be large and difficult to handle for machine learning models.

Scaffold splitting is one solution to this distributional shift. It's first described in "The properties of known drugs. 1. Molecular frameworks:

>A molecular scaffold reduces the chemical structure of a compound to its core components, essentially by removing all side chains and only keeping ring systems and parts that link together ring systems. An additional option for making molecular scaffolds even more general is to “forget” the identities of the bonds and atoms by replacing all atoms with carbons and all bonds with single bonds.

A Murcko scaffold essentially extracts the molecular backbone by extracting the ring structures and the linkers that connect them. These scaffolds collapse molecules into bins based on their generated Murcko scaffold, reducing the number of highly similar structures. When splitting the data by scaffold, molecules sharing a scaffold are in the same split, meaning no similar molecules will be found across different splits. Any bins larger than half of the desired test set size are placed into the training set, to guarantee the scaffold diversity of the validation and test sets. All remaining bins are placed randomly into the training, validation, and test sets until each set has reached its desired size. Clustering molecules by scaffold ensures that the model will be trained on structurally diverse folds when cross-validating.

Compared to a random split, a scaffold split is a more challenging and realistic evaluation setting as it more closely approximates the split present in real-world property prediction data.

A molecule compared to its scaffold below:

### Molecular compound and its generic Bemis-Murcko scaffold

```python
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

# define compound via its SMILES string
smiles = "CN1CCCCC1CCN2C3=CC=CC=C3SC4=C2C=C(C=C4)SC"
# convert SMILES string to RDKit mol object 
mol = Chem.MolFromSmiles(smiles)
# create RDKit mol object corresponding to Bemis-Murcko scaffold of original compound
mol_scaffold = MurckoScaffold.GetScaffoldForMol(mol)
# make the scaffold generic by replacing all atoms with carbons and all bonds with single bonds
mol_scaffold_generic = MurckoScaffold.MakeScaffoldGeneric(mol_scaffold)
# convert the generic scaffold mol object back to a SMILES string format
smiles_scaffold_generic = Chem.CanonSmiles(Chem.MolToSmiles(mol_scaffold_generic))
# display compound and its generic Bemis-Murcko scaffold
display(mol)
print(smiles)
display(mol_scaffold_generic)
print(smiles_scaffold_generic)
```

![prop1](/prop1.png)

```python
CN1CCCCC1CCN2C3=CC=CC=C3SC4=C2C=C(C=C4)SC
```

![prop2](/prop2.png)

```python
C1CCC(CCC2C3CCCCC3CC3CCCCC32)CC1
```

# Generating generic Bemis-Murcko scaffolds

Below I've written code to generate the scaffolds and a few other related molecular representations. But binning and writing the code to split them myself would take too much time so I use a built-in feature of chemprop to split the data before model training.

```python
# import packages

import chemprop
import pandas as pd
import datamol as dm
pd.options.mode.chained_assignment = None  # default='warn'


# loading the datasets
MolNet = pd.read_csv("data/MoleculeNet.csv")
B3DB = pd.read_csv("data/B3DB.csv")
b3_molecules = pd.read_csv("data/b3_molecules.csv")
```

## scaffold_split function

This function generates and sanitizes mols, generates a scaffold from the mol object, generalizes the scaffold, and lastly converts the scaffold to a SMILES string, returning a dataframe.

```python
def scaffold_split(df):
    df["mol"] = [Chem.MolFromSmiles(x) for x in df["standard_smiles"]] # generating mols from the standard_smiles column
    df["mol"] = [dm.sanitize_mol(x, sanifix=True, charge_neutral=True) for x in df['mol']] # sanitize mol objects
    df = df.dropna() # dropping NA values
    df["scaffold"] = [MurckoScaffold.GetScaffoldForMol(x) for x in df["mol"]] # generating scaffolds from mol object
    df["mol_scaffold_generic"] = [MurckoScaffold.MakeScaffoldGeneric(x) for x in df["scaffold"]] # generalizing scaffolds
    # convert the generic scaffold mol object back to a SMILES string format
    df["smiles_scaffold_generic"] = [Chem.CanonSmiles(Chem.MolToSmiles(x)) for x in df["mol_scaffold_generic"]]
    
    return df
```

```python
# the data prior to processing
MolNet.head(1)
```
---

| BBB+/BBB- | SMILES | mol | standard_smiles | selfies | inchi | inchikey |
|:---------|:---------|:---------|:---------|:---------|:---------|:---------|
| 1 | [Cl].CC(C)NCC(O)COc1cccc2ccccc12 | <img data-content="rdkit/molecule" src...> | CC(C)NCC(O)COc1cccc2ccccc12.[Cl-] | [C][C][Branch1][C][C][N][C][C][Branch1][C][O][... | InChI=1S/C16H21NO2.ClH/c1-12(2)17-10-14(18)11-... | ZMRUPTIKESYGQW-UHFFFAOYSA-M |

```python
# results of the scaffold generating function
data_split = scaffold_split(b3_molecules)
data_split.head(1)
```
The three additionally generated columns below:

---

| scaffold | mol_scaffold_generic | smiles_scaffold_generic |
|:---------|:---------|:---------|
| img data>-content="rdkit/molecule" src... | img data-content="rdkit/molecule" src... | C1CCC2CCCCC2C1 |

## chemprop_prep function

This function serves to prepare the data to be used by the chemprop library. Specifically, the input data must contain a SMILES string in the first column and the target value - a binary value in this experiment. The resulting dataframe is saved as a CSV file to later be used by chemprop.

```python
# A function to process the data for use with the chemprop library
def chemprop_prep(df, filename):
    df = df.drop(["mol", "SMILES", "selfies", "inchi", "inchikey"], axis=1) # drop all columns except the smiles and target
    df["smiles"] = df["standard_smiles"] # use standard smiles inplace of smiles
    df = df.drop(["standard_smiles"], axis=1) # drop this column now
    df = df[["smiles", "BBB+/BBB-"]] # reorder the columns with smiles first and target second
    df.to_csv('./data/' + filename + '.csv', index=False) # save the file

    return df
```

```python
# Processing the three different datasets
molnet_chemprop = chemprop_prep(MolNet, 'molnet_chemprop')
B3DB_chemprop = chemprop_prep(B3DB, 'B3DB_chemprep')
b3_mol_chemprop = chemprop_prep(b3_molecules, 'b3_mol_chemprop')
```

```python
# results of processing
molnet_chemprop.head(1)
```

| smiles | BBB+/BBB- |
|:---------|:---------|
| CC(C)NCC(O)COc1cccc2ccccc12.[Cl-] | 1 |

# Model training with chemprop

I ran three experiments with slight differences to test out the predictive ability on the blood-brain barrier molecules:

--- 

1. The split type is set to scaffold balanced; the number of folds is 10, class balancing is set to true, and 200 features are generated with RDKits 2d descriptor generator. The features are scaled when calculated and thus no further scaling is specified.
2. This run sets the number of training epochs within each fold at 100. The rest of the settings mentioned above are kept the same.
3. The last run has the number of folds increased from 10 to 30 while epochs per fold are set to 20.

---

These models are trained using the MoleculeNet data (2035 molecules), further tests will use the combined dataset with ~8100 molecules.

## Model training 1

```python
arguments = [
    '--data_path', './data/molnet_chemprop.csv',
    '--dataset_type', 'classification',
    '--save_dir', './data/chemprop_checkpoints/',
    '--split_type', 'scaffold_balanced',
 #   '--separate_val_path', './data/chemprop_B3DB.csv',
    '--num_folds', '10',
    '--class_balance',
    '--features_generator', 'rdkit_2d_normalized',
    '--no_features_scaling',
    '--quiet'
    

]

args = chemprop.args.TrainArgs().parse_args(arguments)
mean_score, std_score = chemprop.train.cross_validate(args=args, train_func=chemprop.train.run_training)
```

## Model training 2

```python
arguments = [
    '--data_path', './data/molnet_chemprop.csv',
    '--dataset_type', 'classification',
    '--save_dir', './data/chemprop_checkpoints/',
    '--split_type', 'scaffold_balanced',
 #   '--separate_val_path', './data/chemprop_B3DB.csv',
    '--num_folds', '10',
    '--class_balance',
    '--features_generator', 'rdkit_2d_normalized',
    '--no_features_scaling',
    '--quiet',
    '--epochs', '100'
    

]

args = chemprop.args.TrainArgs().parse_args(arguments)
mean_score, std_score = chemprop.train.cross_validate(args=args, train_func=chemprop.train.run_training)
```

## Model training 3


```python
arguments = [
    '--data_path', './data/molnet_chemprop.csv',
    '--dataset_type', 'classification',
    '--save_dir', './data/chemprop_checkpoints/',
    '--split_type', 'scaffold_balanced',
 #   '--separate_val_path', './data/chemprop_B3DB.csv',
    '--num_folds', '30',
    '--class_balance',
    '--features_generator', 'rdkit_2d_normalized',
    '--no_features_scaling',
    '--quiet',
    '--epochs', '20',
    

]

args = chemprop.args.TrainArgs().parse_args(arguments)
mean_score, std_score = chemprop.train.cross_validate(args=args, train_func=chemprop.train.run_training)
```

# Results

### Experiment 1:

10-fold cross validation

```python
Seed 0 ==> test auc = 0.928767
Seed 1 ==> test auc = 0.903879
Seed 2 ==> test auc = 0.955724
Seed 3 ==> test auc = 0.862450
Seed 4 ==> test auc = 0.890246
Seed 5 ==> test auc = 0.892326
Seed 6 ==> test auc = 0.886606
Seed 7 ==> test auc = 0.951849
Seed 8 ==> test auc = 0.906486
Seed 9 ==> test auc = 0.937866
```

Overall test auc = 0.911620 +/- 0.029164 Elapsed time = 0:20:54

This run was the fastest and used the default number of epochs and fold. The overall AUC for the test data was 0.911620 +/- 0.029164.

### Experiment 2:

10-fold cross validation

```python
Seed 0 ==> test auc = 0.929177
Seed 1 ==> test auc = 0.852848
Seed 2 ==> test auc = 0.944112
Seed 3 ==> test auc = 0.862913
Seed 4 ==> test auc = 0.890246
Seed 5 ==> test auc = 0.893487
Seed 6 ==> test auc = 0.881697
Seed 7 ==> test auc = 0.941933
Seed 8 ==> test auc = 0.920437
Seed 9 ==> test auc = 0.917494
```


Overall test auc = 0.903435 +/- 0.030385 Elapsed time = 16:49:39

This run took the most time, nearing almost 17 hours with 100 epochs for each fold. The overall AUC for the test data was 0.903435 +/- 0.030385 which underperformed compared to the model in the first experiment that ran using default parameters.

### Experiment 3:

30-fold cross validation

```python
Seed 0 ==> test auc = 0.928630
Seed 1 ==> test auc = 0.907636
Seed 2 ==> test auc = 0.943526
Seed 3 ==> test auc = 0.882607
Seed 4 ==> test auc = 0.890246
Seed 5 ==> test auc = 0.892559
Seed 6 ==> test auc = 0.870606
Seed 7 ==> test auc = 0.947479
Seed 8 ==> test auc = 0.903720
Seed 9 ==> test auc = 0.930991
Seed 10 ==> test auc = 0.795291
Seed 11 ==> test auc = 0.894284
Seed 12 ==> test auc = 0.936275
Seed 13 ==> test auc = 0.930272
Seed 14 ==> test auc = 0.950706
Seed 15 ==> test auc = 0.919849
Seed 16 ==> test auc = 0.877632
Seed 17 ==> test auc = 0.871467
Seed 18 ==> test auc = 0.888110
Seed 19 ==> test auc = 0.905920
Seed 20 ==> test auc = 0.921946
Seed 21 ==> test auc = 0.928851
Seed 22 ==> test auc = 0.870864
Seed 23 ==> test auc = 0.865779
Seed 24 ==> test auc = 0.902960
Seed 25 ==> test auc = 0.946100
Seed 26 ==> test auc = 0.900540
Seed 27 ==> test auc = 0.923984
Seed 28 ==> test auc = 0.900345
Seed 29 ==> test auc = 0.905261
```

Overall test auc = 0.904481 +/- 0.031870 Elapsed time = 0:38:51

This run took less than an hour (~40 minutes), with 20 training epochs per fold. Increasing the number of folds and epochs seems to bring down the model's accuracy compared to the default settings. The overall AUC for the test data was AUC= 0.904481 +/- 0.031870, which still marginally outperformed the model in experiment 2 that took nearly 17 hours.