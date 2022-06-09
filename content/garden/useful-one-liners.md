---
title: "Useful one-liners"
date: 2022-05-15
lastmod: 2022-05-15
draft: false
garden_tags: ["cli", "general"]
summary: "A growing list of useful single line snippets of code"
status: "seeding"
---

![intro picture](/oneliner.png)

---
# Windows command line
---

#### Extract .tar.gz files

```python
tar -xvzf C:\path\to\file\filename.tar.gz -C C:\path\to\extraction\folder
```

---
# Python
---

#### Return Dict with count of unique list elements

```python
# Import Counter from the collections module
from collections import Counter
# Create a list
feature_list = ['repeat_region', 'snRNA', 'rRNA',
'miRNA', 'CDS', 'snRNA', 'rRNA']
# Create a dictionary of counts
c = Counter(feature_list)
# Results
Counter({'repeat_region': 1, 'snRNA': 1, 'rRNA': 1, 'miRNA': 1, 'CDS': 1})
```

#### Return subset of list with only unique values

```python
# Use set to remove duplicates and convert back to list
u_features = list(set(feature_list))
# Results
['repeat_region', 'miRNA', 'snRNA', 'rRNA', 'CDS']
```