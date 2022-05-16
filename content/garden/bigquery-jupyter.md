---
title: "Connecting BigQuery and Jupyter"
date: 2022-05-16
lastmod: 2022-05-16
draft: false
garden_tags: ["python", "SQL", "data"]
summary: "A reference on how to setup a GCP account to execute BigQuery SQL in Jupyter "
status: "evergreen"
---

### Importing libraries

```python
# for new environments
# pip install --upgrade 'google-cloud-bigquery[bqstorage,pandas]'

import os
import warnings
warnings.filterwarnings('ignore')
```

### Google cloud platform setup

##### Create a project

Go to the google cloud platform console and either choose an existing project or create a new one

![bigquery1](/bigquery1.png)

##### Create a service account

Navigate to the left side menu and proceed to "IAM & Admin" then to "Service Accounts"

![bigquery2](/bigquery2.png)


Set a service account name:

![bigquery3](/bigquery3.png)


Set the "Role" to "Owner":


![bigquery4](/bigquery4.png)


Confirm that the account was created and click on Actions and then "Manage Keys":


![bigquery5](/bigquery5.png)


Create a JSON private key



![bigquery6](/bigquery6.png)



Navigate to API's and Services and enable the BigQuery API:


![bigquery7](/bigquery7.png)



Set the environment variable:


```python
os.environ["GOOGLE_APPLICATION_CREDENTIALS"]="C:/Users/gurka/Downloads/bigquery_key.json"
```


### Executing BigQuery Jupyter cells


##### Loading magic command


The BigQuery client library for Python provides a magic command that lets you run queries with minimal code. To load the magic commands from the client library, paste the following code into the first cell of the notebook.

```python
%load_ext google.cloud.bigquery
```

##### Running a test on public data

The BigQuery client library for Python provides a cell magic, **%%bigquery**, which runs a SQL query and returns the results as a Pandas DataFrame. Enter the following code in the next cell to return total births by year:


```python
%%bigquery
SELECT
    source_year AS year,
    COUNT(is_male) AS birth_count
FROM `bigquery-public-data.samples.natality`
GROUP BY year
ORDER BY year DESC
LIMIT 15

Query complete after 0.02s: 100%|████████████████████████████████████████████████████| 1/1 [00:00<00:00, 999.12query/s]
Downloading: 100%|███████████████████████████████████████████████████████████████████| 15/15 [00:01<00:00,  9.74rows/s]


year	birth_count
0	2008	4255156
1	2007	4324008
2	2006	4273225
3	2005	4145619
4	2004	4118907
5	2003	4096092
6	2002	4027376
7	2001	4031531
8	2000	4063823
9	1999	3963465
10	1998	3945192
11	1997	3884329
12	1996	3894874
13	1995	3903012
14	1994	3956925
```