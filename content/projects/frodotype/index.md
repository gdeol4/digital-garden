---
title: "Frodotype - text generator"
date: 2020-01-17
draft: false
project_tags: ["Web app", "Tensorflow", NLP]
status: "evergreen"
weight: 1
summary: "A GPT2 generator finetuned on fantasy novels. Built using Tensorflow, Chart.js, GCP and Docker"
links:
    external_link:
        text: "Link to live app"
        icon: "fas fa-external-link-alt"
        href: "https://frodotype.net"
        weight: 1
    another_link:
        text: "The github repo"
        icon: "fab alt brands fa-github"
        href: "https://github.com/gdeol4/frodotype"
        weight: 2
---

# A fantasy tuned text generator

Frodotye is a webapp that generates text which grammatically correct and logical most of the time. This project aims to fine-tune the GPT2 model using the gpt-2-simple python package. OpenAI offers the use of 3 of their models, however the size and complexity of these models makes it diffcult to train on consumer hardware. Therefor, the smallest (117m) parameter model was used in conjucntion with a Google Deep Leanring VM to retrain the model.The data used to fine-tuned the model consists of 102 fantasy novels by 20 authors. This data was used for 2 reasons:

1. Fantasy is the genre I read the most.
2. I own the ebooks used in this project.

Frodotype was built using Tensorflow, docker, Google Cloud Platform, Javascript, Chart.js, and Bulma CSS. A bulk of the heavy lifting in python was done using Max Wolf's gpt-2-simple project.

## Gathering the data

The 102 books used were converted from the Amazon Kindle format .azw to plain text files. This processes included stripping all images, formattting, and hyperlinks. The books were then manually stripped of their table of contents, appencices, and glossaries. The final text file used to retrain the model was put together using the following commands:

1.Concatenate the files:

```bash
$ for f in *.txt do (cat "${f}"; echo) >> unprocessed.txt; done
```

2. Deleting all none ASCII charcaters:

```bash
$ LC_ALL=C tr -dc '\0-\177' < unprocessed.txt > processed.txt
```

3. Removing numbers and dashes from the text:

```bash
$ tr -d '[0-9-]' < processed.txt > final.txt
```

Additional processing is done in the text-analysis notebook.

## Training the model

The model was trained on a Google Deep Learning VM using a Tesla K80 GPU, TensorFlow 1.15, and CUDA 10.0.

The model was retrained using gpt-2-simple, a python package that eases the process of tweeking hyperparameters. The model was trained for three differeing lengths. The one used in this app was trianed for 45,000 steps or approximattly 90 hours. Two additional models were trained at 25,000 steps and 80,000 steps. The smaller of the two models had a much higher loss value, while the larger model had a simillar loss that began to increase towards the end.