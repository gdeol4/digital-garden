---
title: "A primer on Apache Airflow"
date: 2021-04-20
lastmod: 2021-04-20
draft: false
garden_tags: ["python", "data"]
summary: "Airflow basics for reference, and some pictures I drew "
status: "evergreen"
---

# What is airflow?
   
Airflow is a platform used to author, schedule, and monitor workflows.
It’s essentially a queuing system that runs on a metadata database and a scheduler that runs tasks. Workflows are written as Directed Acyclic Graphs (DAGs). A workflow and DAG are interchangeable.

# What are DAGs?
A DAG is a collection of tasks you want to run and are organized in a way that illustrates dependencies and relationships between tasks.
The image below shows how a DAG is a unidirectional, acyclic graph, where each node in the graph is a task and edges define dependencies among tasks. There is no case where you should be able to go backwards from a forward node to one that's already been executed.

![Scenario 1: Across columns](/airflow1_1.jpeg)

A DAG can be broken up into smaller and smaller jobs and gives the user full control by generating dynamic pipelines written in code. Airflow DAGs are also extensible and can scale. DAGs are powerful because they allow for collaborative, manageable, and testable workflows. A bonus is that Airflow is developed in python and can interface with any python API.

![Scenario 1: Across columns](/airflow1_2.jpeg)

The image above shows how Airflow divides the tasks into branches so that if one fails, there is still output from the other. Also, the processing time is reduced as parallel computing occurs. The chances of failure should decrease overall as each task is independent.

## How are tasks executed?
An operator represents a single task in a workflow that helps carry out your task (running a python function for example).
Operators determine what actually gets to be done when your dag runs.
A task is an operator when instantiated. It is something on which the worker works upon.

## Airflow Architecture
    
![Scenario 1: Across columns](/airflow1_3.jpeg)
   
- Metadata — is a relational database with info on task state, such as the top ten tasks consuming the most memory, it contains all data pertaining to jobs currently running as well as historical data.
- Scheduler — decides which task to run, when, and in what order.
- Web server— the UI which is essentially a flask app that talks to the metadata.
- Executor — performs the task at ground level. The executor is a message queuing process that figures out which workers will execute which tasks. The default is the sequential executor — which cannot run tasks in parallel — meaning it can’t be used for production level code. The local executor can be used too which will run tasks till all resources on the server are at capacity. This is good for a moderate amount of DAGs. Both of these are used in single node clusters and therefore cannot be used to scaled.
- Multi node clusters — have the same components and only the scheduler and web server are placed in the same node (master), the workers are placed in a separate instance. This set up works well because it allows for scaling by letting you add more multi-node clusters (celery is the executor of choice here for python).

If you're not dealing with terabytes of data then it's better to have the scheduler, web server, and executor together in the master node/cluster. The downside is that this single cluster approach runs everything on the same machine, so if you make a change to a DAG/scheduler, then you need to restart the entire workflow — even tasks that were in the process of executing. Celery avoids this.

![Scenario 1: Across columns](/airflow1_4.jpeg)

If you do build a distributed workflow with celery then a queuing system component is needed (like Redis). For local workflows, the queuing is handled by the system.

## The life cycle of a task

1. The scheduler periodically checks the DAG folder to see if there are any DAGS that need to be run.
2. If any DAGS are found pending execution, the scheduler creates a diagram for it, which is an instantiation of a DAG in real time.
3. The scheduler will update the DAG state to running in the metadata and the tasks will execute.
4. The scheduler then reads the DAG and puts the tasks in order of execution into the queuing system in the form of a message. Each message contains info like DAG ID, TASK ID, and function to be executed.
5. The status of these tasks changes to queued at that point.
6. The executor then begins to execute tasks and sends fail/success messages for the tasks to the metadata.
7. The scheduler finally updates the status of the diagram when all tasks have run to success or failure.
