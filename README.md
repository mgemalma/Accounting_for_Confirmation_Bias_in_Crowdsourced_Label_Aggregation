<p align="center">

# Accounting for Confirmation Bias in Crowdsourced Label Aggregation

This is a repository for the data and code for the paper "Accounting for Confirmation Bias in Crowdsourced Label Aggregation".

## Repository Structure

* `Algorithm` contains the code of the proposed algorithm implemented by C. GSL library version 2.6 is used widely in the code and this library needs to be downloaded.

* `Datasets` contains 1: `worker_answers.txt` the cleaned data that contains the answers of the workers' and is ready to be given to the algorithm 2:`worker_data.csv` contains the raw data of the workers who passed the attention check.
## Algorithm

## INSTALLATION INSTRUCTIONS

1. Install GSL, if not already installed.
2. Modify Makefile inside the `Algorithm` directory to point to the correct locations of the GSL library.
3. Run make.

To run the algorithm;
	./em \<Name of the dataset file\>

The algorithm will produce two files:

`results.txt` which will include the label predictions and
`worker_parameter_results.txt` which will include the inferred labeler and statement parameters.

`em.c` contains the main function and the high level EM code.

`prob_functions.c` contains the EM code and its helper functions.

`data.c` contains the code that manages the dataset, it receives the input dataset and outputs the two prediction txts from above.

## Supplementary Materials

Supplemental Material.docx contains the full derivations for our paper and also provides the political statements we used in our experiment.

## Data Sample

Here attached a data sample from the raw data `worker_data.csv`. Each row is a piece of record on a worker that answered all of the tasks.

| worker_index | stance | age | gender | education | state | party | gc_stance | S1 | S2 | S3 | S4 | S5 | S6 | S7 | S8 | S9 | S10 | S11 | S12
| -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- |
| 0 | 2 | 4.0 | 1.0 | 8.0 | AZ | 1.0 | 3.0 | -1 | 0 | 0 | 0 | 0 | -1 | 0 | 0 | 0 | 0 | 1 | -1 |

In all of the below fields, the values that are not in the specifically defined range are data loss.

* `worker_index` Int. The unique ID for each worker assigned.
* `stance` Int. [1 - 7] Workers’ self-reported political stance. The smaller the stance is, the more worker is holding liberal values and vice-versa for conservative values.
* `age` Int. [1 - 9] The age range of the worker. Under 12 years old (1), 12-17 years old (2), 18-24 years old (3), 25-34 years old (4), 35-44 years old (5), 45-54 years old (6), 55-64 years old (7), 65-74 years old (8) and 75 years or older (9).
* `gender` Int. [1 - 3] The gender of the worker. Male (1), Female (2) and Other (3).
* `education` Int. [1 - 11] The worker's highest degree or level of school that she have completed. No schooling completed (1), Nursery school to 8th grade (2), Some high school, no diploma (3), High school graduate, diploma or the equivalent (for example: GED) (4), Some college credit, no degree (5), Trade/technical/vocational training (6), Associate degree (7), Bachelor’s degree (8), Master’s degree (9), Professional degree (10) and Doctorate degree (11).
* `state` Str. The initials of the state that the worker is completing the task from.
* `party` Int. [1 - 3] The worker's political party of affiliation. Democrat (1), Republican (2) and Independent (3).
* `gc_stance` Int. [1 - 7] The worker's stance for Gun Control, the smaller the stance is, the more worker is strongly against gun control and vice-versa for strongly supporting gun control.
* `S1 - S12` Int. Each of them respectively reports the answer given by the worker to that task. `1` indicates "Opinion" `0` indicates "Factual" and `-1` indicates "I don't know".

The `worker_answers.txt` is in the format that the algorithm is expecting and contains the following.

NOTE: We treat every "I don't know" label as a missing label and do not include it.

The first row contains;
| Total number of Labels | Number of Labelers | Number of Statements |
| -------- | -------- | -------- |
| 1213 | 110 | 12 |
* `Total number of Labels` Int. The total number of labels acquired from the workers.
* `Number of Labelers` Int. The total number of workers who labeled the tasks.
* `Number of Statements` Int. The number of statements or the number of tasks that we expect to find ground-truth value of.

All the other rows;
| Statement ID | Labeler ID | Label |
| -------- | -------- | -------- |
| 0 | 3 | 1 |
* `Statement ID` Int. Index of the statement.
* `Labeler ID` Int. Index of the labeler (worker).
* `Label` Int. The label given by the worker to the statement or the task.

## Licensing & Citing
When using or building upon the data in an academic publication, please consider citing as follows:

Gemalmaz, M., & Yin, M. (2021, August). Accounting for Confirmation Bias in Crowdsourced Label Aggregation. In *Proceedings of the Thirtieth International Joint Conference on Artificial Intelligence*
