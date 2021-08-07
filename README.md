# Accounting_for_Confirmation_Bias_in_Crowdsourced_Label_Aggregation
To run the algorithm
	- make clean
	- make
	./em data.txt

The output will produce two files:
TransCend.txt which will include the label predictions and
TransCend_parameters.txt which will include the labeler and statement parameters.

em.c is the starter file, has the main function and the high level EM code.
prob_functions.c contains the EM function codes and its helper functions.
data.c manages the receiving the input dataset and outputting the two prediction txts from above.

Demo dataset 'data.txt' is given, it's real ground truth is given at the ground_truth.txt
data.txt columns are as follows:

First Line: 'Total number of Labels' 'Number of Labelers' 'Number of Statements' 
All the other lines: 'Statement ID' 'Labeler ID' 'Label'
