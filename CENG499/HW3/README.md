### Alperen Oğuz Çakmak - 2237162
# CENG499 | Homework 3 - Decision Trees - SVM 

## Part 1
For test purposes of this part we have 4 different modes. Also adding any second argument (`plot` for example) will create a tree as dot file and save it as pdf.

`$ python dt.py info` -> Construct a decision tree using information gain
`$ python dt.py gain` -> Construct a decision tree using gain ratio
`$ python dt.py gini` -> Construct a decision tree using gini index
`$ python dt.py gain-c` -> Apply pre-pruning to a tree constructed with gain ratio

## Part 2
SVM takes several arguments such as task number, kernel, c and gamma. First argument is task number from 1 to 4.
##### Task1
`$ python svm.py 1 0.1` -> Trains linear classifier with given C and plots the result

##### Task2
Possible kernels are `linear`, `rbf`, `poly` and `sigmoid`
`$ python svm.py 2 rbf` -> Trains the classifier with given kernel and plots the result

##### Task3
Possible kernels are `linear`, `rbf`, `poly` and `sigmoid`
`$ python svm.py 3 linear 1` -> Trains classifier with linear kernel and given c value
`$ python svm.py 3 rbf 1 0.001` -> Trains classifier with given kernel, c and gamma values

##### Task4
Possible modes are `normal`, `over`, `under` and `cw`
`$ python svm.py 4 normal` -> Trains the svm and prints the accuracy and confusion matrix
`$ python svm.py 4 over` -> Uses oversampling for before training the svm

