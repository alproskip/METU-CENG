### Alperen Oğuz Çakmak - 2237162
# CENG499 | Homework 2 - KNN, K-means, HAC
Hello again! Sorry about the previous homework. This time I managed to make it easy to test so here we go. 

## Task 1
For test purposes of this task, knn.py file has 3 modes.

`$ python knn.py plot` -> Plots accuracies for different k values in range 1 - 199 

`$ python knn.py accs` -> Prints the most accurate k value and it's accuracy

`$ python knn.py best` -> Prints accuracy of test data for most accurate k value

## Task 2
kmeans has 2 modes and takes 2 arguments, first one is mode second one is cluster number

`$ python kmeans.py obj 3` -> Plots obj vs k for elbow method for cluster 3

`$ python kmeans.py cluster 2` -> Plots indicated cluster (Cluster 2)

## Task 3
hac has 4 modes and takes 2 arguments (mode and data number)

`$ python hac.py s 1` -> Use data 1 and single-linkage criterion 

`$ python hac.py c 2` -> Use data 2 and complete-linkage criterion 

`$ python hac.py a 3` -> Use data 3 and average-linkage criterion 

`$ python hac.py cen 4` -> Use data 4 and the centroid criterion 
