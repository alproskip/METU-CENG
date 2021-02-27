import matplotlib.pyplot as plt 
import numpy as np
import collections
import sys
from math import *

def plot(knns, accs):
	plt.plot(knns, accs, label = "Accuracies", color='blue', linestyle='dashed', linewidth = 3, marker='o', markerfacecolor='red', markersize=12)
	plt.xticks(np.arange(min(knns), max(knns)+1, 4.0)) 
	plt.xlabel('k of KNN') 
	plt.ylabel('Accuracies') 
	plt.title('Accuracy of Different k values') 
	plt.legend() 
	plt.show() 

def dist(t, p):
	a = (t-p)**2
	s = a.sum()
	return sqrt(s)

def chunkify(fold_idx,data):
	fold = data[fold_idx*25:(fold_idx+1)*25]
	left = data[:fold_idx*25]
	right = data[(fold_idx+1)*25:]
	rest = np.concatenate((left,right),axis=0)
	return fold, rest

def test_for_best_k(k):
	train_data = np.load('hw2_data/knn/train_data.npy')
	train_labels = np.load('hw2_data/knn/train_labels.npy')
	test_data = np.load('hw2_data/knn/test_data.npy')
	test_labels = np.load('hw2_data/knn/test_labels.npy')

	accs = []

	distances = []
	preds = []
	for i in test_data:
		distances = [dist(x,i) for x in train_data]
		distances = list(enumerate(distances))
		
		k_min = []
		
		for j in range(k):
			curr_min = min(distances,key=lambda t: t[1])
			distances.pop(distances.index(curr_min))
			k_min.append(curr_min[0])

		c_idx = []

		for n in k_min:
			c_idx.append(train_labels[n])
		preds.append(collections.Counter(c_idx).most_common(1)[0][0])			
		
	preds = np.asarray(preds)
	diff = np.asarray(preds) - test_labels
	acc = len([x for x in diff if x==0])/200.

	return acc

def train_for_ks(k):

	train_data = np.load('hw2_data/knn/train_data.npy')
	train_labels = np.load('hw2_data/knn/train_labels.npy')

	f_accs = []
	cv = 10

	for i in range(cv):

		distances = []
		fold, train = chunkify(i,train_data)
		label_fold, labels = chunkify(i,train_labels) 

		preds = []
		for i in fold:
			distances = [dist(x,i) for x in train]
			distances = list(enumerate(distances))
			
			k_min = []
			
			for j in range(k):
				curr_min = min(distances,key=lambda t: t[1])
				distances.pop(distances.index(curr_min))
				k_min.append(curr_min[0])
			
			c_idx = []

			for n in k_min:
				c_idx.append(labels[n])
			preds.append(collections.Counter(c_idx).most_common(1)[0][0])			
			

		preds = np.asarray(preds)

		diff = preds - label_fold

		acc = len([x for x in diff if x==0])/25.
		#print(acc)
		f_accs.append(acc)

	#print(f_accs)
	return sum(f_accs)/len(f_accs)

if __name__ == '__main__':
	
	mode = sys.argv[1]

	knns = [x for x in list(range(200)) if x%2==1]
	
	if mode == "plot":
		accuracies = [train_for_ks(k) for k in knns]
		plot(knns, accuracies)

	elif mode == "accs":
		accuracies = [train_for_ks(k) for k in knns]
		max_idx = accuracies.index(max(accuracies))*2+1
		print("Most accurate k is",max_idx,"with the accuracy of", max(accuracies))

	elif mode == "best":
		k = 11
		acc = test_for_best_k(k)
		acc *= 100
		print("Accuracy of test data for best k is %"+"%d" %acc)

