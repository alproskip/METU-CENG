import draw
import sys
import numpy as np 
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn import metrics
from sklearn.model_selection import cross_val_score

def task1():
	c = float(sys.argv[2])
	train_data = np.load('hw3_data/linsep/train_data.npy')
	train_labels = np.load('hw3_data/linsep/train_labels.npy')

	cf = SVC(c, kernel = 'linear')
	cf.fit(train_data, train_labels)

	x_min, x_max, y_min, y_max = min(train_data[:,0]), max(train_data[:,0]), min(train_data[:,1]), max(train_data[:,1])
	draw.draw_svm(cf, train_data, train_labels, x_min, x_max, y_min, y_max)

def task2():
	kernel = sys.argv[2]
	train_data = np.load('hw3_data/nonlinsep/train_data.npy')
	train_labels = np.load('hw3_data/nonlinsep/train_labels.npy')

	cf = SVC(kernel = kernel)
	cf.fit(train_data, train_labels)

	x_min, x_max, y_min, y_max = min(train_data[:,0]), max(train_data[:,0]), min(train_data[:,1]), max(train_data[:,1])
	draw.draw_svm(cf, train_data, train_labels, x_min, x_max, y_min, y_max)

def task3():
	kernel, c = sys.argv[2], float(sys.argv[3])
	if kernel != "linear":
		gamma = float(sys.argv[4])

	train_data, train_labels, test_data, test_labels = np.load('hw3_data/catdog/train_data.npy'), np.load('hw3_data/catdog/train_labels.npy'), np.load('hw3_data/catdog/test_data.npy'), np.load('hw3_data/catdog/test_labels.npy')
	#NORMALIZE
	n_train_data = np.interp(train_data, (train_data.min(), train_data.max()), (-1, +1))
	n_test_data = np.interp(test_data, (test_data.min(), test_data.max()), (-1, +1))

	if kernel == "linear":
		cf = SVC(C = c, kernel = kernel)	
	else:
		cf = SVC(C = c, kernel = kernel, gamma=gamma)

	cf.fit(n_train_data, train_labels)
	#CROSS VALIDATION
	pred = cf.predict(n_test_data)
	acc = metrics.accuracy_score(test_labels, pred)
	print("accuracy is ", acc)

def task4():
	mode = sys.argv[2]
	train_data, train_labels, test_data, test_labels = np.load('hw3_data/catdogimba/train_data.npy'), np.load('hw3_data/catdogimba/train_labels.npy'), np.load('hw3_data/catdogimba/test_data.npy'), np.load('hw3_data/catdogimba/test_labels.npy')
	
	cf = SVC(C = 1, kernel = 'rbf')

	if mode == 'over':
		os_data = []
		for j in range(len(train_labels)):
			if train_labels[j]==0:
				os_data.append(train_data[j])
		os_data *= 4
		os_labels = np.zeros(len(os_data))
		os_data = np.asarray(os_data)
		train_data = np.append(train_data, os_data, axis=0)
		train_labels = np.append(train_labels, os_labels)

	elif mode == 'under':
		idx = []
		for i in range(len(train_labels)):
			if train_labels[i]==1:
				idx.append(i)
			if len(idx) > 625:
				break
		train_data = np.delete(train_data, idx, 0)
		train_labels = np.delete(train_labels, idx, 0)
	
	elif mode == 'cw':
		cf = SVC(C = 1, kernel= 'rbf', class_weight='balanced')	

	n_train_data = np.interp(train_data, (train_data.min(), train_data.max()), (-1, +1))
	n_test_data = np.interp(test_data, (test_data.min(), test_data.max()), (-1, +1))

	cf.fit(n_train_data, train_labels)
	pred = cf.predict(n_test_data)
	acc = metrics.accuracy_score(test_labels, pred)
	print("accuracy is ", acc)
	print("confusion matrix is \n", metrics.confusion_matrix(test_labels, pred))


if __name__ == '__main__':
	task = sys.argv[1]
	if task == '1':
		task1()
	elif task == '2':
		task2()
	elif task == '3':
		task3()
	elif task == '4':
		task4()