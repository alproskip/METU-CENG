import matplotlib.pyplot as plt 
import numpy as np 
import sys
from math import *

def plot_elbow(xs, ys, cn):
	plt.plot(xs, ys, label = "Accuracies", color='blue', linestyle='dashed', linewidth = 3, marker='o', markerfacecolor='red', markersize=12)
	plt.xticks(np.arange(min(xs), max(xs)+1, 1.0)) 
	plt.xlabel('k-means iterations') 
	plt.ylabel('Objective') 
	plt.title('Cluster ' + cn) 
	plt.legend() 
	plt.show() 

def plot_cluster(centers, data, k):
	class1_x = [p[1][0] for p in data if p[0]==0]
	class1_y = [p[1][1] for p in data if p[0]==0]

	class2_x = [p[1][0] for p in data if p[0]==1]
	class2_y = [p[1][1] for p in data if p[0]==1]

	plt.plot(class1_x, class1_y, 'ro', color = '#cc9900', marker='o', markersize=1)
	plt.plot(class2_x, class2_y, 'ro', color = '#333399', marker='o', markersize=1)
	plt.plot([centers[0][0]], [centers[0][1]], 'ro', color = '#cc9900', marker='^', markersize=20, markeredgewidth=2, markeredgecolor='#66ff33')
	plt.plot([centers[1][0]], [centers[1][1]], 'ro', color = '#333399', marker='^', markersize=20, markeredgewidth=2, markeredgecolor='#66ff33')

	if k>2:
		class3_x = [p[1][0] for p in data if p[0]==2]
		class3_y = [p[1][1] for p in data if p[0]==2]
		plt.plot(class3_x, class3_y, 'ro', color = '#cc00cc', marker='o', markersize=1)
		plt.plot([centers[2][0]], [centers[2][1]], 'ro', color = '#550066', marker='^', markersize=20, markeredgewidth=2, markeredgecolor='#66ff33')
	if k>3:	
		class4_x = [p[1][0] for p in data if p[0]==3]
		class4_y = [p[1][1] for p in data if p[0]==3]
		plt.plot(class4_x, class4_y, 'ro', color = '#339966', marker='o', markersize=1)
		plt.plot([centers[3][0]], [centers[3][1]], 'ro', color = '#206040', marker='^', markersize=20, markeredgewidth=2, markeredgecolor='#66ff33')
	if k>4:
		class5_x = [p[1][0] for p in data if p[0]==4]
		class5_y = [p[1][1] for p in data if p[0]==4]
		plt.plot(class5_x, class5_y, 'ro', color = '#993333', marker='o', markersize=1)
		plt.plot([centers[4][0]], [centers[4][1]], 'ro', color = '#993333', marker='^', markersize=20, markeredgewidth=2, markeredgecolor='#66ff33')

	plt.title('Cluster %d' % (k-1)) 
	plt.show() 

def dist(t, p):
	a = (t-p)**2
	s = a.sum()
	return sqrt(s)

def obj_func(points,centers):
	distances = [dist(centers[p[0]], p[1]) for p in points]
	total = 0
	for d in distances:
		total += (d**2)
	return total

def classify(p, centers):
	if len(centers) == 1:
		return 0
	lst = [dist(p,c) for c in centers]
	return(lst.index(min(lst)))

def init_centroids(clustering, k):
	x_max, x_min = max(clustering[:,0]), min(clustering[:,0])
	y_max, y_min = max(clustering[:,1]), min(clustering[:,1])

	dx = (x_max - x_min)/(k+1)
	dy = (y_max - y_min)/(k+1)

	centroid_x = [(i+1)*dx+x_min for i in range(k)]
	centroid_y = [(i+1)*dy+y_min for i in range(k)]

	return [(x,y) for x,y in zip(centroid_x,centroid_y)]

def kmeans(c, k, ret_p):
	centers = init_centroids(c, k)
	classes = [(classify(p,centers),p) for p in c]
	of_values = []
	for i in range(10):
		for x in range(k):
			temp_cluster = np.asarray([p[1] for p in classes if p[0]==x])
			
			if temp_cluster.shape[0] == 0:
				continue
			
			x_col = temp_cluster[:,0]
			y_col = temp_cluster[:,1]

			new_x = sum(x_col)/len(x_col)
			new_y = sum(y_col)/len(y_col)

			centers[x] = (new_x,new_y)

		classes = [(classify(p,centers),p) for p in c]
		of_values.append(obj_func(classes,centers))

	avg_obj = sum(of_values)/10.

	if ret_p:
		return centers, classes

	return centers, avg_obj

if __name__ == '__main__':

	c1 = np.load("hw2_data/kmeans/clustering1.npy")
	c2 = np.load("hw2_data/kmeans/clustering2.npy")
	c3 = np.load("hw2_data/kmeans/clustering3.npy")
	c4 = np.load("hw2_data/kmeans/clustering4.npy")

	mode = sys.argv[1]
	n = sys.argv[2]
	if mode == "obj":
		k_list = [x+1 for x in range(10)]
		if n=="1":
			cluster_1 = np.asarray([kmeans(c1, k, False) for k in k_list])
			plot_elbow(k_list, cluster_1[:,1],n)
		elif n=="2":
			cluster_2 = np.asarray([kmeans(c2, k, False) for k in k_list])
			plot_elbow(k_list, cluster_2[:,1],n)
		elif n=="3":
			cluster_3 = np.asarray([kmeans(c3, k, False) for k in k_list])
			plot_elbow(k_list, cluster_3[:,1],n)
		elif n=="4":
			cluster_4 = np.asarray([kmeans(c4, k, False) for k in k_list])
			plot_elbow(k_list, cluster_4[:,1],n)

	elif mode == "cluster":
		if n=="1":
			k = 2
			centers, points = kmeans(c1,k, True)
			plot_cluster(centers, points,k)
		elif n=="2":
			k = 3
			centers, points = kmeans(c2,k, True)
			plot_cluster(centers, points,k)
		elif n=="3":
			k = 4
			centers, points = kmeans(c3,k, True)
			plot_cluster(centers, points,k)
		elif n=="4":
			k = 5
			centers, points = kmeans(c4,k, True)
			plot_cluster(centers, points,k)

