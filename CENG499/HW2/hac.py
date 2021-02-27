import matplotlib.pyplot as plt 
import numpy as np
import sys
from math import *

def plot(data, n):
	
	c1_x = [p[0] for p in data[0]]
	c1_y = [p[1] for p in data[0]]
	
	c2_x = [p[0] for p in data[1]]
	c2_y = [p[1] for p in data[1]]

	plt.plot(c1_x, c1_y, 'ro', color = '#cc9900', marker='o', markersize=3)
	plt.plot(c2_x, c2_y, 'ro', color = '#cc00cc', marker='o', markersize=3)

	if n==4:
		c3_x = [p[0] for p in data[2]]
		c3_y = [p[1] for p in data[2]]
		
		c4_x = [p[0] for p in data[3]]
		c4_y = [p[1] for p in data[3]]

		plt.plot(c3_x, c3_y, 'ro', color = '#339966', marker='o', markersize=3)
		plt.plot(c4_x, c4_y, 'ro', color = '#993333', marker='o', markersize=3)

	plt.title('Data %d' % n) 
	plt.show() 

def dist(t, p):
	a = (t-p)**2
	s = a.sum()
	return sqrt(s)

def single_linkage(c1, c2):
	d_min = 999
	for i in range(len(c1)):
		for j in range(len(c2)):
			d = dist(c1[i],c2[j])
			if d<d_min:
				d_min = d
	return d_min

def complete_linkage(c1, c2):
	d_max = 0
	for i in range(len(c1)):
		for j in range(len(c2)):
			d = dist(c1[i],c2[j])
			if d>d_max:
				d_max = d
	return d_max

def average_linkage(c1, c2):
	total = 0
	s1, s2 = len(c1), len(c2)
	for i in range(s1):
		for j in range(s2):
			d = dist(c1[i],c2[j])
			total += d
	final_d = total/(s1*s2)
	return final_d

def centroid(c1, c2):
	s1, s2 = len(c1), len(c2)
	c1_x = [c[0] for c in c1]
	c1_y = [c[1] for c in c1]
	c2_x = [c[0] for c in c2]
	c2_y = [c[1] for c in c2]
	center1 = np.asarray([sum(c1_x)/s1, sum(c1_y)/s1]) 
	center2 = np.asarray([sum(c2_x)/s2, sum(c2_y)/s2]) 
	return dist(center1, center2)

def hac(data, linkage_f, k):
	
	distance_func = linkage_f
	active_set = [[data[x]] for x in range(len(data))]

	while len(active_set)>k:
		d_min = 999
		min_idx = [0,0]
		size = len(active_set)
		for i in range(size):
			for j in range(i+1,size):
				d = distance_func(active_set[i], active_set[j])
				if d<d_min:
					d_min = d 
					min_idx = [i,j]
					
		active_set[min_idx[0]] += active_set[min_idx[1]]
		active_set.pop(min_idx[1]) 

	return active_set


if __name__ == '__main__':

	mode = sys.argv[1]
	data_num = sys.argv[2]
	k = 2
	
	if data_num == "1":
		data = np.load("hw2_data/hac/data1.npy")
	elif data_num == "2":
		data = np.load("hw2_data/hac/data2.npy")
	elif data_num == "3":
		data = np.load("hw2_data/hac/data3.npy")
	elif data_num == "4":
		data = np.load("hw2_data/hac/data4.npy")
		k = 4

	if mode == "s":
		clusters = hac(data, single_linkage, k)
	elif mode == "c":
		clusters = hac(data, complete_linkage, k)
	elif mode == "a":
		clusters = hac(data, average_linkage, k)
	elif mode == "cen":
		clusters = hac(data, centroid, k)


	plot(clusters, int(data_num))
