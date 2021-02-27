import numpy as np
import math

def vocabulary(data):
	vocab = set()
	for sentence in data:
		for word in sentence:
			vocab.add(word)
	return vocab

def train(train_data, train_labels, vocab):
	d = len(vocab)
	pi = dict()
	theta = dict()
	for label in train_labels:
		if label in pi.keys():
			pi[label] += 1
		else:
			pi[label] = 1;
		theta[label] = dict.fromkeys(vocab,0);
	ex_size = sum(pi.values())
	for label in pi:
		pi[label] = pi[label]/ex_size

	denom = dict.fromkeys(pi,0)
	numer = dict.fromkeys(pi,[])
	for i in range(ex_size):
		current_class = train_labels[i]
		denom[current_class] += len(train_data[i])

		for word in train_data[i]:
			theta[current_class][word] += 1
	
	for c in theta.keys():
		for word in theta[c].keys():
			theta[c][word] = (1 + theta[c][word]) / (denom[c] + d)

	#print((theta[train_labels[0]]))
	#print(denom)
	return theta, pi

def test(theta, pi, vocab, test_data):
	scores = []
	test_size = len(test_data)
	classes = list(pi.keys())
	for i in range(test_size):
		example = []
		for k in range(len(classes)):
			total = 0
			for word in test_data[i]:
				if word in vocab:
					theta_est = theta[classes[k]][word]
					total += math.log(theta_est)#np.log([theta_est])[0]
			pi_est = pi[classes[k]]
			total += math.log(pi_est)#np.log([pi_est])[0]
			example.append((total,classes[k]))
		scores.append(example)
	return scores

if __name__ == '__main__':
	train_data, train_labels, test_data, test_labels = [], [], [], []
	
	with open("hw4_data/news/train_data.txt") as f:
		all_text = f.readlines()
		for sentence in all_text:
			sent_list = sentence.split(' ')
			new_sentence = []
			for word in sent_list:
				if word.isalnum():
					new_sentence.append(word)
				else:
					clean_word = ''.join(e for e in word if e.isalnum())
					new_sentence.append(clean_word)
			train_data.append(new_sentence)

	with open("hw4_data/news/test_data.txt") as f:
		all_text = f.readlines()
		for sentence in all_text:
			sent_list = sentence.split(' ')
			new_sentence = []
			for word in sent_list:
				if word.isalnum():
					new_sentence.append(word)
				else:
					clean_word = ''.join(e for e in word if e.isalnum())
					new_sentence.append(clean_word)
			test_data.append(new_sentence)

	with open("hw4_data/news/train_labels.txt") as f:
		all_text = f.read()
		train_labels = all_text.split('\n')
	
	with open("hw4_data/news/test_labels.txt") as f:
		all_text = f.read()
		test_labels = all_text.split('\n')

	vocab = vocabulary(train_data)
	theta, pi = train(train_data, train_labels, vocab)
	scores = test(theta, pi, vocab, test_data)

	test_size = len(test_data)
	correct_count = 0

	for i in range(test_size):
		pmax = -999999999
		for c in scores[i]:
			if c[0] > pmax:
				pmax = c[0]
				pred_class = c[1]
		if pred_class == test_labels[i]:
			correct_count += 1

	acc = correct_count/test_size
	print('Accuracy: %{percent:.f1}'.format(percent=acc*100))

