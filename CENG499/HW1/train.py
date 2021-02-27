import torch
import torchvision.transforms as T
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from dataset import myDataset
from r3_model import R3Model
from r2_model import R2Model
from t3_model import T3Model
from t2_model import T2Model
from torch.utils.data import random_split

def train(model, optimizer, train_loader, valid_loader, epochs, device, lr):

	plot_t = []
	plot_v = []

	tr_loss = []
	val_loss = []
	stop = 0
	min_val_loss = 999
	stop_count = 0
	for epoch_idx in range(epochs):
	
		model.train()
	
		count = 0
		sum_loss = 0
	
		for images, labels in train_loader:
			images = images.to(device)
			labels = labels.to(device)
			optimizer.zero_grad()
			pred = model(images)
			loss = F.nll_loss(pred, labels)
			loss.backward()
			optimizer.step()
			count += 1
			sum_loss += loss.item()

			tr_loss.append(loss.item())

		with torch.no_grad():
			model.eval()
			correct = tot = 0
			for images, labels in valid_loader:
				images = images.to(device)
				labels = labels.to(device)
				pred = model(images)
				_, p_labels = torch.max(pred, 1)
				correct += (p_labels == labels).sum()
				tot += labels.size(0)
				loss = F.nll_loss(pred, labels)
				val_loss.append(loss.item())

		trl = sum(tr_loss)/len(tr_loss)
		vll = sum(val_loss)/len(val_loss)

		if vll < min_val_loss:
			min_val_loss = vll

		else:
			min_val_loss = vll
			stop_count += 1

		if stop_count > 5:
			#print("Early Stop")
			stop = 1
			break
		plot_t.append(trl)
		plot_v.append(vll)

		tr_loss = []
		val_loss = []
	print(plot_t)
	print(plot_v)
	model_name = 'myModel' #'myModel-Relu'+str(epochs)+'-'+str(lr)
	torch.save(model.state_dict(), 'myModel')
	#print(model_name,' saved')


def main():
	lr = 0.0003
	mod = 'R'
	epochs = 25
	layer_num = 3
	batch_size = 64

	if mod == 'R':
		print("Layer Number: %d, Batch_size: %d, Learning Rate: %.4f, Activation Function: ReLu " % (layer_num, batch_size, lr))
	else:
		print("Layer Number: %d, Batch_size: %d, Learning Rate: %.4f, Activation Function: Tanh " % (layer_num, batch_size, lr))

	use_cuda = True
	device = torch.device('cuda' if use_cuda else 'cpu')
	torch.manual_seed(1234)

	transforms = T.Compose([
		T.ToTensor(),
		T.Normalize((0.5),(0.5)),
	])
	train_dataset = myDataset('data', 'train', transforms)
	train_set, valid_set = random_split(train_dataset, [8000,2000])
	
	train_loader = DataLoader(train_set, 
		batch_size=batch_size, 
		shuffle=True,
		num_workers=4)

	valid_loader = DataLoader(valid_set,
		batch_size=batch_size,
		shuffle=True,
		num_workers=4)
	if layer_num == 2 and mod == 'T':
		model = T2Model()
	elif layer_num == 3 and mod == 'T':
		model = T3Model()
	elif layer_num == 2 and mod == 'R':
		model = R2Model()
	else:
		model = R3Model()


	model = model.to(device)

	optimizer = torch.optim.Adam(model.parameters(), lr=lr)

	train(model, optimizer, train_loader, valid_loader, epochs, device, lr)
	

if __name__ == '__main__':
	main()