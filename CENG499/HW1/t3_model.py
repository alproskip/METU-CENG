import torch
import torch.nn as nn
import torchvision.transforms as T
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from dataset import myDataset

class T3Model(nn.Module):
	def __init__(self):
		super(T3Model, self).__init__()
		self.fc1 = nn.Linear(64*32, 1024)
		self.fc2 = nn.Linear(1024, 1024)
		self.fc3 = nn.Linear(1024, 100)

	def forward(self, x):
		x = torch.flatten(x,1)
		x = self.fc1(x)
		x = F.tanh(x)
		x = self.fc2(x)
		x = F.tanh(x)
		x = self.fc3(x)
		x = torch.log_softmax(x, dim=1)
		return x

if __name__ == "__main__":
	transforms = T.Compose([
		T.ToTensor(),
		T.Normalize((0.5),(0.5)),
	])
	dataset = myDataset('data', 'train', transforms)
	dataloader = DataLoader(dataset, 
		batch_size=64, 
		shuffle=True,
		num_workers=4)

	model = T3Model()

	for images, labels in dataloader:
		pred = model(images)
		print(pred)
		exit()