import os
from torch.utils.data import Dataset, DataLoader
from PIL import Image
import torchvision.transforms as T

class myDataset(Dataset):
	def __init__(self, d_path, split, transforms):
		imgs_path = os.path.join(d_path, split)
		self.data = []
		with open(os.path.join(imgs_path, 'labels.txt'), 'r') as f:
			for line in f:
				img_name, label = line.split()
				label = int(label)
				img_path = os.path.join(imgs_path, img_name)
				self.data.append((img_path, label))
		self.transforms = transforms

	def __len__(self):
		return len(self.data)

	def __getitem__(self, index):
		img_path = self.data[index][0]
		label = self.data[index][1]
		image = Image.open(img_path)
		image = self.transforms(image)
		return image, label 

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
	for images, labels in dataloader:
		print(images.size())
		print(labels)
		exit()
	print(len(dataset))
	print(dataset[0])