import os
import torch
from torch.utils.data import Dataset, DataLoader
from PIL import Image
import torchvision.transforms as T
import torch.nn as nn
import torch.nn.functional as F
from model import myModel


class TestDataset(Dataset):
    def __init__(self, d_path, split, transforms):
        split='train'
        self.data = []
        imgs_path = os.path.join(d_path, split)
        with open(os.path.join(imgs_path, "labels.txt"), "r") as f:
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
        T.Normalize((0.5,), (0.5,)),
    ])
    dataset = TestDataset("data", "test", transforms)
    dataloader = DataLoader(dataset, batch_size=1024, shuffle=True, num_workers=1)
    model = myModel()
    model.load_state_dict(torch.load("myModel-triple-Relu25-0.0003"))
    model.eval()

    device = torch.device("cuda" if False else "cpu")
    with torch.no_grad():
        for image_name, images in dataloader:
            images = images.to(device)
            output = model(images)
            _, predicted_labels = torch.max(output, 1)
            for i in range(len(image_name)):
                print(image_name[i], int(predicted_labels[i]))