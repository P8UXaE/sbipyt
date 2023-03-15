import torch

# Define the neural network model
class MyModel(torch.nn.Module):
    def __init__(self):
        super(MyModel, self).__init__()
        self.linear = torch.nn.Linear(2, 1)
        
    def forward(self, x1, x2):
        x = torch.cat((x1, x2), dim=1)
        out = self.linear(x)
        return out

# Define the input matrices
x1 = torch.randn(1, 2)
x2 = torch.randn(1, 2)

# Create an instance of the neural network model
model = MyModel()

# Pass the input matrices through the model to get the output
output = model(x1, x2)

# Print the output
print(output.item())
