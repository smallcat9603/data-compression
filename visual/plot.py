import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

RG = nx.Graph()  
# RG = nx.grid_graph(dim=[4,4,4], periodic=True)
data = pd.read_csv("test1.txt", header=None, sep=", ")
# RG.add_node(1)
# RG.add_node(2)
# RG.add_node(3)
# RG.add_node(4)
# position={1:(1,1),2:(2,2),3:(3,4),4:(4,3)}
position = {}
print(len(data))
for i in len(data):
    RG.add_node(i)
    position[i] = (data[0], data[1])
    i = i + 1
nx.draw(RG, node_size=30, with_labels=False, pos=position)  
import matplotlib.pyplot as plt
plt.show() 