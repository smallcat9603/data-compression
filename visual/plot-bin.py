import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import os

RG = nx.Graph()  
# RG = nx.grid_graph(dim=[4,4,4], periodic=True)
folder = "../impl/dataset/"
base = "bit" # "bit", "bitcomp"
filename = base + ".txt"
filename_output = base + "_output.txt"

file = open(folder+filename)
if os.path.exists(folder+filename_output):
    os.remove(folder+filename_output)  
for line in file.readlines():
    line_output = line[:-1]
    for i in range(32+1-len(line.split(" "))):
        line_output += "-1 "
    line_output += "\n"
    with open(folder+filename_output, "a") as file_output:  
        file_output.write(line_output) 

data = pd.read_csv(folder+filename_output, sep=' ', header=None)

node = 0
position = {}
colorlist = []
num = 1024 #len(data)
for i in range(32):
    for j in range(num):
        if data[i][j] == 1:
           RG.add_node(node) 
           position[node] = (i, j)
           node += 1
           colorlist.append(1)
#         elif data[i][j] == 0:
#            RG.add_node(node) 
#            position[node] = (i, j)
#            node += 1
#            colorlist.append(0)
#         elif data[i][j] == -1:
#            RG.add_node(node) 
#            position[node] = (i, j)
#            node += 1
#            colorlist.append(0.5)           

plt.figure(figsize=(5, 10), dpi=64)
nx.draw(RG, node_size=10, with_labels=False, pos=position, node_color=colorlist, node_shape="s")  
plt.axis('on')
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig(filename+".png", dpi=128)
#plt.show() 