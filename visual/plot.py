import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

RG = nx.Graph()  
# RG = nx.grid_graph(dim=[4,4,4], periodic=True)
base = "input" # "input", "testdouble_8_8_128"
filename = base + ".txt"
filename_output = base + "_output_0.000001.txt"
data = pd.read_csv(filename, header=None)
data_output = pd.read_csv(filename_output, header=None)
# RG.add_node(1)
# RG.add_node(2)
# RG.add_node(3)
# RG.add_node(4)
# position={1:(1,1),2:(2,2),3:(3,4),4:(4,3)}
 
colorlist = []
position = {}
scale = 10
node_shape_list = []
for i in range(0, len(data)/2):
    RG.add_node(i)
    position[i] = (data[0][2*i], data[0][2*i+1])
    colorlist.append((data_output[0][i]*0.9)/scale)
for i in range(len(data)/2, len(data)/2+(len(data_output)-len(data)/2)/2):
    RG.add_node(i)
    position[i] = (data_output[0][len(data)/2+(i-len(data)/2)*2], data_output[0][len(data)/2+(i-len(data)/2)*2+1])
    colorlist.append(((i-len(data)/2)*0.9)/scale)
for i in range(0, len(data)/2):
    RG.add_edge(i, len(data)/2+data_output[0][i])

# print "colorlist: ", colorlist
# print "position: ", position
nx.draw(RG, node_size=30, with_labels=False, pos=position, node_color=colorlist, node_shape=">")  
import matplotlib.pyplot as plt
plt.axis('on')
plt.xlabel('X')
plt.ylabel('Y')
plt.show() 