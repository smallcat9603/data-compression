import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import sys, getopt, os

def printUsage():
    print(f'Usage: python3 {os.path.basename(__file__)} <kmeans_2d_dataset_file> <cluster_result_file>')

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "h")
    except getopt.GetoptError:
        printUsage()
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            printUsage()
            sys.exit()
        else:
            printUsage()
            sys.exit(1)
    if len(args) == 2:
        kmeans_2d_dataset_file = args[0]
        cluster_result_file = args[1] 
    else:
        printUsage()
        sys.exit(1)        

    # read kmeans dataset and clustering result
    data = pd.read_csv(kmeans_2d_dataset_file, header=None)
    means = pd.read_csv(cluster_result_file, header=None)
    
    # assign one color to one cluster
    colorlist = []
    position = {}
    scale = 100
    num_data = int(len(data)/2)
    num_means = int((len(means)-num_data)/2)
    G = nx.Graph()  
    for i in range(num_data):
        G.add_node(i)
        position[i] = (data[0][2*i], data[0][2*i+1])
        colorlist.append((means[0][i]*0.9)/scale)
    for i in range(num_data, num_data+num_means):
        G.add_node(i)
        position[i] = (means[0][num_data+(i-num_data)*2], means[0][num_data+(i-num_data)*2+1])
        colorlist.append(((i-num_data)*0.9)/scale)
    for i in range(num_data):
        G.add_edge(i, num_data+means[0][i])

    # draw or save the figure
    nx.draw(G, node_size=30, with_labels=False, pos=position, node_color=colorlist, node_shape="o")  
    plt.savefig(".".join(cluster_result_file.split(".")[:-1])+".png", dpi=256)
    #plt.show() 

if __name__ == "__main__":
   main(sys.argv[1:])  