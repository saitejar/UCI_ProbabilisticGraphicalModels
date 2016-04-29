import numpy as np
import pyGM as gm
import itertools
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

edges = np.genfromtxt('data/edges.txt',delimiter=None)
loc = np.genfromtxt('data/locations.txt',delimiter=None)

x,y = loc[:,0],loc[:,1]

n = loc.shape[0]
G = nx.Graph()
for node in xrange(n):
	G.add_node(node,pos=(x[node],y[node]))

for e in edges:
	G.add_edge(e[0],e[1])
pos=nx.get_node_attributes(G,'pos')
nx.draw(G,pos)
mp = {}
for i in range(n):
	mp[i]=i
nx.draw_networkx_labels(G,pos,mp)
plt.show()
