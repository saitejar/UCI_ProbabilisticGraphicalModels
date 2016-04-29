import numpy as np
import pyGM as gm
import itertools
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

# Load the data points D, and the station locations (lat/lon)
D = np.genfromtxt('data/data.txt',delimiter=None)
loc = np.genfromtxt('data/locations.txt',delimiter=None)
m,n = D.shape # m = 2760 data points, n=30 dimensional

x,y = loc[:,0],loc[:,1]

# D[i,j] = 1 if station j observed rainfall on day i
noOfFactors = (n*(n-1))/2
X = [gm.Var(i,2) for i in range(noOfFactors)]
factors = [gm.Factor([X[pair[0]],X[pair[1]]],1.0) for pair in itertools.combinations(xrange(n),2)]

for i in range(len(factors)):               # fill the tables with random values
	factors[i].table = np.zeros((2,2))

for fac in xrange(noOfFactors):
	for day in xrange(m):
		factors[fac].table[D[day,factors[fac].vars[0]],D[day,factors[fac].vars[1]]]+=1.0

for fac in xrange(noOfFactors):
	factors[fac].table /= np.sum(factors[fac].table)

def mutualInformation(fac):
	var = factors[fac].vars
	return np.sum(factors[fac].table * np.log(np.divide(factors[fac].table,np.outer(factors[fac].marginal({X[var[0]]}).table,factors[fac].marginal({X[var[1]]}).table))))
print factors[0].marginal({X[0]}).table
print type(X[factors[0].vars[0]])
MI = np.array([(0,0,0.0)]*noOfFactors)
for fac,pair in zip(range(noOfFactors),itertools.combinations(xrange(n),2)):
	v = mutualInformation(fac)

	MI[fac] = (pair[0],pair[1],-v)
print MI
G = nx.Graph()
G.add_weighted_edges_from(MI)
T = nx.minimum_spanning_tree(G)
Tree = nx.Graph()
for node in xrange(n):
	Tree.add_node(node,pos=(x[node],y[node]))
for e in T.edges(data=True):
	Tree.add_edge(e[0],e[1],weight=-e[2]['weight'])
pos=nx.get_node_attributes(Tree,'pos')
nx.draw(Tree,pos)
mp = {}
for i in range(n):
	mp[i]=i
nx.draw_networkx_labels(Tree,pos,mp)
plt.show()
