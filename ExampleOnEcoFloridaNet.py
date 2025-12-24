# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 12:30:01 2025

@author: jafish
"""

from NetworkFragilityClean import fragile_net
from RunNetworkFragilityClean import run_sparse_fn
import networkx as nx
import numpy as np
from datetime import datetime 
from NetworkDestruction import NetworkDestruction
import pymetis
#The follwing two codes (and a third necessary file) can be found at:
#https://github.com/PPNew1/Edge_Collective_Influence/tree/main
from EdgeCollectiveInfluence import *
from DualCompetitivePercolation import *
#NOTE MUST HAVE igraph installed using pip install igraph
import igraph as ig

#USEFUL (AND NECESSARY) FUNCTIONS
def nx_to_igraph(G_nx):
    """
    Convert an undirected NetworkX graph to an igraph Graph.
    Node IDs are stored in the 'name' attribute.
    """
    # Get nodes and edges
    nodes = list(G_nx.nodes())
    edges = [(nodes.index(u), nodes.index(v)) for u, v in G_nx.edges()]

    # Create igraph graph
    G_ig = Graph(n=len(nodes), edges=edges, directed=False)
    G_ig.vs["name"] = nodes  # preserve original node IDs    G_ig.vs["name"] = nodes  # preserve original node IDs
    return G_ig

def igraph_to_networkx(ig_graph):
    # Choose directed/undirected
    if ig_graph.is_directed():
        nx_graph = nx.DiGraph()
    else:
        nx_graph = nx.Graph()

    # --- Add nodes with attributes ---
    for v in ig_graph.vs:
        node_id = v.index
        attrs = dict(v.attributes())
        nx_graph.add_node(node_id, **attrs)

    # --- Add edges with attributes ---
    for e in ig_graph.es:
        source, target = e.tuple
        attrs = dict(e.attributes())
        nx_graph.add_edge(source, target, **attrs)

    return nx_graph


Today = datetime.today().strftime('%d-%m-%Y')
delta = 0.5
NetworkType = 'FloridaEcological'

#Replace the path with your own path to the file
path = 'C:/Users/jafish/Downloads/eco-florida/eco-florida.edges'
fh = open(path, "rb")
#Get the graph
G = nx.read_edgelist(fh,data=(("weight", float),))
fh.close()
#The graph has strings as node labels, which breaks things, so convert
#to integer labels
Gex = nx.convert_node_labels_to_integers(G, label_attribute="old_label")

#Start with pymetis NOTE MUST HAVE PYMETIS INSTALLED, can use 
#pip install pymetis if you have not already

print("Calculating PyMETIS partition...")
A_list1 = []
for l in range(len(Gex)):
    A_list1.append(np.array(list(Gex[l])))
    
#Find the number of edges that must be removed
NumEdges,PartLabels = pymetis.part_graph(2,A_list1)
#NumEdges represents the number of edges removed!
PyMETISFrac = NumEdges/len(Gex.edges())
FN = fragile_net()
FN.compute_Fragility(len(Gex.nodes()),len(Gex.nodes())//2,PyMETISFrac)
Frag_pymetis = FN.return_Fragility()
print("PyMETIS Fragility: ", Frag_pymetis)
print("")

#Now for the Kernighan-Lin bipartition
print("Calculating Kernighan-Lin partition...")
ND = NetworkDestruction(Graph=Gex)
C,New,S1,S2 = ND.Kernighan_Lin_NumEdges()
#C represents the number of edges removed using the Kernighan-Lin algorithm
KLFrac = C/len(Gex.edges())
FN.compute_Fragility(len(Gex.nodes()),len(Gex.nodes())//2,KLFrac)
Frag_KL = FN.return_Fragility()
print("Kernighan-Lin Fragility: ", Frag_KL)
print("")
i = 0
print("Calculating All Rewire Strategies...")
#Run for future iterative add back
Output = run_sparse_fn(Gex,delta,NetworkType+'Net_FullRemoval'+str(i))

#To save the results of this if you wish
np.savez(NetworkType+'Net_FullRemoval_' +str(i) + Today+ '.npz',Output)

#First compute Min degree strategy
G2 = nx.Graph(Output['MinDegree']['Graphs'][-1])
MinFrac = (len(Gex.edges())-len(G2.edges()))/len(Gex.edges())
FN = fragile_net()
FN.compute_Fragility(len(Gex.nodes()),len(Gex.nodes())//2,MinFrac)
Frag_Min = FN.return_Fragility()
print("Min-Degree Fragility: ",Frag_Min)
print("")
FN.add_graph(Gex)
FN.add_graph_new(G2)
G3 = FN.sparse_iterative_add_back(len(Gex.nodes())//2)
MinDegRewFrac = (len(Gex.edges())-len(G3.edges()))/len(Gex.edges())
FN.compute_Fragility(len(Gex.nodes()),len(Gex.nodes())//2,MinDegRewFrac)
Frag_MinRew = FN.return_Fragility()
#Now for edge betweenness
G2 = nx.Graph(Output['EdgeBetweenness']['Graphs'][-1])
EBFrac = (len(Gex.edges())-len(G2.edges()))/len(Gex.edges())
FN = fragile_net()
FN.compute_Fragility(len(Gex.nodes()),len(Gex.nodes())//2,EBFrac)
Frag_EB = FN.return_Fragility()
print("Edge-Betweenness Fragility: ",Frag_EB)
print("")
FN.add_graph(Gex)
FN.add_graph_new(G2)
G3 = FN.sparse_iterative_add_back(len(Gex.nodes())//2)
EBRewFrac = (len(Gex.edges())-len(G3.edges()))/len(Gex.edges())
FN.compute_Fragility(len(Gex.nodes()),len(Gex.nodes())//2,EBFrac)
Frag_EBRew = FN.return_Fragility()
#Now for edge sum
G2 = nx.Graph(Output['EdgeSum']['Graphs'][-1])
ESFrac = (len(Gex.edges())-len(G2.edges()))/len(Gex.edges())
FN = fragile_net()
FN.compute_Fragility(len(Gex.nodes()),len(Gex.nodes())//2,MinFrac)
Frag_ES = FN.return_Fragility()
print("Edge-Sum Fragility: ",Frag_ES)
print("")
FN.add_graph(Gex)
FN.add_graph_new(G2)
G3 = FN.sparse_iterative_add_back(len(Gex.nodes())//2)
ESRewFrac = (len(Gex.edges())-len(G3.edges()))/len(Gex.edges())
FN.compute_Fragility(len(Gex.nodes()),len(Gex.nodes())//2,ESRewFrac)
Frag_ESRew = FN.return_Fragility()

Frag_Rew = np.max([Frag_MinRew,Frag_EBRew, Frag_ESRew])
print("Rewire Fragility: ",Frag_Rew )
print("")

#Now run for ECI (edge collective influence)
g = nx_to_igraph(Gex)
res = IECIR(g,p=0.5)
NumEdges = len(res[2][2])
ECIFrac = NumEdges/len(Gex.edges())
FN = fragile_net()
FN.compute_Fragility(len(Gex.nodes()),len(Gex.nodes())//2,ECIFrac)
Frag_ECI = FN.return_Fragility()
print("ECI Fragility: ", Frag_ECI)
print("")

print("Full Results: ")
print("")
print("PyMETIS Fragility: ", Frag_pymetis)
print("")
print("Kernighan-Lin Fragility: ", Frag_KL)
print("")
print("Min-Degree Fragility: ",Frag_Min)
print("")
print("Edge-Betweenness Fragility: ",Frag_EB)
print("")
print("Edge-Sum Fragility: ",Frag_ES)
print("")
print("Rewire Fragility: ",Frag_Rew )
print("")
print("ECI Fragility: ", Frag_ECI)
print("")