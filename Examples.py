# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 12:43:51 2022

@author: jefis
"""

from NetworkFragilityClean import fragile_net
from RunNetworkFragilityClean import run_sparse_fn,run_sparse_worw
import scipy.io as IO
import networkx as nx
import numpy as np
from datetime import datetime
import copy as COPY

Today = datetime.today().strftime('%d-%m-%Y')

##############################################################################
##############################################################################
################# Start by destroying the networks  ##############
#                 then we will estimate the fragility
##############################################################################
##############################################################################

#Setup the net
FN = fragile_net()

#Number of trials
NumTrials = 1 #This will be ignored if the network is the mall network...

#Network parameters
n = 100

## Set delta (between 1/n and (n-1)/n) ##
delta = 0.5 #This corresponds to a paramter c = 50
c = 50 #Since this represents delta = 0.5

#Barabasi-Albert parameter
m = 4
#Note the this is the total number of edges of Barabasi-Albert graph
# with n =100, m = 4
total_edges = 384

#Watts-Strogatz parameters
k = 8
p = 0.25

Types = ['BA', 'ER', 'WS', 'NorthMall']
NetworkType = Types[2]

if NetworkType =='BA':
    G = nx.barabasi_albert_graph(n,m)
    
elif NetworkType=='ER':
    A = FN.exact_erdos_renyi_graph(n,total_edges)
    G = nx.Graph(A)
    
elif NetworkType =='WS':
    G = nx.watts_strogatz_graph(n,k,p)
    Diff = len(G.edges())-total_edges
    for i in range(Diff):
        Arr = np.array(list(G.edges()))
        np.random.shuffle(Arr)
        G.remove_edge(Arr[0,0],Arr[0,1])
elif NetworkType =='NorthMall':
    LM = IO.loadmat('Amatrix_northmall.mat')
    A = LM['A1']
    G = nx.Graph(A)   


if NetworkType == 'NorthMall':
    i = 0
    Fragilities = np.zeros(1) #Set the fragilities array
    Output = run_sparse_fn(G,delta,NetworkType+'Net_FullRemoval'+str(i))#
    np.savez(NetworkType+'Net_FullRemoval_' +str(i) + '_Edge' + Today+ '.npz',Output)
    
    ##############################################################
    ##       Fragility estimation step                          ##
    ##############################################################
    FN = fragile_net()
    
    #We must check both the minimum degree attack strategy and the edge betweenness
    #strategy because for different networks and even different values of delta
    #one is typically better than the other.
    MDDict = Output['MinDegree']
    #Grab the final graph after the edges have been removed
    G2 = nx.Graph(MDDict['Graphs'][-1])
    G3 = COPY.deepcopy(G2)
    
    
    #Add the original graph and the final graph to compute fragility
    FN.add_graph(G)
    FN.add_graph_new(G2)
    
    #Complete the final step of the algorithm (so far we have only done the 
    #first two steps not the part for adding edges back which do not affect the
    #LCC) c was defined above based on delta
    G2 = FN.sparse_iterative_add_back(c)
    #Fraction of edges which have been removed from G
    MDFrac = (len(G.edges())-len(G2.edges()))/len(G.edges())
    FN.compute_Fragility(len(G.nodes()),c,MDFrac)
    
    #The estimate for the minimum degree attack strategy
    MDEstimatedFrag = FN.return_Fragility()
    
    #Compute minimum degree without rewiring fragility
    MDFrac2 = (len(G.edges())-len(G3.edges()))/len(G.edges())
    FN.compute_Fragility(len(G.nodes()),c,MDFrac2)
    OthMD = FN.return_Fragility()
    
    #Reset FN to eliminate graphs
    FN = fragile_net()
    
    EBDict = Output['EdgeBetweenness']
    #Grab the final graph after the edges have been removed
    G2 = nx.Graph(EBDict['Graphs'][-1])
    G4 = COPY.deepcopy(G2)
    
    #Add the original graph and the final graph to compute fragility
    FN.add_graph(G)
    FN.add_graph_new(G2)
    
    #Complete the final step of the algorithm (so far we have only done the 
    #first two steps not the part for adding edges back which do not affect the
    #LCC) c was defined above based on delta
    G2 = FN.sparse_iterative_add_back(c)
    #Fraction of edges which have been removed from G
    EBFrac = (len(G.edges())-len(G2.edges()))/len(G.edges())
    FN.compute_Fragility(len(G.nodes()),c,EBFrac)
    #The estimate for the edge betweenness attack strategy
    EBEstimatedFrag = FN.return_Fragility()
    
    #Compute edge betweenness without rewiring fragility
    EBFrac2 = (len(G.edges())-len(G4.edges()))/len(G.edges())
    FN.compute_Fragility(len(G.nodes()),c,EBFrac2)
    OthEB = FN.return_Fragility()
    
    #Now we can get the true estimated fragility
    Fragility = np.max(np.array([MDEstimatedFrag,EBEstimatedFrag]))
    Fragilities[i] = Fragility
    print("This is the estimated fragility: ", Fragility)
    
    #NOTE THAT THIS IS HOW FRAGILITY WAS ESTIMATED, HOWEVER ONE COULD LOOK AT
    #THE SITUATION WITHOUT REWIRING, THIS IS DONE BELOW
    #If desired this can be compared to the results without rewiring as part
    #of the process...
    Output2 = run_sparse_worw(G,delta,NetworkType+'Net_WORW_FullRemoval'+str(i))
    np.savez(NetworkType+'Net_WORW_FullRemoval_' +str(i) + '_Edge' + Today+ '.npz',Output2)


else:
    Fragilities = np.zeros(NumTrials)
    for i in range(NumTrials):
        
        ######################################################################
        ######################################################################
        #Make a new G at every step...
        if NetworkType =='BA':
            G = nx.barabasi_albert_graph(n,m)
    
        elif NetworkType=='ER':
            A = FN.exact_erdos_renyi_graph(n,total_edges)
            G = nx.Graph(A)
    
        elif NetworkType =='WS':
            G = nx.watts_strogatz_graph(n,k,p)
            Diff = len(G.edges())-total_edges
            for j in range(Diff):
                Arr = np.array(list(G.edges()))
                np.random.shuffle(Arr)
                G.remove_edge(Arr[0,0],Arr[0,1])
        
        ###########################################################
        Output = run_sparse_fn(G,delta,NetworkType+'Net_FullRemoval'+str(i))#
        np.savez(NetworkType+'Net_FullRemoval_' +str(i) + '_Edge' + Today+ '.npz',Output)


        ##############################################################
        ##       Fragility estimation step                          ##
        ##############################################################
        FN = fragile_net()
        
        #We must check both the minimum degree attack strategy and the edge betweenness
        #strategy because for different networks and even different values of delta
        #one is typically better than the other.
        MDDict = Output['MinDegree']
        #Grab the final graph after the edges have been removed
        G2 = nx.Graph(MDDict['Graphs'][-1])
        G3 = COPY.deepcopy(G2)
        
        #Add the original graph and the final graph to compute fragility
        FN.add_graph(G)
        FN.add_graph_new(G2)
        
        #Complete the final step of the algorithm (so far we have only done the 
        #first two steps not the part for adding edges back which do not affect the
        #LCC) c was defined above based on delta
        G2 = FN.sparse_iterative_add_back(c)
        #Fraction of edges which have been removed from G
        MDFrac = (len(G.edges())-len(G2.edges()))/len(G.edges())
        FN.compute_Fragility(len(G.nodes()),c,MDFrac)
        #The estimate for the minimum degree attack strategy
        MDEstimatedFrag = FN.return_Fragility()
        
        #Compute minimum degree without rewiring fragility
        MDFrac2 = (len(G.edges())-len(G3.edges()))/len(G.edges())
        FN.compute_Fragility(len(G.nodes()),c,MDFrac2)
        OthMD = FN.return_Fragility() #Fragility without rewiring
        
        #Reset FN to eliminate graphs
        FN = fragile_net()
        
        EBDict = Output['EdgeBetweenness']
        #Grab the final graph after the edges have been removed
        G2 = nx.Graph(EBDict['Graphs'][-1])
        G4 = COPY.deepcopy(G2)
        
        #Add the original graph and the final graph to compute fragility
        FN.add_graph(G)
        FN.add_graph_new(G2)
        
        #Complete the final step of the algorithm (so far we have only done the 
        #first two steps not the part for adding edges back which do not affect the
        #LCC) c was defined above based on delta
        G2 = FN.sparse_iterative_add_back(c)
        #Fraction of edges which have been removed from G
        EBFrac = (len(G.edges())-len(G2.edges()))/len(G.edges())
        FN.compute_Fragility(len(G.nodes()),c,EBFrac)
        #The estimate for the edge betweenness attack strategy
        EBEstimatedFrag = FN.return_Fragility()
        
        #Compute edge betweenness without rewiring fragility
        EBFrac2 = (len(G.edges())-len(G4.edges()))/len(G.edges())
        FN.compute_Fragility(len(G.nodes()),c,EBFrac2)
        OthEB = FN.return_Fragility() #Fragility without rewiring
        
        #Now we can get the true estimated fragility
        Fragility = np.max(np.array([MDEstimatedFrag,EBEstimatedFrag]))
        print("This is the estimated fragility: ", Fragility)
        Fragilities[i] = Fragility


    
        Output2 = run_sparse_worw(G,delta,NetworkType+'Net_WORW_FullRemoval'+str(i))
        np.savez(NetworkType+'Net_WORW_FullRemoval_' +str(i) + '_Edge' + Today+ '.npz',Output2)
    

##############################################################################
##############################################################################
        