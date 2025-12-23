#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 12:20:07 2022

@author: jeremie
"""

from NetworkFragilityClean import fragile_net
import numpy as np
import networkx as nx
import copy as COPY


def run_sparse_fn(G,Delta,NetworkType=''):
    Output = {} #There is going to be a lot of output so get ready!
    NestedOutput = {}
    NestedOutput2 = {}
    NestedOutput3 = {}
    FN = fragile_net()
    FN.add_graph(G)
    
    CG = COPY.deepcopy(G)
    #Start the first method...
    Method = 'MinDegree'
    print('Starting ' + Method + ' Now!')
    
    Message = 'Starting: ' + NetworkType + ' Runs Now!'
    print(Message)
    fracs = []
    removals = []
    fracwithoutrewire = []
    min_frac = 1
    numRemovals = 1
    Glist = [nx.adjacency_matrix(G)]
    Glist2 = [nx.adjacency_matrix(G)]
    while min_frac>Delta:
        print('Number of Removals: ',numRemovals)
        removals.append(numRemovals)
        #print('This is the list of removals: ', removals)
        Gnew,RL = FN.sparse_gg_removal(1,AttackStrategy=Method)
        AAANew = nx.adjacency_matrix(Gnew)
        Glist.append(AAANew)
        #Anew1 = COPY.deepcopy(Anew)
        Size,frac = FN.size_LCC(Gnew,return_fraction=True)
        
        fracwithoutrewire.append(frac)
        print(frac)
        SUM = 1
        while SUM>0:
            Gnew2 = FN.sparse_rwtf_greedy()
            SUM = FN.size_LCC(Gnew2)-FN.size_LCC(Gnew)
        
        Size,frac = FN.size_LCC(Gnew2,return_fraction=True)
        min_frac = frac
        print(frac)
        fracs.append(frac)
        Glist2.append(nx.adjacency_matrix(Gnew2))

        numRemovals = numRemovals+1
        #Reset Graph to previous one without rewiring...
        FN.add_graph_new(None)
        FN.add_graph(nx.Graph(Glist[len(Glist)-1]))
    
    #print("This is removals: ", removals)
    NestedOutput['Removals'] = np.array(removals)
    NestedOutput['Fractions'] = np.array(fracs)
    NestedOutput['FinalAdjacency'] = Gnew2
    NestedOutput['FinalNumRemoved'] = numRemovals-1
    NestedOutput['FractionsWithoutRewire'] = np.array(fracwithoutrewire)
    NestedOutput['GraphsWORW'] = Glist
    NestedOutput['Graphs'] = Glist2
    Output[Method] = NestedOutput
    
    #################################################################

    #################################################################    
    Method = 'EdgeBetweenness'
    print('Starting ' + Method + ' Now!')
    
    #Reset for calculational purposes
    FN.add_graph_new(None)
    FN.add_graph(CG)
    
    Message = 'Starting: ' + Method + ' Runs Now!'
    print(Message)
    fracs = []
    removals = []
    fracwithoutrewire = []
    min_frac = 1
    numRemovals = 1
    Glist = [nx.adjacency_matrix(CG)]
    Glist2 = [nx.adjacency_matrix(CG)]
    while min_frac>Delta:
        print('Number of Removals: ',numRemovals)
        removals.append(numRemovals)
        Gnew,RL = FN.sparse_gg_removal(1,Method)
        AAANew = nx.adjacency_matrix(Gnew)
        Glist.append(AAANew)
        #Glist.append(Gnew)
        Size,frac = FN.size_LCC(Gnew,return_fraction=True)
        
        fracwithoutrewire.append(frac)
        SUM=1
        while SUM>0:
            Gnew2 = FN.sparse_rwtf_greedy()
            SUM = FN.size_LCC(Gnew2)-FN.size_LCC(Gnew)
        Size,frac = FN.size_LCC(Gnew2,return_fraction=True)
        min_frac = frac
        fracs.append(frac)
        Glist2.append(nx.adjacency_matrix(Gnew2))
        numRemovals = numRemovals+1
        #Reset Graph to previous one without rewiring...
        FN.add_graph_new(None)
        FN.add_graph(nx.Graph(Glist[len(Glist)-1]))
        
    NestedOutput2['Removals'] = np.array(removals)
    NestedOutput2['Fractions'] = np.array(fracs)
    NestedOutput2['FinalAdjacency'] = Gnew2
    NestedOutput2['FinalNumRemoved'] = numRemovals-1
    NestedOutput2['FractionsWithoutRewire'] = np.array(fracwithoutrewire)
    NestedOutput2['GraphsWORW'] = Glist
    NestedOutput2['Graphs'] = Glist2
    Output[Method] = NestedOutput2

    Method = 'EdgeSum'
    print('Starting ' + Method + ' Now!')
    
    #Reset for calculational purposes
    FN.add_graph_new(None)
    FN.add_graph(CG)
    
    Message = 'Starting: ' + Method + ' Runs Now!'
    print(Message)
    fracs = []
    removals = []
    fracwithoutrewire = []
    min_frac = 1
    numRemovals = 1
    Glist = [nx.adjacency_matrix(CG)]
    Glist2 = [nx.adjacency_matrix(CG)]
    while min_frac>Delta:
        print('Number of Removals: ',numRemovals)
        removals.append(numRemovals)
        Gnew,RL = FN.sparse_gg_removal(1,Method)
        AAANew = nx.adjacency_matrix(Gnew)
        Glist.append(AAANew)
        #Glist.append(Gnew)
        Size,frac = FN.size_LCC(Gnew,return_fraction=True)
        
        fracwithoutrewire.append(frac)
        SUM=1
        while SUM>0:
            Gnew2 = FN.sparse_rwtf_greedy()
            SUM = FN.size_LCC(Gnew2)-FN.size_LCC(Gnew)
        Size,frac = FN.size_LCC(Gnew2,return_fraction=True)
        min_frac = frac
        fracs.append(frac)
        Glist2.append(nx.adjacency_matrix(Gnew2))
        numRemovals = numRemovals+1
        #Reset Graph to previous one without rewiring...
        FN.add_graph_new(None)
        FN.add_graph(nx.Graph(Glist[len(Glist)-1]))
        
    NestedOutput3['Removals'] = np.array(removals)
    NestedOutput3['Fractions'] = np.array(fracs)
    NestedOutput3['FinalAdjacency'] = Gnew2
    NestedOutput3['FinalNumRemoved'] = numRemovals-1
    NestedOutput3['FractionsWithoutRewire'] = np.array(fracwithoutrewire)
    NestedOutput3['GraphsWORW'] = Glist
    NestedOutput3['Graphs'] = Glist2
    Output[Method] = NestedOutput3
    
    return Output



def run_sparse_worw(G,Delta,NetworkType=''):
    Output = {} #There is going to be a lot of output so get ready!
    NestedOutput = {}
    NestedOutput2 = {}
    NestedOutput3 = {}
    FN = fragile_net()
    GG = COPY.deepcopy(G)
    FN.add_graph(GG)
    
    
    #Start the first method...
    Method = 'MinDegree'
    print('Starting ' + Method + ' Now!')
    
    Message = 'Starting: ' + NetworkType + ' Runs Now!'
    print(Message)
    fracs = []
    removals = []
    fracwithoutrewire = []
    Glist = [nx.adjacency_matrix(G)]
    min_frac = 1
    numRemovals = 1
    while min_frac>Delta:
        print('Number of Removals: ',numRemovals)
        removals.append(numRemovals)
        
        Gnew,L = FN.sparse_gg_removal(1,AttackStrategy=Method)
        AAANew = nx.adjacency_matrix(Gnew)
        Glist.append(AAANew)
       
        Size,frac = FN.size_LCC(Gnew,return_fraction=True)
        print(frac)
        fracwithoutrewire.append(frac)
        #FN.add_graph(Gnew)
        min_frac = frac
        

        numRemovals = numRemovals+1
    
    #print("This is removals: ", removals)
    NestedOutput['Removals'] = np.array(removals)
    NestedOutput['Fractions'] = np.array(fracs)
    NestedOutput['FinalAdjacency'] = Gnew
    NestedOutput['FinalNumRemoved'] = numRemovals-1
    NestedOutput['FractionsWithoutRewire'] = np.array(fracwithoutrewire)
    NestedOutput['Graphs'] = Glist
    Output[Method] = NestedOutput
    
    #################################################################

    #################################################################    
    Method = 'EdgeBetweenness'
    print('Starting ' + Method + ' Now!')
    
    #Reset for calculational purposes
    
    FN.add_graph_new(None)
    
    Message = 'Starting: ' + NetworkType + ' Runs Now!'
    print(Message)
    fracs = []
    removals = []
    fracwithoutrewire = []
    min_frac = 1
    Glist = [nx.adjacency_matrix(G)]
    
    numRemovals = 1
    while min_frac>Delta:
        print('Number of Removals: ',numRemovals)
        removals.append(numRemovals)
        Gnew,L = FN.sparse_gg_removal(1,AttackStrategy=Method)
        AAANew = nx.adjacency_matrix(Gnew)
        Glist.append(AAANew)
        
        Size,frac = FN.size_LCC(Gnew,return_fraction=True)
        print(frac)
        fracwithoutrewire.append(frac)
        
        min_frac =  frac
        numRemovals = numRemovals+1
        
        
    NestedOutput2['Removals'] = np.array(removals)
    NestedOutput2['Fractions'] = np.array(fracs)
    NestedOutput2['FinalAdjacency'] = Gnew
    NestedOutput2['FinalNumRemoved'] = numRemovals-1
    NestedOutput2['FractionsWithoutRewire'] = np.array(fracwithoutrewire)
    NestedOutput2['Graphs'] = Glist
    Output[Method] = NestedOutput2

    #################################################################    
    Method = 'EdgeSum'
    print('Starting ' + Method + ' Now!')
    
    #Reset for calculational purposes
    
    FN.add_graph_new(None)
    
    Message = 'Starting: ' + NetworkType + ' Runs Now!'
    print(Message)
    fracs = []
    removals = []
    fracwithoutrewire = []
    min_frac = 1
    Glist = [nx.adjacency_matrix(G)]
    
    numRemovals = 1
    while min_frac>Delta:
        print('Number of Removals: ',numRemovals)
        removals.append(numRemovals)
        Gnew,L = FN.sparse_gg_removal(1,AttackStrategy=Method)
        AAANew = nx.adjacency_matrix(Gnew)
        Glist.append(AAANew)
        
        Size,frac = FN.size_LCC(Gnew,return_fraction=True)
        print(frac)
        fracwithoutrewire.append(frac)
        
        min_frac =  frac
        numRemovals = numRemovals+1
        
        
    NestedOutput3['Removals'] = np.array(removals)
    NestedOutput3['Fractions'] = np.array(fracs)
    NestedOutput3['FinalAdjacency'] = Gnew
    NestedOutput3['FinalNumRemoved'] = numRemovals-1
    NestedOutput3['FractionsWithoutRewire'] = np.array(fracwithoutrewire)
    NestedOutput3['Graphs'] = Glist
    Output[Method] = NestedOutput3
    
    return Output


