#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 10:16:35 2022

@author: jeremie
"""
import numpy as np
import networkx as nx
import copy as COPY


class fragile_net:
    """A class for determining network fragility based on the paper by
       Fish, Banavar and Bollt..."""
    def __init__(self):
        self.A = np.array([0])
        self.Aold = np.array([0])
        self.Anew = np.array([0])
        self.Graph = None
        self.CompFrac = None
        self.Anew2 = np.array([0])
        self.maxLength = None
        self.Gnew = None
        self.Fragility=None

    
    def return_adjacency_matrix(self):
        return self.A
    
    def return_old_adjacency_matrix(self):
        return self.Aold
    
    def return_graph(self):
        return self.Graph
    
    def return_size_LCC(self):
        return self.maxLength
    
    def return_new_adjacency_matrix(self):
        return self.Anew
    
    def return_complete_graph_fraction(self):
        return self.CompFrac
    
    def return_Fragility(self):
        return self.Fragility
    

    def complete_graph_fraction(self,n,c):
        """Assumption is that the graph is undirected
           Input: n - (pos. int.) The number of nodes in the network
                  c - (pos. int. with c<n) The number of nodes in the LCC"""
        m = n*(n-1)/2
        b = np.mod(n,n-c)
        flr = np.floor(n/(c))
        Kc = (c)*(c-1)/2
        Kb = b*(b-1)/2
        numerator = m - (flr*Kc+Kb)
        self.CompFrac = numerator/m
    

        
    def compute_Fragility(self,n,c,TopFrac):
        """Compute the estimated fragility of the network
           Input: n - (pos. int.)The number of nodes of the network (and corresponding 
                      complete network) 
                  c - (pos. int with c< n) The size of the largest connected 
                      component.
           """
        self.complete_graph_fraction(n,c)
        CompFraction = self.return_complete_graph_fraction()
        
        self.Fragility = 1-(TopFrac/CompFraction)
        
    def add_adjacency_matrix(self,A):
        """Add the adjacency matrix, this will be used for most things
           in the fragile_net class. However note that some may require
           Aold and Anew rather than just A...This is used only to add
           an external adjacency matrix (there are methods for creating
                                         some built in)."""
        self.A = A
    
    def add_graph_new(self,G):
        """Method for addding a new graph"""
        self.Gnew = G
    
    def add_graph_old(self,G):
        """Method for saving a an original copy of the graph"""
        self.Gold = COPY.deepcopy(G)
    def add_graph(self,G):
        """Add the graph to the class object. This is needed for some things 
           in the class, so it is not always necessary to add this. This is meant
           so you can provide your own if you so choose...
           """
        self.Graph = G
    
    def add_old_adjacency_matrix(self,Aold):
        """Method to add the old adjacency matrix"""
        self.Aold = Aold
    
    def add_new_adjacency_matrix(self,Anew):
        """Method to add the new adjacency matrix"""
        self.Anew = Anew
    

    def sparse_gg_removal(self,numRemovals,AttackStrategy="EdgeBetweenness",AttackLCC=True):
        """Performs gready removal on a networkx Graph
           -numRemovals - (pos. int.)The number of edges to remove 
            AttackStrategy - (string) The method of attack, currently available options are MinDegree to attack the 
                             minimum degree nodes and EdgeBetweenness to attack 
                             edges with the largest value of edge betweenness
            AttackLCC - (Bool) Set to True if the attack should be performed on
                        the largest connected component
            
            returns Gnew- the new networkx graph with edges removed base on AttackStrategy
                    RemovalList - The list of edges removed from the graph in 
                                  the order they were removed."""
        Gold = COPY.deepcopy(self.Graph)
        if self.Gnew ==None:
            G = COPY.deepcopy(Gold)
        else:
            G = self.Gnew
        RemovalList = []
        #G = nx.Graph(Anew)
        if AttackStrategy=='EdgeBetweenness':
            if AttackLCC:
                LCCcomps = self.components_LCC(G)
                print("This is LCC: ",LCCcomps)
                S = G.subgraph(LCCcomps)
                for i in range(numRemovals):
                    print(i)
                    Dict = nx.edge_betweenness_centrality(S)
                    L = list(sorted(Dict.items(), key=lambda item: item[1]))
                    L = L[len(L)-1]
                    G.remove_edge(*L[0])
                    RemovalList.append(L[0])
                    LCCcomps = self.components_LCC(G)
                    #print("This is LCC: ",LCCcomps)
                    S = G.subgraph(LCCcomps)
            else:
                for i in range(numRemovals):
                    print(i)
                    Dict = nx.edge_betweenness_centrality(G)
                    L = list(sorted(Dict.items(), key=lambda item: item[1]))
                    L = L[len(L)-1]
                    G.remove_edge(*L[0])
                    RemovalList.append(L[0])
            
            
        elif AttackStrategy=='MinDegree':
            if AttackLCC:
                LCCcomps = self.components_LCC(G)
                S = G.subgraph(LCCcomps)
                Start = 0
                Num=0
                
                for i in range(numRemovals):
                    print(i)
                    Degs = np.array([val for (node, val) in S.degree()])
                    Inds = np.argsort(Degs)
                    Nodes = list(S.nodes())
                    
                    Attack = Nodes[Inds[Start]]
                   
                            
                    Edges = np.array([val for (node, val) in S.edges(Attack)])
                    G.remove_edge(Attack,Edges[0])
                    LCCcomps = self.components_LCC(G)
                    S = G.subgraph(LCCcomps)
            else:
                Start = 0
                Num=0
                #print("Hello I am here!")
                for i in range(numRemovals):
                    Degs = np.array([val for (node, val) in G.degree()])
                    Inds = np.argsort(Degs)
                    #print(Degs[Inds[Start]])
                    if Degs[Inds[Start]]>0:
                        Attack = Inds[Start]
                    else:
                        while Num<1:
                            Start=Start+1
                            Attack = Inds[Start]
                            Num = Degs[Attack]
                            
                    Edges = np.array([val for (node, val) in G.edges(Attack)])
                    G.remove_edge(Attack,Edges[0])
                    RemovalList.append(Edges[0])
                
     
                    
        else:
            raise ValueError('Only EdgeBetweenness and MinDegree Currently Allowed for AttackStrategy')
        
        self.Gnew = G
        return G,RemovalList    
    
    def greedy_global_removal(self,numRemovals,AttackStrategy="EdgeBetweenness",AttackLCC=True):
        """
           numRemovals    - (pos. int.) The number of edges to remove.
           AttackStrategy - (String) Currently only Edge Betweenness is supported
                                     Min Degree and Max degree strategies are
                                     allowed.
           AttackLCC      - (default True) If true then the LCC will always be
                            attacked for the greedy portion of the algorithm,
                            otherwise the greedy algorithm will operate on the 
                            entire network.
         returns Anew     - The new adjacency matrix after removals
                 RemovalList - The list of edge removals. Currently only implemented
                               for EdgeBetweenness but is on todo list for others.
                                     """
        A = self.A
        if len(A)==1:
            raise ValueError("Missing adjacency matrix, must be added!")
        Anew = COPY.deepcopy(A)
        RemovalList = []
        G = nx.Graph(Anew)
        if AttackStrategy=='EdgeBetweenness':
            if AttackLCC:
                LCCcomps = self.components_LCC(G)
                S = G.subgraph(LCCcomps)
                for i in range(numRemovals):
                    print(i)
                    Dict = nx.edge_betweenness_centrality(S)
                    L = list(sorted(Dict.items(), key=lambda item: item[1]))
                    L = L[len(L)-1]
                    G.remove_edge(*L[0])
                    RemovalList.append(L[0])
                    LCCcomps = self.components_LCC(G)
                    S = G.subgraph(LCCcomps)
            else:
                for i in range(numRemovals):
                    print(i)
                    Dict = nx.edge_betweenness_centrality(G)
                    L = list(sorted(Dict.items(), key=lambda item: item[1]))
                    L = L[len(L)-1]
                    G.remove_edge(*L[0])
                    RemovalList.append(L[0])
            
            
        elif AttackStrategy=='MinDegree':
            if AttackLCC:
                LCCcomps = self.components_LCC(G)
                S = G.subgraph(LCCcomps)
                Start = 0
                Num=0
                
                for i in range(numRemovals):
                    Degs = np.array([val for (node, val) in S.degree()])
                    Inds = np.argsort(Degs)
                    Nodes = list(S.nodes())
                    
                    Attack = Nodes[Inds[Start]]
                   
                            
                    Edges = np.array([val for (node, val) in S.edges(Attack)])
                    G.remove_edge(Attack,Edges[0])
                    LCCcomps = self.components_LCC(G)
                    S = G.subgraph(LCCcomps)
            else:
                Start = 0
                Num=0
                #print("Hello I am here!")
                for i in range(numRemovals):
                    Degs = np.array([val for (node, val) in G.degree()])
                    Inds = np.argsort(Degs)
                    #print(Degs[Inds[Start]])
                    if Degs[Inds[Start]]>0:
                        Attack = Inds[Start]
                    else:
                        while Num<1:
                            Start=Start+1
                            Attack = Inds[Start]
                            Num = Degs[Attack]
                            
                    Edges = np.array([val for (node, val) in G.edges(Attack)])
                    G.remove_edge(Attack,Edges[0])
                    RemovalList.append(Edges[0])
                
       
                    
        else:
            raise ValueError('Only EdgeBetweenness and MinDegree Currently Allowed for AttackStrategy')
        Anew = nx.adjacency_matrix(G)
        self.Anew = Anew.todense()
        #print("This is the LCC stuff: ", self.size_LCC(nx.Graph(Anew)))
        return self.Anew, RemovalList    
    def components_LCC(self,G):
        """Method to find the components of the largest connected component of 
        the graph. 
        Input: G- The graph whose largest connected component should be found
        
        Output: maxList - The list of nodes in the largest connected component"""
        self.G = G
        
    
        lst = list(nx.connected_components(G))
        maxList = max(lst, key = lambda i: len(i)) 
        
        return maxList
    def size_LCC(self,G,return_fraction=False):
        """G               - A networkx Graph.
           return_fraction - If True then return the fraction of the LCC in addition
                             to the size of the LCC"""
        
        self.G = G
        if return_fraction == False:
            lst = list(nx.connected_components(G))
            maxList = max(lst, key = lambda i: len(i)) 
            
            self.maxLength = len(maxList) 
            return self.maxLength
        else:
            lst = list(nx.connected_components(G))
            maxList = max(lst, key = lambda i: len(i)) 
            self.maxLength = len(maxList)
            
            return self.maxLength,self.maxLength/len(G.nodes)
        
    def bottle_neck_graph(self,n,p):
        """n - The number of nodes in the network, should be even. If not it will be
               changed to be even... it is assumed that n>=16
           p - The Erdos-Renyi connection probability.
           
           returns A - The adjacency matrix..
           
           Note this network has a lot of symmetry, on both ends of the fragile edge
           are hubs which are connected to identical ER graphs and then in the center
           are two nodes which are connected to the same set of nodes on their respective
           sides of the ER graphs..."""
        if np.floor(n/2)!=n/2:
            print("Warning Changing n to: ",n," please remember to set n to an even number!")
            n = n+1
        if n<16:
            raise ValueError("n must be greater than or equal to 8!")
        A = np.zeros((n,n))
        #split into 3 sets of nodes, one set contains two highly connected nodes
        #another set contains the nodes across which the fragile edge will go
        #and a third is the Erdos-Renyi graph which will be the same on both sides
        #of the fragile edge.
        k = n-4
        k = np.int32(k/2)
        
        #Find the random connections for the central nodes (across which the fragile
        #edge will be formed).
        NumConnections =np.int32(np.max([1,np.floor(0.25*k)]))
        Rng = np.arange(1,k+1)
        np.random.shuffle(Rng)
        RandConnections = Rng[0:NumConnections]
        
        #Start creating the adjacency matrix combining above elements
        A[0,1:k+1] = np.ones(k)
        A[1:k+1,0] = np.ones(k)
        A[n-1,-1-k:-1] = np.ones(k)
        A[-1-k:-1,n-1] = np.ones(k)
        #Erdos-Renyi graph...
        R = np.random.rand(k,k)
        R = np.triu(R,1)
        R = R+R.T
        R[R<=p] = 1
        R[R<1] = 0
        #Make ER graph on BOTH sides!
        A[1:k+1,1:k+1] = R
        A[-1-k:-1,-1-k:-1] = R
        #Now add in the central nodes
        A[k+1,RandConnections] = 1
        A[RandConnections,k+1] = 1
        A[-k-2,-1-RandConnections] = 1
        A[-1-RandConnections,-k-2] = 1
        #Now make the one fragile connection...
        A[k+1,-k-2]=1
        A[-k-2,k+1]=1
        np.fill_diagonal(A,0)
        self.A = np.int32(A)
        return self.A
    
    def sparse_rwtf_greedy(self):
        """Method to perform rewiring on a networkx graph (if available). This 
        method is performed after the greedy removal process has been performed
        
        Output - The new graph after rewiring (may be the same as the original)"""
        Gold = COPY.deepcopy(self.Graph)
        Gnew = COPY.deepcopy(self.Gnew)
        M = self.size_LCC(Gnew)
        
        if M == len(Gold.nodes):
            return Gnew
        else:
            #Find the nodes that live in the LCC
            M = self.components_LCC(Gold)#max(nx.connected_components(Gold))
            M2 = self.components_LCC(Gnew)#max(nx.connected_components(Gnew))
            
            #Find the node which could potentially be rewired
            Degs = np.array([val for (node, val) in Gold.degree(M2)])
            Degs2 = np.array([val for (node, val) in Gnew.degree(M2)])
            M2 = list(M2)
            M = list(M)
            PossibleSet = Degs-2*Degs2
            Wh = np.where(PossibleSet>=0)
            Wh = Wh[0]
            
            OthNodes = np.setdiff1d(M,M2)
            for i in range(len(Wh)):
                CurrentNode = M2[Wh[i]]
                Gpos = COPY.deepcopy(Gnew)
                PEdges1 = np.array([val for (node, val) in Gold.edges(CurrentNode)])
                PEdges2 = np.array([val for (node, val) in Gnew.edges(CurrentNode)])
                EdgeDiff = np.setdiff1d(PEdges1,PEdges2)
                Intersect = np.intersect1d(OthNodes,EdgeDiff)
                if len(Intersect)>=len(PEdges2):
                    Removed = []
                    Added = []
                    for j in range(len(PEdges2)):
                        Removed.append((CurrentNode,PEdges2[j]))
                        Added.append((CurrentNode,Intersect[j]))
                        Gpos.remove_edge(CurrentNode,PEdges2[j])
                        Gpos.add_edge(CurrentNode,Intersect[j])
                    
                    if self.size_LCC(Gpos)<self.size_LCC(Gnew):
                        Gnew.add_edges_from(Added)
                        Gnew.remove_edges_from(Removed)
                        
                    
            
            
            
            return Gnew
    
    def rewire_to_fix_greedy(self):
        """
        Aold and Anew must have been added to the class...
        Aold - The original (n x n) adjacency matrix
           Anew - The "new" (n x n) adjacency matrix after edges have been removed.
        
        NOTE: Currently it is assumed that if greedy removal is unable to break up
        the LCC that the rewiring step could not be used to break up the net. 
        It is not clear if the above assumption is true, but we will assume it so
        for now.
        """
        if len(self.Aold)==1:
            raise ValueError("Must add Aold using the method add_old_adjacency_matrix")
        
        if len(self.Anew) ==1:
            raise ValueError("Must add Anew using the method add_new_adjacency_matrix (or by destroying the network using the greedy removal method")
        
        
        Aoldcopy = COPY.deepcopy(self.Aold)
        Anewcopy = COPY.deepcopy(self.Anew)
        
        Gold = nx.Graph(Aoldcopy)
        Gnew = nx.Graph(Anewcopy)
        M = self.size_LCC(Gnew)
        
        if M == len(Gold.nodes):
            self.Anew=Anewcopy
            return Anewcopy
        else:
           
            #Find the nodes that live in the LCC
            M = self.components_LCC(Gold)
            M2 = self.components_LCC(Gnew)
            
            #Find the node which could potentially be rewired
            Degs = np.array([val for (node, val) in Gold.degree(M2)])
            Degs2 = np.array([val for (node, val) in Gnew.degree(M2)])
            M2 = list(M2)
            M = list(M)
            PossibleSet = Degs-2*Degs2
            Wh = np.where(PossibleSet>=0)
            Wh = Wh[0]
            
            OthNodes = np.setdiff1d(M,M2)
            for i in range(len(Wh)):
                CurrentNode = M2[Wh[i]]
                Gpos = COPY.deepcopy(Gnew)
                PEdges1 = np.array([val for (node, val) in Gold.edges(CurrentNode)])
                PEdges2 = np.array([val for (node, val) in Gnew.edges(CurrentNode)])
                EdgeDiff = np.setdiff1d(PEdges1,PEdges2)
                Intersect = np.intersect1d(OthNodes,EdgeDiff)
                if len(Intersect)>=len(PEdges2):
                    Removed = []
                    Added = []
                    for j in range(len(PEdges2)):
                        Removed.append((CurrentNode,PEdges2[j]))
                        Added.append((CurrentNode,Intersect[j]))
                        Gpos.remove_edge(CurrentNode,PEdges2[j])
                        Gpos.add_edge(CurrentNode,Intersect[j])
                    
                    if self.size_LCC(Gpos)<self.size_LCC(Gnew):
                        Gnew.add_edges_from(Added)
                        Gnew.remove_edges_from(Removed)
                        Anew = nx.adjacency_matrix(Gnew)
                        self.Anew = Anew.todense()
                    
            
            
            
            return self.Anew
    
    def sparse_iterative_add_back(self,c):
        """A sparse (networkx graph) method for performing the final stage of the fragility 
           estimation algorithm. Adding back in any edges which do not increase
           the size of the largest connected components.
           """
        Gold = COPY.deepcopy(self.Graph)
        Gnew = COPY.deepcopy(self.Gnew)

        Gnew2 = COPY.deepcopy(Gnew)
        #SizeLCC = self.size_LCC(Gnew)
        Edges = list(Gold.edges())
        for edge in Edges:
            Gnew2.add_edge(edge[0], edge[1])
            if self.size_LCC(Gnew2)<c+1:
                Gnew.add_edge(edge[0],edge[1])
            
            Gnew2 = COPY.deepcopy(Gnew)
        self.Gnew = Gnew
        return Gnew
   
    def exact_erdos_renyi_graph(self,n,m):
        """n - (positive integer) The number of nodes in the ER network
           m - (positive integer) The exact number of undirected edges...
           
           returns- A - (n x n) Adjacency Matrix, containing an undirected ER
                        network with exactly m edges.
                        """
        A = np.zeros((n,n))
        Arr = np.arange(n)
        for i in range(m):
            np.random.shuffle(Arr)
            source = Arr[0]
            target = Arr[1]
            if A[source,target] ==1:
                while A[source,target] == 1:
                    np.random.shuffle(Arr)
                    source = Arr[0]
                    target = Arr[1]
                
            A[source,target] = 1
            A[target,source] = 1
        
        self.A = A
        self.Aold = A
        return self.A
