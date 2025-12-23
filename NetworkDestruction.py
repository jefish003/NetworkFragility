# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 13:26:54 2022

@author: fishja
"""

import networkx as nx
import numpy as np
from copy import deepcopy as deepcopy
from networkx.algorithms import community
import pymetis

class NetworkDestruction:
    
    def __init__(self,A=None,Graph=None,Undirected=True,LCC_Size=None):
        """A -(optional, but 1 of A or Graph must be present) The adjacency matrix, may be in sparse form
       Graph - A networkx graph, nodes are assumed labeled from 0 to n-1 (behavior will change if this is not the case!)
  Undirected - A bool indicating whether or not the graph is undirected
    LCC_Size - A positive integer between 1 and n-1, the size of the largest connected component after attack"""
        self.A = deepcopy(A)
        self.Graph = deepcopy(Graph)
        self.OldGraph = None
        self.Undirected = Undirected
        self.LCC_Size = LCC_Size
        self.n = None
        self.relabel = True
        #Assume a bipartition for METIS
        self.numparts = 2
    
    def set_numparts(self,NumPartitions):
        self.numparts = NumPartitions
        
    def set_relabel(self,TruthVal):
        self.relabel = TruthVal
        
    def set_A(self,A):
        #Set the adjacency matrix, it may be sparse
        self.A = deepcopy(A)
    
    def set_Graph(self,Graph):
        #Set networkx graph
        self.Graph = deepcopy(Graph)
    
    def set_LCC_Size(self,LCC_Size):
        self.LCC_Size = LCC_Size
        
    def clear_Graph(self):
        self.OldGraph = self.Graph
        self.Graph = None
        
    def return_LCC_Size(self):
        return self.LCC_Size
    
    def shuffle_Graph(self):
        if self.A is not None:
            index = np.arange(self.A.shape[0])
            np.random.shuffle(index)
            A = self.A[index,:]
            A = self.A[:,index]
            return A
        else:
            A = nx.adjacency_matrix(deepcopy(self.Graph))
            index = np.arange(A.shape[0])
            np.random.shuffle(index)
            A = A[index,:]
            A = A[:,index]
            if self.Undirected:
                Graph = nx.Graph(A)
            else:
                Graph = nx.DiGraph(A)
            return deepcopy(Graph)
    def sort_Graph(self):
        if self.Graph is not None:
            G = deepcopy(self.Graph)
            Deg =np.array(list(G.degree()))
            #Deg[:,1] = -Deg[:,1]
            Index = Deg[:,1].argsort()
            A = deepcopy(nx.adjacency_matrix(G))
            A = A[:,Index]
            A = A[Index,:]
            Gnew = deepcopy(nx.Graph(A))
            return Gnew
            self.Graph = deepcopy(Gnew)
        else:
            raise ValueError("Sorry only networkx graphs are supported")
    
    def Find_Max_Degree_Node(self):
        Graph = deepcopy(nx.Graph(self.Graph))
        Deg = np.array(list(Graph.degree()))
        Deg[:,1] = -Deg[:,1]
        Index = Deg[:,1].argsort() 
        return Index[0]
    def Set_Checker(self,n):
        Set1 = []
        Set2 = []
        Start = self.Find_Max_Degree_Node()
        Graph = nx.Graph(self.Graph)
        Set1.append(Start)
        Neighbors = np.array(list(Graph[Start]))
        Tester = False
        Diff = np.setdiff1d(np.arange(n),Start)
        for i in range(n-1):
            if len(np.intersect1d(Neighbors,Diff[i])) == 0:
                Set2.append(Diff[i])
            
            if len(Set2) == self.LCC_Size:
                Tester = True
                break
        
        if not Tester:
            Arange = np.arange(n)
            Arange = np.setdiff1d(Arange,np.array(Set1))
            Arange = np.setdiff1d(Arange,np.array(Set2))
            TotLength = len(Arange)
            D = np.zeros((TotLength,2))
            for i in range(TotLength):
                Node = Arange[i]
                D[i,0] = len(np.intersect1d(np.array(Set1),np.array(list(Graph[Node]))))
                D[i,1] = len(np.intersect1d(np.array(Set2),np.array(list(Graph[Node]))))
            
            D = D[:,1]-D[:,0]
            AS = D.argsort()
            Arange = Arange[AS]
            for i in range(TotLength):
                Node = Arange[i]
                Set2.append(Node)
                if len(Set2) == self.LCC_Size:
                    break
            Set2 = np.array(Set2)
            Set1 = np.setdiff1d(np.arange(n),Set2)
            self.Set1 = Set1
            self.Set2 = Set2
        else:
            Set2 = np.array(Set2)
            Set1 = np.setdiff1d(np.arange(n),Set2)
            self.Set1 = Set1
            self.Set2 = Set2
            
        
    def Check_Sets(self,Set1,Set2,n):
        #If the graph is directed, it will be converted to 
        #undirected for this step (helps in minimization stage)
        Graph = deepcopy(nx.Graph(self.Graph))
        for i in range(n):
            #print(i)
            Edges = np.array(list(Graph[i]))
            if i in Set1:
                IntersectIn = np.intersect1d(Edges,Set1)
                IntersectOut = np.intersect1d(Edges,Set2)
                DiffInt = len(IntersectOut)-len(IntersectIn)
                #Lint = len(IntersectOut)
                ChangeInLinks = np.zeros(len(Set2))
                for j in range(len(Set2)):
                    Edges2 = np.array(list(Graph[Set2[j]]))
                    IntersectIn2 = np.intersect1d(Edges2,Set2)
                    IntersectOut2 = np.intersect1d(Edges2,Set1)
                    DiffInt2 = len(IntersectOut2)-len(IntersectIn2)
                    if Set2[j] in IntersectOut:
                        ChangeInLinks[j] = DiffInt+DiffInt2-2
                    else:
                        ChangeInLinks[j] = DiffInt+DiffInt2
                #print(ChangeInLinks,ChangeInLinks.argmax())
                if np.any(ChangeInLinks>0):
                    Index = ChangeInLinks.argmax()
                    Index = Set2[Index]
                    #print(i,Index)
                    Set1 = np.setdiff1d(Set1,i)
                    Set2 = np.setdiff1d(Set2,Index)
                    Set1 = np.union1d(Set1,Index)
                    Set2 = np.union1d(Set2,i)               
                    
            else:
                IntersectIn = np.intersect1d(Edges,Set2)
                IntersectOut = np.intersect1d(Edges,Set1)
                DiffInt = len(IntersectOut)-len(IntersectIn)
                #Lint = len(IntersectOut)
                ChangeInLinks = np.zeros(len(Set1))
                for j in range(len(Set1)):
                    Edges2 = np.array(list(Graph[Set1[j]]))
                    IntersectIn2 = np.intersect1d(Edges2,Set1)
                    IntersectOut2 = np.intersect1d(Edges2,Set2)
                    DiffInt2 = len(IntersectOut2)-len(IntersectIn2)
                    if Set1[j] in IntersectOut:
                        ChangeInLinks[j] = DiffInt+DiffInt2-2
                    else:
                        ChangeInLinks[j] = DiffInt+DiffInt2
                #print(ChangeInLinks,ChangeInLinks.argmax())
                if np.any(ChangeInLinks>0):
                    Index = ChangeInLinks.argmax()
                    Index = Set1[Index]
                    #print(i,Index)
                    Set2 = np.setdiff1d(Set2,i)
                    Set1 = np.setdiff1d(Set1,Index)
                    Set2 = np.union1d(Set2,Index)
                    Set1 = np.union1d(Set1,i)
            
        self.Set1 = deepcopy(Set1)
        self.Set2 = deepcopy(Set2)
                
    def Check_Sets_Old(self,Set1,Set2,n):
        Graph = deepcopy(self.Graph)
        for i in range(n):
            Edges = np.array(list(Graph[i]))
            if i in Set1:
                Intersect = np.intersect1d(Edges,Set2)
                Lint = len(Intersect)
                #Set1IntersecNum[i] = len(Intersect)
                Diff = len(list(Graph[i]))-Lint
                if Lint > Diff:
                    #Try to swap with a node out of set, only accept swap
                    #if it reduces the total number of cross set edges
                    Set2Diff = np.setdiff1d(Set2,Intersect)
                    #TotDiffs = np.zeros(len(Set2Diff))
                    for j in range(len(Set2Diff)):
                        #Test if this node can be swapped
                        NodeNum = Set2Diff[j]
                        Edges2 = np.array(list(Graph[NodeNum]))
                        Intersect2 = np.intersect1d(Edges2,Set1)
                        Lint2 = len(Intersect2)
                        Diff2 = len(list(Graph[NodeNum]))-Lint2
                        #TotDiffs[j] = Lint2-Diff2
                        if Lint2>=Diff2:
                                Set1 = np.setdiff1d(Set1,i)
                                Set2 = np.setdiff1d(Set2,NodeNum)
                                Set1 = np.union1d(Set1,NodeNum)
                                Set2 = np.union1d(Set2,i)
                                break
                    #if np.any(TotDiffs>0):
                    #    BestNode = TotDiffs.argmax()
                    #    BestNode = Set2Diff[BestNode]
                    #    Set1 = np.setdiff1d(Set1,i)
                    #    Set2 = np.setdiff1d(Set2,BestNode)
                    #    Set1 = np.union1d(Set1,BestNode)
                    #    Set2 = np.union1d(Set2,i)
                        #For now the first time we find an acceptable swap for both
                        #we stop looking and do the swap. This may change in the future..
                        
                        
                
            else:
                Intersect = np.intersect1d(Edges,Set1)
                Lint = len(Intersect)
                #Set1IntersecNum[i] = len(Intersect)
                Diff = len(list(Graph[i]))-Lint
                if Lint > Diff:
                    #Try to swap with a node out of set, only accept swap
                    #if it reduces the total number of cross set edges
                    Set1Diff = np.setdiff1d(Set1,Intersect)
                    #TotDiffs = np.zeros(len(Set1Diff))
                    for j in range(len(Set1Diff)):
                        #Test if this node can be swapped
                        NodeNum = Set1Diff[j]
                        Edges2 = np.array(list(Graph[NodeNum]))
                        Intersect2 = np.intersect1d(Edges2,Set2)
                        Lint2 = len(Intersect2)
                        Diff2 = len(list(Graph[NodeNum]))-Lint2
                        if Lint2>=Diff2:
                            Set2 = np.setdiff1d(Set2,i)
                            Set1 = np.setdiff1d(Set1,NodeNum)
                            Set2 = np.union1d(Set2,NodeNum)
                            Set1 = np.union1d(Set1,i)
                            break                          
                        #TotDiffs[j] = Lint2-Diff2
                    #if np.any(TotDiffs>0):
                    #    BestNode = TotDiffs.argmax()
                    #    BestNode = Set1Diff[BestNode]
                    #    Set2 = np.setdiff1d(Set2,i)
                    #    Set1 = np.setdiff1d(Set1,BestNode)
                    #    Set2 = np.union1d(Set2,BestNode)
                    #    Set1 = np.union1d(Set1,i)
                        #For now the first time we find an acceptable swap for both
                        #we stop looking and do the swap. This may change in the future..
                     
        self.Set1 = deepcopy(Set1)
        self.Set2 = deepcopy(Set2)
    def NumEdgesToDestruction(self):
        if self.Graph is None:
            if self.A is None:
                raise ValueError("One of Graph or A cannot be None, please set one or the other")
            else:
                if self.Undirected:
                    Graph = nx.Graph(self.A)
                else:
                    Graph = nx.DiGraph(self.A)
                
        else:
            Graph = self.Graph
        
        if self.LCC_Size is None:
            raise ValueError("LCC_Size must be a positive integer between 1 and n-1")
       
        self.n = Graph.number_of_nodes()
        n = self.n
        if self.LCC_Size < 2 or self.LCC_Size>n-1:
            raise ValueError("LCC_Size must be a positive integer between 2 and n-1")
        
        LCC_Size = self.LCC_Size
        #self.sort_Graph()
        Graph = deepcopy(Graph)
        Deg = np.array(list(Graph.degree()))
        #Convert in case...
        Deg.astype(np.int64)
        Deg[:,1] = -Deg[:,1]
        Index = Deg[:,1].argsort()
        #Arange = np.append(Index[0:-1:2],np.setdiff1d(Index,Index[0:-1:2]))
        if self.LCC_Size >=n/2:
            Set1_NumNodes = n-LCC_Size
            Arange = np.array(list(nx.dfs_postorder_nodes(nx.Graph(deepcopy(Graph)),Index[0])))
            #Initialize the sets based on a tree from highest
            #degree node. If the graph is directed it is converted
            #to undirected for the bfs tree search
            #if multiple nodes have highest degree then one of
            #them is chosen "at random"
            #Arange = np.array(list(nx.bfs_tree(nx.Graph(Graph),Index[0])))
            #Arange =np.append(np.arange(0,n,2),np.arange(1,n,2))#np.arange(n)
            Set1 = Arange[0:Set1_NumNodes]
            Set2 = Arange[Set1_NumNodes:n]
            #del Arange
            #Need to check how many edges are in set vs out of set
            #Will need one more pass to determine the number of nodes...
            #Set1 = Index[0:Set1_NumNodes]
            #Set2 = Index[Set1_NumNodes:n]
            #self.Set_Checker(n)
            #Set1 = deepcopy(self.Set1)
            #Set2 = deepcopy(self.Set2)
            self.Check_Sets(Set1,Set2,n)
            #self.Set_Checker(n)
            #self.Check_Sets(deepcopy(self.Set1),deepcopy(self.Set2),n)
            #self.Check_Sets(deepcopy(self.Set1),deepcopy(self.Set2),n)
            Set1 = deepcopy(self.Set1)
            Set2 = deepcopy(self.Set2)
            #self.Check_Sets(Set1,Set2,n)
            #Set1 = deepcopy(self.Set1)
            #Set2 = deepcopy(self.Set2)
            #Now count all the cross edges and break them!
            Count = 0
            for i in range(n):
                Edges = np.array(list(Graph[i]))
                if i in Set1:
                    Intersect = np.intersect1d(Edges,Set2)
                    Count = Count + len(Intersect)
                    Tuples = [(i,Intersect[index]) for index in range(len(Intersect))]
                    Graph.remove_edges_from(Tuples)
                    #NewGraph
                    
                else:
                    Intersect = np.intersect1d(Edges,Set1)
                    Count = Count+len(Intersect)
                    Tuples = [(i,Intersect[index]) for index in range(len(Intersect))]
                    Graph.remove_edges_from(Tuples)
                
            return Count,Graph,Set1,Set2                 
                        
        else:
            raise ValueError("Sorry this functionality is not yet available.")
    
    def Metis_NumEdges_Multi(self,delta):
        """Find the number of edges which need to be removed so that the fraction of nodes in the LCC is at least delta. 
        NOTE: If delta is set to larger than 1-1/n it will be assumed that you want delta = 1-1/n. Similarly if delta < 2/n,
        it will be assumed that you want delta = 2/n (obviously you need to remove all edges otherwise...)
        Also Metis assumes undirected graphs...
        """
        Graph = deepcopy(nx.Graph(self.Graph))
        n = len(Graph)
        Val1 = 1/n
        Val2 = 1-Val1
        #Ensure that the value of delta makes sense, automatically convert if
        #delta does not make sense...
        if delta>Val2 or delta<2*Val1:
            if delta>Val2:
                print("Warning, delta is greater than 1-1/n, assuming you want to set delta to 1-1/n!")
                delta = Val2
            else:
                print("Warning, delta is less than 2/n, assuming you want to set delta to 2/n!")
                delta = 2*Val1
        
        #The maximum component size in number of nodes 
        K = int(np.ceil(delta*n))
        if K>= n/2:
            #The number of components which will be needed
            self.numparts = 2
        else:
            self.numparts = int(np.ceil(1/delta))
            
        
        #If the nodes are not labeled properly there will be issues. 
        #so the nodes are relabeled from 0 to n-1. If self.relabel is False
        #we can assume the nodes are properly labeled and skip this step.
        if self.relabel:
            NodeLabels = np.array(list(Graph.nodes))
            Map = dict(zip(NodeLabels,np.arange(n)))
            G = nx.relabel_nodes(Graph,Map)    
        A_list = []
        for k in range(len(G)):
            A_list.append(np.array(list(G[k])))
        
        #Perform the l-section (i.e for bisection l = 2...) using METIS
        #Note METIS is fast and generally does quite a good job relative
        #to other bisection methods like Kernighan-Lin.
        NumEdges,PartLabels = pymetis.part_graph(self.numparts,A_list)
        PartLabels = np.array(PartLabels)
        Sizes = np.zeros(self.numparts)
        
        for k in range(self.numparts):
            Sizes[k] = len(np.where(PartLabels==k)[0])
            
        
        #If FindEdges is False, then we can simply use the NumEdges calculated
        #by METIS, otherwise we will need to recalculate
        FindEdges = False
        #Check to see that all components are small enough. If not handle this
        #problem by greedily rewiring (i.e moving nodes to other sets)
        
        
        
        #####################
        # It seems you have this backward here....Should be moving from smaller
        # component to the bigger component. This will need to be fixed !!!!!!
        # please fix the next time you get a chance
        #####################
        if np.max(Sizes)>K:
            FindEdges = True
            Check = np.max(Sizes)
            MaxCompArg = Sizes.argmax()
            MaxComp = Sizes[MaxCompArg]
            while Check>K:
                Wh = np.where(PartLabels==MaxComp)[0]
                InternalEdges = np.zeros(len(Wh))
                for ijk in range(len(Wh)):
                    Edges = np.array(list(G[Wh[ijk]]))
                    NumIntersect = len(np.intersect1d(Wh,Edges))
                    InternalEdges[ijk] = NumIntersect
                
                #Find the node with the smallest internal edges. The idea being
                #that this will increase the number of edges which need to be removed
                #by the least amount.
                print(MaxComp,Wh)
                MinNode = Wh[InternalEdges.argmin()]
                #Now Search over the other sets
                OtherComps = np.setdiff1d(np.arange(self.numparts),PartLabels[MinNode])
                #Need to know who are the neighbors of the node which will be moved
                Edges = np.array(list(G[MinNode]))
                #Need to figure out how many intersections the node has with the
                #other sets.
                NumInter = np.zeros(self.numparts-1)
                #Keep track of all of the components
                #ComponentsDict ={}
                #ComponentsDict[MaxComp] = np.where(PartLabels == MaxComp)[0]
                for k in range(self.numparts-1):
                    Wh = np.where(PartLabels==OtherComps[k])[0]
                    NumInter[k] = len(np.intersect1d(Wh,Edges))
                
                #Want to move to the component with the maximum intersections...
                OthComp = OtherComps[NumInter.argmax()]
                
                #Now update the Partition label
                PartLabels[MinNode] = OthComp
                
                #Update the sizes of each component
                Sizes[MaxComp] = Sizes[MaxComp]-1   
                Sizes[OthComp] = Sizes[OthComp]+1
                
                Check = np.max(Sizes)
            
            #Now if neccessary, recalculate the number of edges removed
            print("Find Edges: ", FindEdges)
            if FindEdges:
                #Need to recount the number of edges
                NumEdges = 0
                #Find the components
                ComponentsDict = {}
                #Find the relevant subgraphs
                SubGraphsDict = {}
                CurrentNodes = np.arange(n)
                #Now loop to count the total number of edges removed
                for k in range(self.numparts-1):
                    Comp = np.where(PartLabels==k)[0]
                    ComponentsDict[k] = Comp
                    SubG = deepcopy(G.subgraph(Comp))
                    SubGraphsDict[k] = SubG
                    
                    
                    if k == 0:
                        CurrentNodes = np.setdiff(CurrentNodes,Comp)
                        SubG2 = deepcopy(G.subgraph(CurrentNodes))
                        Num1 = len(G.edges())
                        Num2 = len(SubG.edges())+len(SubG2.edges())
                        NumEdges = NumEdges + (Num1-Num2)
                    
                    else:
                        SubGMain = deepcopy(G.subgraph(CurrentNodes))
                        CurrentNodes = np.setdiff(CurrentNodes,Comp)
                        SubG2 = deepcopy(G.subgraph(CurrentNodes))
                        Num1 = len(SubGMain.edges())
                        Num2 = len(SubG.edges())+len(SubG2.edges())
                        NumEdges = NumEdges+(Num1-Num2)
                        
                        
                
        return NumEdges,PartLabels
    
    def Metis_NumEdges(self):
        Graph = deepcopy(nx.Graph(self.Graph))
        n = len(Graph)
        if self.relabel:
            NodeLabels = np.array(list(Graph.nodes))
            Map = dict(zip(NodeLabels,np.arange(n)))
            G = nx.relabel_nodes(Graph,Map)    
        A_list = []
        for k in range(len(G)):
            A_list.append(np.array(list(G[k])))
        NumEdges,PartLabels = pymetis.part_graph(self.numparts,A_list)
        return NumEdges,PartLabels
                
    def Kernighan_Lin_NumEdges(self):
        Graph = deepcopy(nx.Graph(self.Graph))
        n = len(Graph)
        if self.relabel:
            NodeLabels = np.array(list(Graph.nodes))
            Map = dict(zip(NodeLabels,np.arange(n)))
            Graph = nx.relabel_nodes(Graph,Map)
        
        Deg = np.array(list(Graph.degree()))
        #convert in case, currently only int degree supported
        #Deg = Deg.astype(np.int64)
        Deg[:,1] = -Deg[:,1]
        Index = Deg[:,1].argsort()        
        #Arange = np.array(list(nx.bfs_tree(Graph,Index[0])))
        Arange = np.array(list(nx.dfs_postorder_nodes(nx.Graph(deepcopy(Graph)),Deg[Index[0],0])))
        if len(Arange) >= int(np.ceil(n/2)):
        #Arange =np.append(np.arange(0,n,2),np.arange(1,n,2))#np.arange(n)
            Set1 = Arange[0:int(np.ceil(n/2))].astype('int32')
            Set2 = np.setdiff1d(np.array(list(Graph.nodes())),Set1).astype('int32')
            print(len(Set1),len(Set2))
            #Set2 = Arange[int(np.ceil(n/2)):n] 
        else:
            Diff = np.setdiff1d(np.array(list(Graph.nodes())),Arange)
            C = int(np.ceil(n/2))
            D = C-len(Arange)
            Set1 = np.union1d(Arange,Diff[0:D])
            Set2 = np.setdiff1d(np.arange(n),Set1)
                
        
        print(Set1)
        
        Part = community.kernighan_lin_bisection(Graph,[np.array(Set1),np.array(Set2)])
        Set1 = np.array(list(Part[0]))
        Set2 = np.array(list(Part[1]))
         #Now count all the cross edges and break them!
        Count = 0
        for i in range(n):
            Edges = np.array(list(Graph[Deg[i,0]]))
            if Deg[i,0] in Set1:
                Intersect = np.intersect1d(Edges,Set2)
                Count = Count + len(Intersect)
                Tuples = [(Deg[i,0],Intersect[index]) for index in range(len(Intersect))]
                Graph.remove_edges_from(Tuples)
                #NewGraph
                
            else:
                Intersect = np.intersect1d(Edges,Set1)
                Count = Count+len(Intersect)
                Tuples = [(Deg[i,0],Intersect[index]) for index in range(len(Intersect))]
                Graph.remove_edges_from(Tuples)
            
        return Count,Graph,Set1,Set2 

    def Girvan_Newman_NumEdges(self):
        Graph = deepcopy(nx.Graph(self.Graph))
        comp = community.girvan_newman(Graph)
        T = tuple(sorted(c) for c in next(comp))
        Set1 = np.array(T[0])
        Set2 = np.array(T[1])
        n = len(Graph)
        #Deg = np.array(list(Graph.degree()))                      
        Count = 0
        for i in range(n):
            Edges = np.array(list(Graph[i]))
            if i in Set1:
                Intersect = np.intersect1d(Edges,Set2)
                Count = Count + len(Intersect)
                Tuples = [(i,Intersect[index]) for index in range(len(Intersect))]
                Graph.remove_edges_from(Tuples)
                #NewGraph
                
            else:
                Intersect = np.intersect1d(Edges,Set1)
                Count = Count+len(Intersect)
                Tuples = [(i,Intersect[index]) for index in range(len(Intersect))]
                Graph.remove_edges_from(Tuples)
            
        return Count,Graph,Set1,Set2             
            
            
        
        