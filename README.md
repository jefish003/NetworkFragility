# NetworkFragility
Code for "How fragile is your network...."

There are two files containing a number of functions meant to make it easy and convenient to estimate the fragility of various networks, including Barabasi-Albert (BA), Erdos-Renyi (ER) and Watts-Strogatz (WS) networks. We will walk through the Examples.py file which is also provided. 

First we will start by importing all of the functions we need to estimate fragility and for saving purposes.

```
from NetworkFragilityClean import fragile_net
from RunNetworkFragilityClean import run_sparse_fn,run_sparse_worw
import scipy.io as IO
import networkx as nx
import numpy as np
from datetime import datetime 
```

Now we can get started. First get the setup out of the way:

```
Today = datetime.today().strftime('%d-%m-%Y')
#Setup the net
FN = fragile_net()
#Number of trials
NumTrials = 50 #This will be ignored if the network is the mall network...

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
NetworkType = Types[3]
```
Notice that in the above section we chose what we want delta to be. Since delta is directly related to c, we can also supply c very easily. Remember that delta is the ratio of the number of nodes in the Largest Connected Component (LCC) to the number of nodes in the network. Also the NorthMall graph is supplied as a mat file. We will finish out the setup now.

```
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
```

Now we will start destroying the network and estimating the fragility

```
if NetworkType == 'NorthMall':
    i = 0
    Fragilities = np.zeros(1) #Set the fragilities array
    Output = run_sparse_fn(G,delta,NetworkType+'Net_FullRemoval'+str(i))#
    np.savez(NetworkType+'Net_FullRemoval_' +str(i) + '_Edge' + Today+ '.npz',Output)
```

Yes it is that easy to destroy the network. Some things will print to the screen as a sanity check since for large networks this can take a really long time to run. We are not done however, we have yet to estimate the fragility, fortunately there is code provided for this as well.

```
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
```

Now the fragility has been estimated using rewiring and the minimum degree attack method. Note that the third step was implemented in the FN.sparse_iterative_add_back method, which completes the algorithm. However as noted in the paper we are not done for estimating fragility, since this attack method only works well for certain networks and certain deltas. The other attack method used in the paper is the edge betweenness attack method and we implement that method with rewiring and the iterative_add_back algorithm. 

```
  #Reset FN to eliminate graphs
    FN = fragile_net()
    
    EBDict = Output['EdgeBetweenness']
    #Grab the final graph after the edges have been removed
    G2 = nx.Graph(EBDict['Graphs'][-1])
    
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
```

Now we have completed the process for both the minimum degree and edge betweenness attack methods. To get the final estimate of fragility we need to take the maximum between the two (i.e the minimum number of edges which needed to be removed). 

```
#Now we can get the true estimated fragility
    Fragility = np.max(np.array([MDEstimatedFrag,EBEstimatedFrag]))
    Fragilities[i] = Fragility
    print("This is the estimated fragility: ", Fragility)
```

Note that these steps have been computed with rewiring. However if we wish we can also figure out which set of edges would be removed without rewiring as is done below.

```
#NOTE THAT THIS IS HOW FRAGILITY WAS ESTIMATED, HOWEVER ONE COULD LOOK AT
    #THE SITUATION WITHOUT REWIRING, THIS IS DONE BELOW
    #If desired this can be compared to the results without rewiring as part
    #of the process...
    Output2 = run_sparse_worw(G,delta,NetworkType+'Net_WORW_FullRemoval'+str(i))
    np.savez(NetworkType+'Net_WORW_FullRemoval_' +str(i) + '_Edge' + Today+ '.npz',Output2)
```
That is it! The rest of the Examples.py file simply does all of these steps, but for the other types of networks (BA,ER and WS). 
