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

