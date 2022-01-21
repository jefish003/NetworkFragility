# NetworkFragility
Code for "How fragile is your network...."

There are two files containing a number of functions meant to make it easy and convenient to estimate the fragility of various networks, including Barabasi-Albert (BA), Erdos-Renyi (ER) and Watts-Strogatz (WS) networks. 

First we will start by importing all of the functions we need to estimate fragility and for saving purposes

```from NetworkFragilityClean import fragile_net
from RunNetworkFragilityClean import run_sparse_fn,run_sparse_worw
import scipy.io as IO
import networkx as nx
import numpy as np
from datetime import datetime '''

Now we can get started.
