# Meta-Tree Labeling

The program implements the proposed Meta-Tree Indexing framework in [**Structure-Based Shortest-Path Distance Queries on Large Networks**](http://courses.cecs.anu.edu.au/courses/CSPROJECTS/17S1/Final_presentations/Presentations%20by%20student/Hao_Zhang_Final.pdf).

### Abstract
The shortest-path distance has been one of the most important concepts in graphs and has wide range of applications. To query distances in large networks, a subset of distances between vertex-pairs are usually precomputed and stored as an index, e.g. 2-hop index. However, finding a proper subset of vertex-pairs efficiently remains a challenge for large networks. Thus, several methods have been proposed recently by exploiting the existence of highly central vertices and core-fringe structure of networks. As distance queries are often used in social networks and other scenarios where structures like communities are inborn, it is worth exploring the possibility of improving indexing methods by making use of network structures.
In my paper, a novel distance indexing method is proposed. The method extracts a meta-tree from a network through partitioning the network by cutting articulation points and taking each induced subgraph as a meta-node. For a pair of vertices in a network, the meta-tree of the network provides a meta-path consisting of a sequence of articulation points on the shortest path between the pair. The sequence of points divides the shortest path between the vertex-pair into segments each of which can be considered as an independent sub-problem in corresponding subgraph therefore can be evaluated in parallel.
Extensive experiments on real-world datasets show the proposed method outperforms the state-of-the-art baseline method published on SIGMOD with on average 30% less index size and up to 50% reduced query time.
Moreover, the meta-tree structure can be easily extended to solve (multi) shortest path and many other related problems in graphs.

### Meta Trees
Meta tree structures extracted from real-world datasets

![Meta-Trees](https://github.com/zenithanu/MetaTree/blob/master/ScreenShots/MetaTrees.png)


### Classes:

	- meta_tree.h
	Meta-Tree class that is able to construct a tree structure from a graph

	- graph.h
	A graph class that is able to find articulation points of a graph and rank them by number of pair of vertices they disconnect

	- pruned_landmark_labeling.h
	Implementation of pruned landmark labeling that indexes a large graph based on the notion of 2-hop cover

### Programs:

	- query_main.cpp
	Index a graph using Meta-Tree Labeling and let user query distances using the produced index

	- benchmark_main.cpp
	Compare results, query time, index size of Meta-Tree Labeling to that of Pruned Landmark Labeling on random vertex-pairs and summarise the performance

	- output_metatree_main.cpp
	Store edge list and node weight of Meta-Tree of a graph to a file for visualisation
	
### Compile:
	make clean && make all
* Note: For MacOS, substitute <malloc.h> with <malloc/malloc.h> in pruned_landmark_labeling.h.

### Examples:

	- Build index and query distances on sample dataset
	./cpll Datasets/Email-Enron.txt

	- Run 100 tests on sample dataset for the top-500 cut
	./benchmark Datasets/Email-Enron.txt 100 500

	- Output Meta-Tree of sample dataset to a file
	./output_metatree Datasets/Email-Enron.txt Email_EdgeList.txt
