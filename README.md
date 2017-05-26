# Meta-Tree Labeling

The program implements the proposed Meta-Tree Indexing framework in [**Structure-Based Shortest-Path Distance Queries on Large Networks**](http://courses.cecs.anu.edu.au/courses/CSPROJECTS/17S1/Final_presentations/Presentations%20by%20student/Hao_Zhang_Final.pdf).

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
	
	* For MacOS, substitute <malloc.h> with <malloc/malloc.h> in pruned_landmark_labeling.h.

### Examples:

	- Build index and query distances on sample dataset
	./cpll Datasets/Email-Enron.txt

	- Run 100 tests on sample dataset for the top-500 cut
	./benchmark Datasets/Email-Enron.txt 100 500

	- Output Meta-Tree of sample dataset to a file
	./output_metatree Datasets/Email-Enron.txt Email_EdgeList.txt
  ![Screenshot](https://github.com/zenithanu/MetaTree/blob/master/ScreenShots/Benchmark.png)

* Note: Detail parameters descriptions are on top of each file
