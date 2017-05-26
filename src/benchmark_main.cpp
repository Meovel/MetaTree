// @author: Hao Zhang, 2017
//
// Test program for Meta-Tree Lebling
// The program generates random distance queries on a graph 
// check correctness of the results 
// and compare performance with Pruned Landmark Labeling
//
// @params: filename     - path to the graph file (edge list)
//          n (optional) - number of pairs of vertices to be tested
//          k (optional) - top-k articulation points to be cut, default all 

#include "meta_tree.h"
#include <math.h>

int main(int argc, char* argv[]) {
    
  // Input Sanity Check 
  if (argc < 2) {
      printf("No dataset provided");
      exit(EXIT_FAILURE);
  }

  int k = 0;
  if (argc > 3) k = atoi(argv[3]);

  // disable stdout buffer for imediate result
  setbuf(stdout, NULL);

  // Construct index using different methods
  MetaTree metaTree;
  metaTree.ConstructMetaTree(argv[1], k);
  metaTree.IndexMetaTree();
  
  PrunedLandmarkLabeling<> pll;
  printf("- Benchmark Index ...\n");
  pll.ConstructIndex(argv[1], false, false);
  printf("- FINISH -\n\n");

  // get vertices vector for random access
  int nNodes = metaTree.GetNodes();
  std::vector<int> vVec = metaTree.GetNodesVector();

  
  // test on n random pairs of vertices
  int n = 10;
  if (argc > 2) n = atoi(argv[2]);
  int passed = 0;
  double mtTime = 0;
  double mtTimeP = 0;
  double pllTime = 0;
  // format settings
  int ndigit = ceil(log10(n+1));

  srand((unsigned)time(0));

  for (int i = 0; i < n; i++) {
  	int v = vVec[rand()%(nNodes-1)];
  	int w = vVec[rand()%(nNodes-1)];
    printf("Test %*d - D(%d, %d) = ", ndigit, i+1, v, w);

    mtTime -= GetCurrentTimeSec();
    int mtDistance = metaTree.QueryDistance(v, w);
    mtTime += GetCurrentTimeSec();
    mtTimeP += metaTree.GetLastQueryTime();

    pllTime -= GetCurrentTimeSec();
  	int pllDistance = pll.QueryDistance(v, w);
    pllTime += GetCurrentTimeSec();

    // output test result
  	printf("%d ", mtDistance);
  	if (mtDistance == pllDistance) {
  		printf("... PASS\n");
  		++ passed;
  	}
  	else 
  		printf("... FAIL (Correct %d)\n", pllDistance);
  }

  printf("\n%d / %d Tests Passed\n\n", passed, n);

  // Compare index size
  printf("Average Label Size for Pruned Landmark Labeling: %.1lf\n", pll.GetAvgLabelSize());
  printf("Average Label Size for Meta-Tree Labeling:       %.1lf (%.1lf%%)\n\n", 
  metaTree.GetAvgLabelSize(), metaTree.GetAvgLabelSize()/pll.GetAvgLabelSize()*100);

  // Compare preprocessing time
  printf("Preprocessing time for Pruned Landmark Labeling: %f s\n", pll.GetIndexTime());
  printf("Preprocessing time for Meta-Tree Labeling:       %f + %f + %f s\n", 
    metaTree.GetMetaTreeTime(), metaTree.GetSubgraphTime(), metaTree.GetIndexTime());
  printf("                                                 MetaTree + Subgraph + Index\n\n");

  // Compare query time
  printf("Average Query Time for Pruned Landmark Labeling: %f ms\n", pllTime/n*1E6);
  // Sequential and estimated parallel query time
  printf("Average Query Time for Meta-Tree Labeling      : %f ms | %f ms (%.1lf%%)\n",
   mtTime/n*1E6, mtTimeP/n*1E6, mtTimeP/pllTime*100);
  printf("                                                 Sequential  | Parallel\n");

  return 0;  
}