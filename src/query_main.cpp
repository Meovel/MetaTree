// @author: Hao Zhang, 2017
//
// The program builds index on dataset using Meta-Tree Labeling
// and let user query distance between vertex-pairs
//
// @params: filename     - path to the graph file (edge list)
//          k (optional) - top-k articulation points to be cut, default all

#include "meta_tree.h"


int main(int argc, char* argv[]) {
	int k;
    
    MetaTree t;
    
    if (argc > 2) k = atoi(argv[2]);
    else k = 0;

    t.ConstructMetaTree(argv[1], k);
    t.IndexMetaTree();
    printf("- FINISH -\n");

    printf("\n=== Ready for Query ===\n");

    for (int u, v; std::cin >> u >> v; ) 
        printf("%d\n", t.QueryDistance(u, v));

    return 0;
}