// @author: Hao Zhang, 2017
//
// The program produces Meta-Tree for a graph
// and store it as edge list in a file
// for the convenience of visulization
//
// @params: graph_file    - path to the graph file
//          edgelist_file - path to the output file that stores meta-tree (edge list)
//          nodeweight_file (optional) - path to the output file that stores meta-nodes and sizes

#include "meta_tree.h"


int main(int argc, char* argv[]) {
  // Input Sanity Check 
  if (argc < 3) {
      printf("Error: Expect two params - graph_file output_file");
      exit(EXIT_FAILURE);
  }
  // File sanity check
  std::ofstream file;
  file.open(argv[2], ios_base::out | ios_base::in); // will not create file
  if (file.is_open()) {
    cout << "Error - output file already exists";
    exit(EXIT_FAILURE);
  }

  
  // Construct Meta-Tree
  MetaTree t;
  t.ConstructMetaTree(argv[1]);
  EdgeList es = t.GetEdgeList();


  // Output Edge List
  file.clear();
  file.open(argv[2], ios_base::out); 

  for (unsigned int e = 0; e < es.size(); e++) {
    file << es[e].first << "\t" << es[e].second << endl;
  }

  file.close();
  printf("- EdgeList Output Finish -\n");

  // Output Node Wight List
  if (argc > 3) {
    file.open(argv[3], ios_base::out | ios_base::in);
    if (file.is_open()) {
      cout << "Error - output file already exists";
      exit(EXIT_FAILURE);
    }

    es = t.GetNodeSizes();

    file.clear();
    file.open(argv[3], ios_base::out); 

    file << "Id" << "\t" << "Size" << endl;
    for (unsigned int e = 0; e < es.size(); e++) {
      file << es[e].first << "\t" << es[e].second << endl;
    }

    file.close();
  }

  printf("- NodeWeight Output Finish -\n");
  
  return 0;  
}