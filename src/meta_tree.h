/* @author: Hao Zhang, 2017
 * All rights reserved.
 *
 * Implementation of proposed algorithm in 
 * <<Structure-Based Shortest-Path Distance Queries on Large Networksrk>>
 */

#include "graph.h"
#include "pruned_landmark_labeling.h"
#include <unordered_set>

typedef std::vector<std::pair<int, int> > EdgeList;
typedef std::map<int, std::vector<int> > MetaNodeMap;
typedef std::vector<PrunedLandmarkLabeling<> > PLLV;


class MetaTree {
public:
  // Construct and index meta-tree from a graph
  // Returns |true| when successful
  bool ConstructMetaTree(std::vector<std::pair<int, int> > &es, int k, const bool nodeIdsContinuous=true);
  bool ConstructMetaTree(std::istream &ifs, int k, const bool nodeIdsContinuous=true);
  bool ConstructMetaTree(const char *filename, int k, const bool nodeIdsContinuous=true);
  bool IndexMetaTree();

  // Calculate the shortest-path distance between vertex v and w
  int QueryDistance(int v, int w);

  // Calculate statistics
  double GetAvgLabelSize();

  // Access propertities
  int GetNodes() { return n; };
  int GetEdges() { return m; };
  int GetMaxNId() { return maxNId; };
  double GetMetaTreeTime() { return time_metatree_; };
  double GetSubgraphTime() { return time_subgraph_; };
  double GetIndexTime() { return time_indexing_; };
  double GetLastQueryTime() { return time_query_; };

  // Produce edge list of meta-tree
  // can be used for visualization
  EdgeList GetEdgeList();
  EdgeList GetNodeSizes();
  // output Meta-Tree as edge list
  bool OutputMetaTree(std::ostream &ofs);
  bool OutputMetaTree(const char *filename);

  std::vector<int> GetNodesVector();

  std::vector<std::pair<int, int> > get_subgraph(std::vector<int> nodes);

  void Free();

  MetaTree()
      : nodeIdBV(new bool[0]), metaNodeBV(new int[0]), n(0), m(0), maxNId(0), 
      time_metatree_(0), time_subgraph_(0), time_indexing_(0), time_query_(0) {}
  virtual ~MetaTree() {
    Free();
  }

private:
  // graph
  Graph G;
  std::set<int> nodes;
  // adjacency matrix
  std::map<int, std::vector<int> > adjm;
  // articulation points
  std::vector<int> aps;
  // bit vectors
  bool* nodeIdBV;
  int* metaNodeBV;
  // meta tree
  MetaNodeMap metaNodeMap;
  // labels
  PLLV metaNodeLabels;
  std::map<int, int> metaNodeLabelMap;

  // propertities
  int n, m;
  int maxNId;
  int k_ap;

  // time
  double time_metatree_, time_subgraph_, time_indexing_, time_query_;
};


double GetCurrentTimeSec() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}


// Read graph from edge list/file/input stream
bool MetaTree::ConstructMetaTree(const char *filename, int k=0, const bool nodeIdsContinuous) {
  std::ifstream ifs(filename);
  printf("- %s ", filename);
  return ifs && ConstructMetaTree(ifs, k, nodeIdsContinuous);
}

bool MetaTree::ConstructMetaTree(std::istream &ifs, int k=0, const bool nodeIdsContinuous) {
  EdgeList es;
  n = 0;
  m = 0;
  maxNId = 0;

  // produce edge list, calculate number of nodes and edges
  for (int v, w; ifs >> v >> w; ) {
    es.push_back(std::make_pair(v, w));
    adjm[v].push_back(w);
    adjm[w].push_back(v);
    nodes.insert(v);
    nodes.insert(w);
    ++m;
    if (std::max(v,w) > maxNId) maxNId = std::max(v,w);
  }

  if (ifs.bad()) return false;

  maxNId++;
  n = nodes.size();
  printf("contains %d nodes\n", n);

  // bit vector record existence of node ids for multiple usage 
  // e.g. handle non-continuous ides
  nodeIdBV = new bool[maxNId];
  memset(nodeIdBV, false, maxNId*sizeof(bool));
  for (set<int>::iterator i = nodes.begin(); i != nodes.end(); i++) 
     nodeIdBV[*i] = true;
  
  return ConstructMetaTree(es, k, nodeIdsContinuous);
}

bool MetaTree::ConstructMetaTree(std::vector<std::pair<int, int> > &es, int k=0, const bool nodeIdsContinuous) {
  G.read_edge_list(maxNId, es);
  setbuf(stdout, NULL);

  printf("- Analyse Articulation Points ... ");
  time_metatree_ = -GetCurrentTimeSec();
  // get edge list for traversal tree
  EdgeList tes = G.get_tree();
  // get score for articulation points
  int ap_count = G.ap_count;

  // get top-k articulation points
  if (k>ap_count || k==0) k = ap_count;
  k_ap = k;
  time_metatree_ += GetCurrentTimeSec();

  // // TODO: O(nlogn) sort
  for (int i=0; i<k; i++) {
    long double max = 0;
    int maxi = 0;
    for (int j = 0; j < maxNId; j++)
      if (nodeIdBV[j] && max<G.rank[j]) {
        maxi = j;
        max = G.rank[j];
      }
    G.rank[maxi] *= (long double) (-1);
    aps.push_back(maxi);
  }
  printf("%d Found\n", k);
  time_metatree_ -= GetCurrentTimeSec();

  
  // Construct a traversal tree
  Graph GTree;
  GTree.set_size(maxNId+1, maxNId);
  for (unsigned int i = 0; i < tes.size(); i++)
    GTree.add_edge(tes[i].first, tes[i].second);


  //////////////////////// Begin Meta-Tree
  printf("- Construct Meta-Tree ...\n");

  // init bit vector
  // -1 is used as a virtual articulation point of the root meta-node 
  metaNodeBV = new int[maxNId];
  memset(metaNodeBV, -1, maxNId*sizeof(int));

  for (int i = 0; i < k; i++) {
    int ap = aps[i];
    // printf("Processing Articulation point %d\n", ap);
 
    // BFS on subtrees to find elements for the new meta-node
    std::queue<int> candidates;
    // candidates.push(ap);

    int metaNodeId = metaNodeBV[ap];
    int disc = G.disc[ap];
    for(List *e=GTree.nodes[ap].first; e!=NULL; e=e->next) {
      int child = e->index;
      if (G.low[child] >= disc) {
        candidates.push(child);
      }
    }

    while (!candidates.empty()) {
      int candidate = candidates.front();
      candidates.pop();
      metaNodeBV[candidate] = ap;

      // iterate out-edges to visit children
      for(List *e2=GTree.nodes[candidate].first; e2!=NULL; e2=e2->next) {
        int nbr = e2->index;
        if (metaNodeBV[nbr]==metaNodeId) {
          // annotate visited nodes as new meta-node
          metaNodeBV[nbr] = ap;
          candidates.push(nbr);
        }
      }
    }
    metaNodeBV[ap] = metaNodeId;
  }
  time_metatree_ += GetCurrentTimeSec();

  // Get elements of each Meta-Node
  for (int i = 0; i < maxNId; i++) 
    if (nodeIdBV[i]) metaNodeMap[metaNodeBV[i]].push_back(i);
  // copy articulation points to their descendant meta-nodes
  for (MetaNodeMap::iterator it = metaNodeMap.begin(); it != metaNodeMap.end(); it++) 
    if (it->first != -1) (it->second).push_back(it->first);
  
  //////////////////////// End Meta-Tree

  return true;
}


bool MetaTree::IndexMetaTree() {
  printf("- Meta-Tree Index ... ");
  // Meta-Tree contains k+1 nodes and k edges
  int metaTreeSize = k_ap+1;
  metaNodeLabels.resize(metaTreeSize);
  // add virtual articulation point (meta-tree root)
  aps.push_back(-1);

  // Use pruned landmark labeling to index each Meta-Node
  time_subgraph_ -= GetCurrentTimeSec();
  unsigned int maxMetaNode = 0;
  for (int i = 0; i < metaTreeSize; i++) {
    int ap = aps[i];
    std::vector<int> nodesV = metaNodeMap[ap];
    metaNodeLabelMap[ap] = i;
    // pruned landmark labeling
    EdgeList es = G.get_subgraph(nodesV);
    // change 2nd param to true to enable deg-1 nodes pruning
    time_indexing_ -= GetCurrentTimeSec();
    metaNodeLabels[i].ConstructIndex(es, false, true);
    time_indexing_ += GetCurrentTimeSec();
    if (nodesV.size() > maxMetaNode) maxMetaNode = nodesV.size();
  }
  time_subgraph_ = time_subgraph_ + GetCurrentTimeSec() - time_indexing_;

  printf("Max Meta-Node Size: %d\n", maxMetaNode);
  
  return true;
}


int MetaTree::QueryDistance(int v, int w) {
  time_query_ = -GetCurrentTimeSec();

  // -1 distance means not connected
  // check if vertices exists
  if (v>maxNId || w>maxNId || !nodeIdBV[v] || !nodeIdBV[w]) {
    time_query_ += GetCurrentTimeSec();
    return -1;
  }

  int metaNodeV = metaNodeBV[v];
  int metaNodeW = metaNodeBV[w];

  // vertices in the same Meta-Node
  if (metaNodeV == metaNodeW) {
    int d = metaNodeLabels[metaNodeLabelMap[metaNodeV]].QueryDistance(v, w);
    time_query_ += GetCurrentTimeSec();
    return d;
  }

  
  //////////////////// Find articulation point sequence for the shortest path
  // method greatly reduces computation time for tree of small depth but big size
  // detail refer to thesis
  // set to -2 because -1 is used for virtual articulation point of root meta-node
  int commonParentId = -2;

  // Path from Meta-Node of v to Meta-Tree root
  std::vector<int> vToRoot;
  std::unordered_set<int> vToRootS;
  vToRoot.push_back(metaNodeV);
  int currentMetaNode = metaNodeV;
  vToRootS.insert(metaNodeV);
  while (currentMetaNode != -1) {
    currentMetaNode = metaNodeBV[currentMetaNode];
    vToRoot.push_back(currentMetaNode);
    vToRootS.insert(currentMetaNode);
  }

  // Path from Meta-Node of w to Meta-Tree root
  currentMetaNode = metaNodeW;
  if (std::find(vToRoot.begin(), vToRoot.end(), metaNodeW) != vToRoot.end()) 
    commonParentId = metaNodeW;
  else
    do {
      currentMetaNode = metaNodeBV[currentMetaNode];
      // hashset for O(1) element search - O(n) set intersection
      std::unordered_set<int>::const_iterator it = vToRootS.find(currentMetaNode);
      if (it != vToRootS.end()) {
        commonParentId = currentMetaNode;
        break;
      }
    } while (currentMetaNode != -1);

  time_query_ += GetCurrentTimeSec();

  // vertices not in the same commnected component
  if (commonParentId == -2) return -1;
  
  // printf("Nearest Common Parent: %d\n", commonParentId);

  //////////////////// Calculate shortest-path distance
  // the distance is sum of 3 parts
  // detail refer to thesis
  int distance = 0;
  double maxTime = 0.0;

  // 1. distance from v to the nearest common parent
  currentMetaNode = metaNodeV;
  int lastAPv = v;
  while (currentMetaNode != commonParentId) {
    int t = -GetCurrentTimeSec();
    int hopDistance = metaNodeLabels[metaNodeLabelMap[currentMetaNode]].QueryDistance(lastAPv, currentMetaNode);
    distance += hopDistance;
    lastAPv = currentMetaNode;
    currentMetaNode = metaNodeBV[currentMetaNode];
    t += GetCurrentTimeSec();
    if (t > maxTime) maxTime = t;
  }

  // 2. distance from v to the nearest common parent
  currentMetaNode = metaNodeW;
  int lastAPw = w;
  while (currentMetaNode != commonParentId) {
    int t = -GetCurrentTimeSec();
    int hopDistance = metaNodeLabels[metaNodeLabelMap[currentMetaNode]].QueryDistance(lastAPw, currentMetaNode);
    distance += hopDistance;
    lastAPw = currentMetaNode;
    currentMetaNode = metaNodeBV[currentMetaNode];
    t += GetCurrentTimeSec();
    if (t > maxTime) maxTime = t;
  }

  // 3. distance arcross the nearest common parent
  int t = -GetCurrentTimeSec();
  int connectD = metaNodeLabels[metaNodeLabelMap[commonParentId]].QueryDistance(lastAPv, lastAPw);
  t += GetCurrentTimeSec();
  if (t > maxTime) maxTime = t;
  time_query_ += maxTime;
  if (connectD >= 0) distance += connectD;
  else return -1;
  
  return distance;
}


double MetaTree::GetAvgLabelSize() {
  double totalLabelSize = 0.0;

  for (int i = 0; i <= k_ap; i++) {
    totalLabelSize += metaNodeLabels[i].GetLabelSize();
  }

  return totalLabelSize/n;
}


EdgeList MetaTree::GetEdgeList() {
  EdgeList es;

  for (int i = 0; i < k_ap; i++) {
    int ap = aps[i];
    int ancestor = metaNodeBV[ap];
    if (ancestor == -1) ancestor = maxNId;
    es.push_back(std::make_pair(ap, ancestor));
  }

  return es;
}


EdgeList MetaTree::GetNodeSizes() {
  EdgeList es;

  for (int i = 0; i < k_ap; i++) {
    int ap = aps[i];
    es.push_back(std::make_pair(ap, metaNodeMap[ap].size()));
  }
  es.push_back(std::make_pair(maxNId, metaNodeMap[-1].size()));

  return es;
}


std::vector<int> MetaTree::GetNodesVector() {
  std::vector<int> nodesV;

  for (int i = 0; i < maxNId; i++) 
    if (nodeIdBV[i]) nodesV.push_back(i);
  
  return nodesV;
}


void MetaTree::Free() {
  delete[] nodeIdBV;
  delete[] metaNodeBV;
}
