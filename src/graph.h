/* @author: Mojtaba Rezvani, Hao Zhang
*
 * A graph class with capability of finding articulation pointts
 * and rank by number of vertex-pairs disconnected
 */

#ifndef _UTILITY_H_
#define _UTILITY_H_
#endif
#define _LINUX_
#define _PRINT_NPR
#define _PRINT_ACS

#ifdef _LINUX_
	#include <sys/time.h>
#endif
#include <iostream>
#include <cstring>
#include <fstream>
#include <stack>
#include <vector>

using namespace std;


/////////////////////////////////////// Begin Struct
struct List 
{
	int index;
	struct List *next;
	List(): next(NULL)
	{}
};

struct Node
{
	public:
	struct List *first, *last;
	unsigned long int lbound;

	Node(): first(NULL), last(NULL), lbound(0)
	{
	}

	~Node()
	{
		if (first == NULL) return ;
		List *e = first;
		while(e != NULL) {
			List *f = e-> next;
			delete e;
			e = f;
		}
	}
};

struct Graph
{
	public:
	char* dir;
	int n, m;
	Node *nodes;
	long double ZETA;
	long double* rank;
	int *disc;
	int *low;
	int ap_count;
	int *edge_labels;
	int *subtree_size;

	Graph();
	~Graph();
	void read_file(char *input);
	void read_edge_list(int k, std::vector<std::pair<int, int> > es);
	void add_edge(int a, int b);

	// Articulation Points Ranking
	int articulation_ranking_holes(int k);
	int analyze_articulation_ranking_holes(int k);
	int articulation_ranking(long double* r);
	int print_isolated_vertices(bool* removed);
	long double print_component_size(int* size);

	void set_size(int _n, int _m);
	int get_disc(int i);
	vector<pair<int, int> > get_tree();
	vector<pair<int, int> > get_subgraph(std::vector<int> nodes);
};
/////////////////////////////////////// End Struct


/////////////////////////////////////// Begin Graph Construction Classes
Graph::Graph()
{
	n = 0;
	m = 0;
	nodes = NULL;
	rank = NULL;
	edge_labels = NULL;
	subtree_size = NULL;
	disc = NULL;
	low = NULL;
}

Graph::~Graph()
{
	if (nodes != NULL) delete[] nodes;
	if (rank != NULL) delete[] rank;
	if (edge_labels != NULL) delete[] edge_labels;
	if (subtree_size != NULL) delete[] subtree_size;
	if (disc != NULL) delete[] disc;
	if (low != NULL) delete[] low;
}

void Graph::set_size(int _n, int _m)
{
	n = _n;
	m = _m;
	if (nodes != NULL) delete[] nodes;
	nodes = new Node[n];
}

// Construct the graph (undirected) from edge list
void Graph::read_edge_list(int k, std::vector<std::pair<int, int> > es)
{
	n = k;
	if (nodes != NULL) delete[] nodes;
	nodes = new Node[n];

	for(unsigned int i=0; i<es.size(); i++) {
		if (es[i].first == es[i].second) continue;
		int a = es[i].first;
		int b = es[i].second;
		add_edge(a,b);
		add_edge(b,a);
	}

	ZETA = ((long double) n)*((long double) n)*((long double) n);

	return ;
}

void Graph::read_file(char *input)
{
	char* dir = new char[100];
	strcpy(dir, input);

	ifstream in;
	in.open(input);

	in >> n >> m;
	nodes = new Node[n];

	for(int i=0; i<m; i++)
	{
		int a, b;
		in >> a >> b;
		add_edge(a,b);
		add_edge(b,a);
	}

	ZETA = ((long double) n)*((long double) n)*((long double) n);
	//ZETA = ((long double) n);

	return ;
}

void Graph::add_edge(int a, int b)
{
	List *aa;	
	if (nodes[a].first == NULL) 
	{
		nodes[a].first = new List;
		aa = nodes[a].first;
		nodes[a].last = aa;
	}
	else
	{
		nodes[a].last->next = new List;
		aa = nodes[a].last->next;
		nodes[a].last = aa;
	}
	aa->index = b;
	return ;
}

// Get a subgraph induced on nodes in the form of edge list
std::vector<std::pair<int, int> > Graph::get_subgraph(std::vector<int> nodeIds)
{
	std::vector<std::pair<int, int> > es;

	for (unsigned int i = 0; i < nodeIds.size(); i++)
	{
		int nodeId = nodeIds[i];
		for(List *e=nodes[nodeId].first; e!=NULL; e=e->next) {
			int nbr = e->index;
			// avoid duplication
			bool f = false;
			for  (unsigned int j = i+1; j < nodeIds.size(); j++)
				if (nbr == nodeIds[j]) {
					f = true;
					break;
				}
			if (f) es.push_back(std::make_pair(nodeId, e->index));
		}
	}

	return es;
}
/////////////////////////////////////// End Graph Construction Classes



/////////////////////////////////////// Begin Articulation Points Ranking

int Graph::articulation_ranking_holes(int k)
{
	//Rank
	long double* r = new long double[n];
	memset(r, 0, n*sizeof(long double));

	ap_count = articulation_ranking(r);

	//Find top-k
	if (k > ap_count) k = ap_count;
    for(int i=0; i<k; i++)
    {
        long double max = 0;
        int maxi = 0;
        for(int j=0; j<n; j++)
            if (max < r[j])
            {
                maxi = j;
                max = r[j];
            }
        r[maxi] *= (long double) (-1);
        cout << maxi << endl;
    }
            
	delete[] r;
	return ap_count;
}

int Graph::analyze_articulation_ranking_holes(int k)
{
	//Rank
	long double* r = new long double[n];
	memset(r, 0, n*sizeof(long double));

    bool *removed = new bool[n];
    memset(removed, false, n*sizeof(bool));
    cout << print_isolated_vertices(removed) << endl;

	ap_count = articulation_ranking(r);

	//Find top-k
	if (k > ap_count) k = ap_count;
	for(int i=0; i<k; i++)
	{
		long double max = 0;
		int maxi = 0;
		for(int j=0; j<n; j++)
			if (max < r[j])
			{
				maxi = j;
				max = r[j];
			}
		r[maxi] *= (long double) (-1);
        removed[i] = true;
        cout << i << ": " << maxi << " " << max << " " << print_isolated_vertices(removed) << endl;
	}

	delete[] r;
    delete[] removed;
	return ap_count;
}

int Graph::articulation_ranking(long double* r)
{
	// Init traversal tree variables
	disc = new int[n];
	low = new int[n];
	edge_labels = new int[2*n];
	memset(edge_labels, -1, 2*n*sizeof(int));
	subtree_size = new int[n];
	memset(subtree_size, 0, n*sizeof(int));

	// Mark all the vertices as not visited
	bool *visited = new bool[n];
	int *low_v = new int[n];
	int *parent = new int[n];
	bool *ap = new bool[n]; // To store articulation points
	List* *whereami = new List*[n];
	List* *whereitwas = new List*[n];
	 
	long long int *num_des = new long long int[n];
	long double *num_c = new long double[n];
	int *root = new int[n];
	int *size = new int[n];
	memset(num_des, 0, n*sizeof(long long int));
	memset(num_c, 0, n*sizeof(long double));
	memset(r, 0, n*sizeof(long double));
	long double bigNum = print_component_size(size);
	//print_components();

    // All the data needed for constructing the new semi-traversal tree
    // Number of edges in the tree
	int *rdisc = new int[n];
	int *redge = new int[n];
    int num_edges_in_tree = 0;


	// Initialize parent and visited, and ap(articulation point) arrays
	for (int i = 0; i < n; i++)
	{
		root[i] = -1;
		parent[i] = -1;
		visited[i] = false;
		ap[i] = false;
		whereitwas[i] = nodes[i].first;
		whereami[i] = nodes[i].first;
	}

	// Call the recursive helper function to find articulation points
	// in DFS tree rooted with vertex 'i'
	for (int i = 0; i < n; i++)
		if (!visited[i])
		{
			stack<unsigned int> s;
			int u = i;
			// A static variable is used for simplicity, we can avoid use of static
			// variable by passing a pointer.
			int time = 0;
			 
			// Count of children in DFS Tree
			int children = 0;
			 
			// Mark the current node as visited
			s.push(u);
			
			while (!s.empty())
			{
				u = s.top();

				// Initialize discovery time and low value
				if (!visited[u])
				{
					disc[u] = low[u] = ++time;
					low_v[u] = u;
                    rdisc[disc[u]] = u;
					visited[u] = true;
				}

				// do the post recursion operation
				for (; whereitwas[u] != whereami[u]; whereitwas[u]=whereitwas[u]->next)
				{
					int v = whereitwas[u]->index;  // v is current adjacent of u

					// Check if the subtree rooted with v has a connection to
					// one of the ancestors of u
					if (visited[v] && (parent[v]==u))
					{
						//low[u]  = min(low[u], low[v]);
                        if (low[u] > low[v])
                        {
                            low_v[u] = low_v[v];
                            low[u] = low[v];
                        }

						num_des[u] += (num_des[v]+1);
						int temp = num_des[v];
						num_des[v] = -1;
					 
						// u is an articulation point in following cases
					 
						// (1) u is root of DFS tree and has two or more chilren.
						if (parent[u] == -1 && children > 1)
						{
							ap[u] = true;
							if (temp != -1)
							{
								r[u] += (n-temp-2)* (temp+1);
								num_c[u] += (temp+1);
							}
						}
					 
						// (2) If u is not root and low value of one of its child is more
						// than discovery value of u.
						if (parent[u] != -1 && low[v] >= disc[u])
						{
							ap[u] = true;
							if (temp != -1)
							{
								r[u] += (n-temp-2)* (temp+1);
								num_c[u] += (temp+1);
							}
						}
					}

				}
				if (whereitwas[u] == NULL) s.pop();

				// do the pre recursion operation
				// Go through all vertices adjacent to this (that haven't been visited)
				bool flag = false;
				for (;whereami[u] != NULL; whereami[u]=whereami[u]->next)
				{
					int v = whereami[u]->index;  // v is current adjacent of u
						 
					// If v is not visited yet, then make it a child of u
					// in DFS tree and recur for it
					if (!visited[v])
					{
						if (u==i) children++;
						if (v!=i) root[v] = i;
						parent[v] = u;
						s.push(v);
						flag = true;
						// not necessary whereami[u] = whereami[u]->next;
                        //
                        // This is a tree edge and put it in the edge-list array of the traversal tree
                        edge_labels[2*num_edges_in_tree] = u;
                        edge_labels[2*num_edges_in_tree+1] = v;
                        redge[v] = 2*num_edges_in_tree;
                        num_edges_in_tree++;
						break ;
					}
					// Update low value of u for parent function calls.
					else if (v != parent[u])
						low[u]  = min(low[u], disc[v]);
				}
				if (flag) {
					continue;
				}
			}
		}

	// n should be replaced by the number of nodes in that connected component; but it doesn't really matter
	// that can be achieved by root[i] <- u && num_des[root[i]]
	for(int i=0; i<n; i++)
    {
        subtree_size[i] = (int) num_c[i];
		if ((ap[i]) && (root[i]!=-1)) 
		{
			r[i] += (n-num_des[root[i]]+num_c[i]-1)*(num_des[root[i]]-num_c[i]);
			r[i] += bigNum - size[i]*(n-size[i]-1);
		}
		else if (ap[i])
		{
			r[i] += bigNum - size[i]*(n-size[i]-1);
		}
    }


	delete[] visited;
	delete[] low_v;
	delete[] parent;
	delete[] whereami;
	delete[] whereitwas;

	// ranking helpers
	delete[] num_des;
	delete[] num_c;
	delete[] root;
	delete[] size;

	ap_count = 0;
	for (int i=0; i<n; i++)
		if (ap[i])
			ap_count++;

	delete[] ap;

	delete[] rdisc;
	delete[] redge;

	return ap_count;
}

int Graph::print_isolated_vertices(bool *removed)
{
	bool *visited = new bool[n];
	memset(visited, false, n*sizeof(bool));
	unsigned int *q = new unsigned int[n];
	
	unsigned int qq, qqq, c;
    int size_1 = 0;

	for(int i=0; i<n; i++)
	{
		if ((removed[i]) || (visited[i])) continue;
		
		c = 1;
		qq = 0;
		qqq = 0;
		q[qqq++] = i;
		visited[i] = true;
		while(qqq > qq)
		{
			int u = q[qq++];
			for(List *e=nodes[u].first; e != NULL; e=e->next)
				if ((!removed[i]) && (!visited[e->index]))
				{
					q[qqq++] = e->index;
					visited[e->index] = true;
					c++;
				}
		}
        if (c==1) size_1++;
	}

	delete[] visited;
	delete[] q;

	return size_1;
}

long double Graph::print_component_size(int* size)
{
	bool *visited = new bool[n];
	memset(visited, false, n*sizeof(bool));
	unsigned int *q = new unsigned int[n];
	
	unsigned int qq, qqq, c;
	long double sum = 0;

	for(int i=0; i<n; i++)
	{
		if (visited[i]) continue;
		
		c = 1;
		qq = 0;
		qqq = 0;
		q[qqq++] = i;
		visited[i] = true;
		while(qqq > qq)
		{
			int u = q[qq++];
			for(List *e=nodes[u].first; e != NULL; e=e->next)
				if (!visited[e->index])
				{
					q[qqq++] = e->index;
					visited[e->index] = true;
					c++;
				}
		}
		sum += c * (n-c-1);
		for (int j=qqq-1; j>=0; j--)
			size[q[j]] = c;
	}

	delete[] visited;
	delete[] q;

	return sum;
}
/////////////////////////////////////// End Articulation Points Ranking



/////////////////////////////////////// Begin Traversal Tree @Hao
vector<pair<int, int> > Graph::get_tree()
{
	vector<std::pair<int, int> > tes;
	// Init articulation point scores
	rank = new long double[n];
	ap_count = articulation_ranking(rank);
    
    // Construct edge list of traversal tree
    int i = 0;
    while(edge_labels[i] != -1)
    {
    	tes.push_back(std::make_pair(edge_labels[i], edge_labels[i+1]));
        i+=2;
    }

    return tes;
}
/////////////////////////////////////// End Traversal Tree