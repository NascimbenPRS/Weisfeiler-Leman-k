// File WL_Graph.hh
#ifndef WL_GRAPH_HH
#define WL_GRAPH_HH

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <set>
#include <tuple>
#include <stdlib.h>
#include <time.h>
#include <math.h>	/* ceil */

using namespace std;


class WL_Graph
{
public:
  WL_Graph(){}; // empty constructor
  WL_Graph(string file_name); // constructor from DIMACS formatted file
  WL_Graph(const vector<vector<bool>> adjMat1); // constructor from adjacency matrix (uncolored graph)
  WL_Graph(const vector<unsigned> nodeColor1, const vector<list<pair<unsigned,unsigned>>> adjList1);
  // constructor from vector of node colors and colored adjacency list (colored graph)
  // -- the resulting graph's adjacency matrix is empty, for memory usage minimisation
  WL_Graph(const WL_Graph& g1, const vector<unsigned> nodes); // return subgraph of g1 induced by 'nodes' (uncoloured)

  unsigned Order() const { return n; }
  unsigned Edges() const { return m; }
  bool IsColored() const { return isColored;}
  bool AdjMat(unsigned u, unsigned v) const { return adjMat[u][v];}
  list<pair<unsigned,unsigned>> AdjList(unsigned u) const { return adjList[u];}
  vector<unsigned> SortedDegrees() const { return sortedDegrees;}
  unsigned NodeColor(unsigned u) const { return nodeColor[u];}
  unsigned EdgeColor(unsigned u, unsigned v) const { return edgeColor[u][v];}

 private:
  unsigned n, m; // # nodes, edges
  bool isColored; // true iff graph has an initial coloring (on either nodes or edges)
  vector<vector<bool>> adjMat; // adjacency matrix of the graph
  vector<list<pair<unsigned, unsigned>>> adjList; // adjList[u]= a list of pairs < color(u,v), v >, for each neighbour v of u
  vector<unsigned> sortedDegrees; // degrees of all nodes, in ascending order
  vector<unsigned> nodeColor; // color of each node
  vector<vector<unsigned>> edgeColor; // color of each edge

};


#endif // WL_GRAPH_HH
