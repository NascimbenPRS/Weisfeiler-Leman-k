// File WL_Graph.cc
#include "WL_Graph.hh"
#include <fstream>

WL_Graph::WL_Graph(string file_name1)
{
    // constructor from DIMACS file
    // 1. Parse instance input file
      const unsigned MAX_DIM = 1000;
      unsigned v1, v2; // temporary nodes
      char ch, buffer[MAX_DIM];

      // 1.1 Read graph
      ifstream is(file_name1);
      if(!is)
      {
        cerr << "Cannot open input file " <<  file_name1 << endl;
        exit(1);
      }

      is >> ch; // read first character
      while (ch == 'c'){ // comment line
        is.ignore(MAX_DIM,'\n'); // skip line
        is >> ch; // read first character of next line
      }
      is >> buffer >> n; // read number of nodes
      is >> m; // read number of edges
      // -- initialise data structures
      isColored= false;
      adjList.resize(n, list<pair<unsigned,unsigned>>());
      adjMat.resize(n, vector<bool>(n,false));
      sortedDegrees.resize(n,0);
      nodeColor.resize(n,0);
      edgeColor.resize(n, vector<unsigned> (n,0));

      // -- read nodes/edges
      while (is >> ch){ // 'e' found
          // read line (edge v1-v2) and populate adjacency lists and matrix
          is >> v1 >> v2;
          adjList[v1-1].push_back(make_pair(0,v2-1));
          adjList[v2-1].push_back(make_pair(0,v1-1));
          adjMat[v1-1][v2-1]= true;
          adjMat[v2-1][v1-1]= true;
          sortedDegrees[v1-1]++;
          sortedDegrees[v2-1]++;
          // !!!  MUST IMPLEMENT READING OF NODE AND EDGE COLORS  !!!!!


          // !!!
      }

      sort(sortedDegrees.begin(), sortedDegrees.end());
}

WL_Graph::WL_Graph(const WL_Graph& g1, const vector<unsigned> nodes)
{
    // constructor from graph G1 (subgraph induced by 'nodes')
    n= nodes.size();
    isColored= g1.isColored;
    adjMat.resize(n,vector<bool> (n,false));
    adjList.resize(n, list<pair<unsigned,unsigned>>());
    sortedDegrees.resize(n,0);
    nodeColor.resize(n,0);
    edgeColor.resize(n, vector<unsigned> (n,0));
    for (unsigned v1=0; v1 < nodes.size(); v1++)
        for (unsigned v2=v1; v2 < nodes.size(); v2++){
            if (g1.adjMat[nodes[v1]][nodes[v2]]){
                adjMat[v1][v2]= adjMat[v2][v1] = true;
                adjList[v1].push_back(make_pair(0,v2));
                adjList[v2].push_back(make_pair(0,v1));
                sortedDegrees[v1]++;
                sortedDegrees[v2]++;
                m++;
            }
        }
    sort(sortedDegrees.begin(), sortedDegrees.end());
    // MUST IMPLEMENT COLORS
}

WL_Graph::WL_Graph(const vector<vector<bool>> adjMat1)
{
    // graph constructor from adjacency matrix (uncoloured graph)
    n= adjMat1.size();
    m= 0;
    isColored= false;
    adjMat= adjMat1;
    adjList.resize(n, list<pair<unsigned,unsigned>>());
    sortedDegrees.resize(n,0);
    nodeColor.resize(n,0);
    edgeColor.resize(n, vector<unsigned> (n,0));
    for (unsigned v1=0; v1 < n; v1++)
        for (unsigned v2=v1; v2 < n; v2++){
            if (adjMat[v1][v2]){
                adjList[v1].push_back(make_pair(0,v2));
                adjList[v2].push_back(make_pair(0,v1));
                sortedDegrees[v1]++;
                sortedDegrees[v2]++;
                m++;
            }
        }
     sort(sortedDegrees.begin(), sortedDegrees.end());
}

WL_Graph::WL_Graph(const vector<unsigned> nodeColor1, const vector<list<pair<unsigned,unsigned>>> adjList1)
{
   // constructor from node colors vector and colored adjList (colored g.)
   n= nodeColor1.size();
   m= 0;
   isColored= true;
   adjMat.clear(); // adjMat of the resulting graph is empty, to minimise memory usage
   nodeColor= nodeColor1;
   adjList= adjList1;
   sortedDegrees.clear();
   // -- edgeColor
}

