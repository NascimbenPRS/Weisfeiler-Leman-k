#include "WL_Isomorph.hh"
#include <map>

// Auxiliary functions

// output: id of the tuple obtained by replacing the i-th element of tuple t with node w
unsigned getNeighbourId(unsigned n, unsigned tId, const vector<unsigned> t, unsigned w, unsigned i){
    return tId +(w- (int)t[i]) * pow (n,t.size()-1-i);
}

// output: adj. Matrix of subgraph of g1 induced by 'nodes'
const vector<vector<bool>> extractSubMat(const WL_Graph& g1, const vector<unsigned> nodes){
    unsigned n= nodes.size();
    vector<vector<bool>> adjMat1(n,vector<bool> (n,false));
    for (unsigned v1=0; v1 < n; v1++)
        for (unsigned v2=v1; v2 < n; v2++){
            if (g1.AdjMat(nodes[v1],nodes[v2])){
                adjMat1[v1][v2]= adjMat1[v2][v1] = true;
            }
        }
    return adjMat1;
}


//----------------------------------------------------------------------------------------

// Weisfeiler-Lehman Isomorphism Test (tuples)
bool WL_IsoTester::WLK_Tuples(const WL_Graph& g1,const WL_Graph& g2, const unsigned k){
    // -- auxiliary variables
    unsigned v, i;
    // 2. Isomorphism test
    bool maybeIsomorphic = (g1.Order() == g2.Order() && g1.Edges()==g2.Edges()
                            && g1.SortedDegrees() == g2.SortedDegrees()
                            ); // false iff G1 and G2 are certainly not isomorphic
    if (maybeIsomorphic){
      typedef unsigned t_color; // type color
      const unsigned n= g1.Order(), m=g1.Edges();
      t_color newColor= 0; // next color, increased whenever a class is split
      bool wasColorSplit= false; // true if currently examined class has been split
      bool isStableColouring = false;

      if (k==1){
        // COLOR REFINEMENT/WL-1

          vector<t_color> colorsG1(n,0), colorsG2(n,0); // colors for each node
          for (v=0; v< n; v++){
            colorsG1[v]= g1.NodeColor(v);
            colorsG2[v]= g2.NodeColor(v);
            if (colorsG1[v] > newColor) newColor= colorsG1[v]; // newColor= largest node color so far
            if (colorsG2[v] > newColor) newColor= colorsG2[v];
          }

          vector<tuple<t_color,multiset<pair<t_color,t_color>>,unsigned>> aggrG1(n), aggrG2(n); // aggregation maps <old c., c. multiset, node>
          // note: multiset contains pairs <c. edge (v,w), c. node w>
          multiset<pair<t_color,t_color>> empty_mset= {};

          vector<t_color> oldColorsG1 = colorsG1, oldColorsG2= colorsG2;
          while( (!isStableColouring) && maybeIsomorphic){ // run actual CR procedure
                  oldColorsG1 = colorsG1; oldColorsG2= colorsG2;
                  for(v=0; v < n; v++){
                      // reset aggrMap of node v and aggregate neighbouring colours
                      aggrG1[v]= make_tuple(oldColorsG1[v], empty_mset, v);
                      aggrG2[v]= make_tuple(oldColorsG2[v], empty_mset, v);
                      list<pair<unsigned,unsigned>> listG1= g1.AdjList(v),listG2= g2.AdjList(v);
                      for (auto iter= listG1.begin(); iter != listG1.end(); iter++)
                        get<1>(aggrG1[v]).insert(make_pair((*iter).first, colorsG1[(*iter).second]));
                      for (auto iter= listG2.begin(); iter != listG2.end(); iter++)
                        get<1>(aggrG2[v]).insert(make_pair((*iter).first, colorsG2[(*iter).second]));
                  } // end aggregation for one node
                  // end aggregation for all nodes

                  // Sort aggregation maps (by previous color and by multiset, in this order)
                  sort(aggrG1.begin(), aggrG1.end());
                  sort(aggrG2.begin(), aggrG2.end());

                  if (get<0>(aggrG1[0]) != get<0>(aggrG2[0]) || get<1>(aggrG1[0]) != get<1>(aggrG2[0]))
                      maybeIsomorphic= false; // // found a non-matching tuple-pair, G1 and G2 are not isomorphic
                  isStableColouring= true;
                  wasColorSplit= false;
                  for (v= 1; v < n && maybeIsomorphic; v++){
                      if (get<0>(aggrG1[v]) == get<0>(aggrG2[v]) && get<1>(aggrG1[v]) == get<1>(aggrG2[v])){ // matching previous colors and multisets
                          if (get<0>(aggrG1[v]) != get<0>(aggrG1[v-1])){ // scan has moved on to a block of nodes with a new previous color
                              wasColorSplit= false;
                          }
                          else {
                              if (get<1>(aggrG1[v]) != get<1>(aggrG1[v-1])){ // new multiset for same c. class -> split class, add new c.
                                  wasColorSplit= true;
                                  newColor++;
                                  isStableColouring= false; // a split has occurred, iterate again
                                  colorsG1[get<2>(aggrG1[v])] = newColor; // update color of corresponding nodes
                                  colorsG2[get<2>(aggrG2[v])] = newColor;
                              }
                              else{ // same colour and multiset as previous tuple-pair
                                  if (wasColorSplit){ // must assign this tuple-pair a new color
                                      colorsG1[get<2>(aggrG1[v])] = newColor; // update color of corresponding nodes
                                      colorsG2[get<2>(aggrG2[v])] = newColor;
                                  }
                              }
                          }
                      }
                      else maybeIsomorphic= false; // found a non-matching tuple-pair, G1 and G2 are not isomorphic
                  }
            }

      } // ----------     end: color refinement/WL-1

      else{  // --------  WL-K, K >= 2
        typedef vector<unsigned> t_atp; // type of an atp
        // ATP(v1,..,vk) = lower triangular half (diagonal excluded) of the matrix M on (v1,..,vk)^2 s.t.:
        //   M_t[vi,vj]= 2 if vi=vj, 1 if edge (vi,vj) exists, 0 else
        unsigned t, temp; // tuple indexes
        const unsigned numTuples= pow(n,k); // total # of tuples (n^k)
        vector<vector<unsigned>> tuples(numTuples,vector<unsigned>(k)); // for each index, corresponding k-tuple
        vector<t_color> colorsG1(numTuples), colorsG2(numTuples); // colors for each tuple
        vector<pair<t_atp, unsigned>> atpG1(numTuples), atpG2(numTuples); // atomic type; for each tuple t, a pair <atp of t,t>

        // 2.1 Compute tuples vector
        // !!! INEFFICIENT, MUST BE REVISED !!!
        unsigned weight; const unsigned maxWeight= numTuples / n;
        for (t=0; t < numTuples; t++){
            temp = t;
            weight= maxWeight; // fist digit: n^(k-1)
            for(unsigned i=0; i < k; i++){
                tuples[t][i]= temp / weight; // i-th node of tuple t
                temp %= weight;
                weight /= n;
            }
        }

        // 2.2 Compute atomic types
        unsigned atpLength= (k*(k-1))/2;
        for (t= 0; t < numTuples; t++){
          vector<unsigned> &tG= tuples[t]; // tuple associated to index t
          atpG1[t]= make_pair(t_atp(),t); // initialise atp of t with empty vector
          atpG2[t]= make_pair(t_atp(),t);
          for (i= 1; i < k; i++)
            for (unsigned j= 0; j < i; j++){ // scan lower half of the tuple's adj matrix, minus the diagonal
              // insert 2 in atp of t if t_i=t_j, 1 if edge (t_i,t_j) exist, 0 otherwise
              atpG1[t].first.push_back((tG[i]==tG[j]) ? 2 : (g1.AdjMat(tG[i],tG[j])));
              atpG2[t].first.push_back((tG[i]==tG[j]) ? 2 : (g2.AdjMat(tG[i],tG[j])));
            }
        }

        // 2.3 Sort atomic types, check if they match across G1 and G2 and map them into a canonical initial coloring
        sort(atpG1.begin(), atpG1.end());
        sort(atpG2.begin(), atpG2.end());
        if (atpG1[0].first == atpG2[0].first){ // matching atp, same color
            colorsG1[atpG1[0].second]= 0;
            colorsG2[atpG2[0].second]= 0;
        }
        else maybeIsomorphic= false;
        for (t= 1; (t < numTuples) && maybeIsomorphic; t++){
          if (atpG1[t].first == atpG2[t].first){ // matching atp
              if (atpG1[t].first != atpG1[t-1].first) newColor++; // atp differs from the previous one
              colorsG1[atpG1[t].second]= newColor; // update color of corresponding tuples
              colorsG2[atpG2[t].second]= newColor;
          }
          else maybeIsomorphic= false;
        }

        // 2.4 If atomic types match, run WL-k until mismatch or until a stable colouring is reached
        if (maybeIsomorphic){
            vector<t_color> oldColorsG1 = colorsG1, oldColorsG2= colorsG2;
            vector<t_color> neighColorsG1(k), neighColorsG2(k);
            vector<tuple<unsigned, multiset<vector<t_color>>, unsigned>> aggrG1(numTuples), aggrG2(numTuples); // : <old c., c. multiset, tuple index>
            multiset<vector<t_color>> empty_mset= {};
            while( (!isStableColouring) && maybeIsomorphic){  // run actual WL-k procedure
                oldColorsG1 = colorsG1; oldColorsG2= colorsG2;
                for(t=0; t < numTuples; t++){
                    // reset aggrMap of t
                    aggrG1[t]= make_tuple(oldColorsG1[t], empty_mset, t);
                    aggrG2[t]= make_tuple(oldColorsG2[t], empty_mset, t);
                    for (v= 0; v < n; v++){ // for each node v
                      for (i= 0; i < k; i++){ // for each index i
                          // compute color tuple of (v,i) neighbours of t
                          temp= getNeighbourId(n, t, tuples[t], v, i); // temp= t(v,i)
                          neighColorsG1[i]= oldColorsG1[temp];
                          neighColorsG2[i]= oldColorsG2[temp];
                      }
                    // add color tuple to color multiset of t
                    get<1>(aggrG1[t]).insert(neighColorsG1);
                    get<1>(aggrG2[t]).insert(neighColorsG2);
                    } // end aggregation for one tuple
                  } // end aggregation for all tuples

                  // Sort aggregation maps (by previous color and by multiset, in this order)
                  sort(aggrG1.begin(), aggrG1.end());
                  sort(aggrG2.begin(), aggrG2.end());

                  if (get<0>(aggrG1[0]) != get<0>(aggrG2[0]) || get<1>(aggrG1[0]) != get<1>(aggrG2[0]))
                      maybeIsomorphic= false; // found a non-matching tuple-pair, G1 and G2 are not isomorphic
                  isStableColouring= true;
                  wasColorSplit= false;
                  for (t= 1; t < numTuples && maybeIsomorphic; t++){
                      if (get<0>(aggrG1[t]) == get<0>(aggrG2[t]) && get<1>(aggrG1[t]) == get<1>(aggrG2[t])){ // matching previous colors and multisets
                          if (get<0>(aggrG1[t]) != get<0>(aggrG1[t-1])){ // scan has moved on to a block of tuple-pairs with a new previous color
                              wasColorSplit= false;
                          }
                          else {
                              if (get<1>(aggrG1[t]) != get<1>(aggrG1[t-1])){ // new multiset for same c. class -> split class, add new c.
                                  wasColorSplit= true;
                                  newColor++;
                                  isStableColouring= false; // a split has occurred, iterate again
                                  colorsG1[get<2>(aggrG1[t])] = newColor; // update color of tuple-pair
                                  colorsG2[get<2>(aggrG2[t])] = newColor;
                              }
                              else{ // same colour and multiset as previous tuple-pair
                                  if (wasColorSplit){ // must assign this tuple-pair a new color
                                      colorsG1[get<2>(aggrG1[t])] = newColor; // update color of corresponding tuples
                                      colorsG2[get<2>(aggrG2[t])] = newColor;
                                  }
                              }
                          }
                      }
                      else maybeIsomorphic= false; // found a non-matching tuple-pair, G1 and G2 are not isomorphic
                  } // end for loop on tuples
            } // end while loop (WL iterations)
        } // end if (maybeIsomorphic)
      } // end else (k >= 2)

    } // end if (n1 == n2, etc...)
    return (maybeIsomorphic);
};



// Weisfeiler-Lehman Isomorphism Test (sets)
bool WL_IsoTester::WLK_Sets(const WL_Graph& g1,const WL_Graph& g2, const unsigned k){
    // -- auxiliary variables
    unsigned i;
    // 2. Isomorphism test
    bool maybeIsomorphic = (g1.Order() == g2.Order() && g1.Edges()==g2.Edges()
                            && g1.SortedDegrees() == g2.SortedDegrees()
                            ); // false iff G1 and G2 are certainly not isomorphic
    if (maybeIsomorphic){
      typedef unsigned t_color; // type color
      const unsigned n= g1.Order(), m=g1.Edges();
      t_color newColor= 0; // next color, increased whenever a new graph iso. class is found

      typedef set<unsigned> t_nodeSet; // type: set of nodes
      typedef vector<vector<bool>> t_adjMat; // type: adjacency matrix
      map<t_nodeSet, unsigned> setsMap; // <set s, index s>
      vector <t_nodeSet> setsVect; // [s]= set s
      vector <t_color> setsColorsG1, setsColorsG2; // color of k-sets
      vector <pair<t_adjMat, unsigned>> kMatG1, kMatG2; // sub-matrices of size k: <submat of set s, index s>
      vector <pair<WL_Graph, t_color>> isoClasses; // a specimen for each different iso. class: <graph, color>

      // 2.1 Compute all k-subsets and sub-matrices
      unsigned numTuples= pow(n,k), t=0, temp=0;
      unsigned weight; unsigned maxWeight= numTuples / n;
      unsigned s=0;
      t_nodeSet temp_set;
      vector<unsigned> temp_setVect(k,0);
      for (t=0; t < numTuples; t++){
        temp = t;
        weight= maxWeight; // fist digit: n^(k-1)
        temp_set.clear();
        for(i= 0; i < k; i++){
                temp_setVect[i]= temp / weight; // i-th node of set s
                temp_set.insert(temp_setVect[i]);
                temp %= weight;
                weight /= n;
        }
        if (temp_set.size() == k && setsMap.find(temp_set) == setsMap.end()){
            // new k-set found
            setsMap.insert({temp_set, s});
            setsVect.push_back(temp_set);
            // insert adj. matrices of subgraphs induced by set s on graphs G1-2
            kMatG1.push_back(make_pair(extractSubMat(g1,temp_setVect),s));
            kMatG2.push_back(make_pair(extractSubMat(g2,temp_setVect),s));
            s++;
        }
      }
      // Color k-subgraphs
      sort(kMatG1.begin(),kMatG1.end());
      sort(kMatG2.begin(),kMatG2.end());
      WL_Graph temp_graph;
      WL_IsoTester isoTest;
      unsigned numSets= s; // number of sets
      bool foundIso= false; // true iff found a saved graph which is iso. to the current one
      t_color lastColorG1= 0, lastColorG2 =0;
      setsColorsG1.resize(numSets); setsColorsG2.resize(numSets);
      // color first k-subset-pair
      isoClasses.push_back(make_pair(WL_Graph(kMatG1[0].first), 0)); // save first graph, with c. 0
      setsColorsG1[kMatG1[0].second]= 0;
      if (kMatG1[0].first == kMatG2[0].first)
        setsColorsG2[kMatG2[0].second]= 0; // first sub-mat. are equal -> same color
      else{ // different sub-mat., check iso.
        temp_graph= WL_Graph(kMatG2[0].first); // make graph from sub-mat.
        for(auto iter= isoClasses.begin(); !foundIso && iter != isoClasses.end(); iter++){
            // check temp_graph against previously saved graphs
            maybeIsomorphic= true;
            for (unsigned k1= 1; k1 < k && maybeIsomorphic; k1++){
              // run WL-1, WL-2.. until k-1 OR confirmed non-isomorphism
              if (!isoTest.WLK_Tuples(temp_graph, (*iter).first, k1))
                maybeIsomorphic= false; // definitely not iso.
            }
            if (maybeIsomorphic){
                // definitely iso., assign same color
                foundIso= true;
                setsColorsG2[kMatG2[0].second]= 0;
                lastColorG2 = 0;
            }
        } // end for (iso classes)
        if (!foundIso){
            // No iso. graph found, save graph with new color
            newColor++;
            isoClasses.push_back(make_pair(temp_graph, newColor));
            setsColorsG2[kMatG2[0].second]= newColor;
            lastColorG2= newColor;
        }
      }

      // color remaining k-subset-pairs
      bool eqMatG1= false, eqMatG2= false;
      for (s= 1; s < numSets; s++){
          eqMatG1= eqMatG2= false;
          if ((kMatG1[s].first == kMatG1[s-1].first)){
            // G1: same as previous matrix -> same color
            setsColorsG1[kMatG1[s].second]= lastColorG1;
            eqMatG1= true;
          }
          if ((kMatG2[s].first == kMatG2[s-1].first)){
            // G2: same as previous matrix -> same color
            setsColorsG2[kMatG2[s].second]= lastColorG2;
            eqMatG2= true;
          }

          if (!eqMatG1){
            // G1: different sub-mat., check iso.
            temp_graph= WL_Graph(kMatG1[s].first); // make graph from sub-mat.
            foundIso= false;
            for(auto iter= isoClasses.begin(); !foundIso && iter != isoClasses.end(); iter++){
                // check temp_graph against previously saved graphs
                maybeIsomorphic= true;
                for (unsigned k1= 1; k1 < k && maybeIsomorphic; k1++){
                  // run WL-1, WL-2.. until k-1 OR confirmed non-isomorphism
                  if (!isoTest.WLK_Tuples(temp_graph, (*iter).first, k1))
                    maybeIsomorphic= false; // definitely not iso.
                }
                if (maybeIsomorphic){
                    // definitely iso., assign same color
                    foundIso= true;
                    setsColorsG1[kMatG1[s].second]= (*iter).second;
                    lastColorG1 = (*iter).second;
                }
            } // end for (iso classes)
            if (!foundIso){
                // No iso. graph found, save graph with new color
                newColor++;
                isoClasses.push_back(make_pair(temp_graph, newColor));
                setsColorsG1[kMatG1[s].second]= newColor;
                lastColorG1= newColor;
            }
          }
          if (!eqMatG2){
            // G2: different sub-mat., check iso.
            temp_graph= WL_Graph(kMatG2[s].first); // make graph from sub-mat.
            foundIso= false;
            for(auto iter= isoClasses.begin(); !foundIso && iter != isoClasses.end(); iter++){
                // check temp_graph against previously saved graphs
                maybeIsomorphic= true;
                for (unsigned k1= 1; k1 < k && maybeIsomorphic; k1++){
                  // run WL-1, WL-2.. until k-1 OR confirmed non-isomorphism
                  if (!isoTest.WLK_Tuples(temp_graph, (*iter).first, k1))
                    maybeIsomorphic= false; // definitely not iso.
                }
                if (maybeIsomorphic){
                    // definitely iso., assign same color
                    foundIso= true;
                    setsColorsG2[kMatG2[s].second]= (*iter).second;
                    lastColorG2 = (*iter).second;
                }
            } // end for (iso classes)
            if (!foundIso){
                // No iso. graph found, save graph with new color
                newColor++;
                isoClasses.push_back(make_pair(temp_graph, newColor));
                setsColorsG2[kMatG2[s].second]= newColor;
                lastColorG2= newColor;
            }
          }
      } // end: k-sets coloring


      // All k-sets are now colored
      // Define edges/k+1 sets and their colors in the new graphs
      // 2.1 Compute all k+1-subsets and sub-matrices
      unsigned oldNumSets= s;
      s=0;
      maxWeight = pow(n,k);; // most significant digit has weight n^k
      numTuples = pow(n,k+1); // n^k+1 tuples
      temp_setVect.resize(k+1,0);
      temp_set.clear();
      setsMap.clear();
      kMatG1.clear(); kMatG2.clear();
      vector <t_color> edgesColorsG1, edgesColorsG2; // color of k+1-sets
      for (t=0; t < numTuples; t++){
        temp = t;
        weight= maxWeight; // fist digit: n^k
        temp_set.clear();
        for(i= 0; i < k+1; i++){
                temp_setVect[i]= temp / weight; // i-th node of set s
                temp_set.insert(temp_setVect[i]);
                temp %= weight;
                weight /= n;
        }
        if (temp_set.size() == k+1 && setsMap.find(temp_set) == setsMap.end()){
            // new k+1-set found
            setsMap.insert({temp_set, s});
            // insert adj. matrices of subgraphs induced by set s on graphs G1-2
            kMatG1.push_back(make_pair(extractSubMat(g1,temp_setVect),s));
            kMatG2.push_back(make_pair(extractSubMat(g2,temp_setVect),s));
            s++;
        }
      }

      // Color k+1-subgraphs
      sort(kMatG1.begin(),kMatG1.end());
      sort(kMatG2.begin(),kMatG2.end());
      numSets= s; // number of k+1-sets
      foundIso= false; // true iff found a saved graph which is iso. to the current one
      lastColorG1= lastColorG2= newColor= 0;
      edgesColorsG1.resize(numSets); edgesColorsG2.resize(numSets);
      // color first k+1-subset-pair
      isoClasses.clear();
      isoClasses.push_back(make_pair(WL_Graph(kMatG1[0].first), 0)); // save first graph, with c. 0
      edgesColorsG1[kMatG1[0].second]= 0;
      if ((kMatG1[0].first == kMatG2[0].first))
        edgesColorsG2[kMatG2[0].second]= 0; // first sub-mat. are equal -> same color
      else{ // different sub-mat., check iso.
        temp_graph= WL_Graph(kMatG2[0].first); // make graph from sub-mat.
        for(auto iter= isoClasses.begin(); !foundIso && iter != isoClasses.end(); iter++){
            // check temp_graph against previously saved graphs
            maybeIsomorphic= true;
            for (unsigned k1= 1; k1 < k+1 && maybeIsomorphic; k1++){
              // run WL-1, WL-2.. until k OR confirmed non-isomorphism
              if (!isoTest.WLK_Tuples(temp_graph, (*iter).first, k1))
                maybeIsomorphic= false; // definitely not iso.
            }
            if (maybeIsomorphic){
                // definitely iso., assign same color
                foundIso= true;
                edgesColorsG2[kMatG2[0].second]= (*iter).second;
                lastColorG2 = (*iter).second;
            }
        } // end for (iso classes)
        if (!foundIso){
            // No iso. graph found, save graph with new color
            newColor++;
            isoClasses.push_back(make_pair(temp_graph, newColor));
            edgesColorsG2[kMatG2[0].second]= newColor;
            lastColorG2= newColor;
        }
      }

      // color remaining k+1-subset-pairs
      eqMatG1= eqMatG2= false;
      for (s= 1; s < numSets; s++){
          eqMatG1= eqMatG2= false;
          if ((kMatG1[s].first == kMatG1[s-1].first)){
            // G1: same as previous matrix -> same color
            edgesColorsG1[kMatG1[s].second]= lastColorG1;
            eqMatG1= true;
          }
          if ((kMatG2[s].first == kMatG2[s-1].first)){
            // G2: same as previous matrix -> same color
            edgesColorsG2[kMatG2[s].second]= lastColorG2;
            eqMatG2= true;
          }

          if (!eqMatG1){
            // G1: different sub-mat., check iso.
            temp_graph= WL_Graph(kMatG1[s].first); // make graph from sub-mat.
            foundIso= false;
            for(auto iter= isoClasses.begin(); !foundIso && iter != isoClasses.end(); iter++){
                // check temp_graph against previously saved graphs
                maybeIsomorphic= true;
                for (unsigned k1= 1; k1 < k+1 && maybeIsomorphic; k1++){
                  // run WL-1, WL-2.. until k OR confirmed non-isomorphism
                  if (!isoTest.WLK_Tuples(temp_graph, (*iter).first, k1))
                    maybeIsomorphic= false; // definitely nsecondot iso.
                }
                if (maybeIsomorphic){
                    // definitely iso., assign same color
                    foundIso= true;
                    edgesColorsG1[kMatG1[s].second]= (*iter).second;
                    lastColorG1 = (*iter).second;
                }
            } // end for (iso classes)
            if (!foundIso){
                // No iso. graph found, save graph with new color
                newColor++;
                isoClasses.push_back(make_pair(temp_graph, newColor));
                edgesColorsG1[kMatG1[s].second]= newColor;
                lastColorG1= newColor;
            }
          }
          if (!eqMatG2){
            // G2: different sub-mat., check iso.
            temp_graph= WL_Graph(kMatG2[s].first); // make graph from sub-mat.
            foundIso= false;
            for(auto iter= isoClasses.begin(); !foundIso && iter != isoClasses.end(); iter++){
                // check temp_graph against previously saved graphs
                maybeIsomorphic= true;
                for (unsigned k1= 1; k1 < k+1 && maybeIsomorphic; k1++){
                  // run WL-1, WL-2.. until k OR confirmed non-isomorphism
                  if (!isoTest.WLK_Tuples(temp_graph, (*iter).first, k1))
                    maybeIsomorphic= false; // definitely not iso.
                }
                if (maybeIsomorphic){
                    // definitely iso., assign same color
                    foundIso= true;
                    edgesColorsG2[kMatG2[s].second]= (*iter).second;
                    lastColorG2 = (*iter).second;
                }
            } // end for (iso classes)
            if (!foundIso){
                // No iso. graph found, save graph with new color
                newColor++;
                isoClasses.push_back(make_pair(temp_graph, newColor));
                edgesColorsG2[kMatG2[s].second]= newColor;
                lastColorG2= newColor;
            }
          }
      } // end: k+1-sets coloring

      // -- free memory
      kMatG1.clear(); kMatG2.clear();
      isoClasses.clear();

      // -- k+1 subsets are colored, stored in edgesColorsG1/G2
      t_adjMat newAdjMat (oldNumSets, vector<bool> (oldNumSets, false));
      vector<list<pair<t_color,t_color>>> newAdjList(oldNumSets);
      t_nodeSet temp_set2;
      // Create new graph G1'
      for (unsigned s1=0; s1 < oldNumSets; s1++){
        for(unsigned s2=s1+1; s2 < oldNumSets; s2++){
            // merge s1-s2
            temp_set= setsVect[s1];
            temp_set2= setsVect[s2];
            for(auto iter= temp_set2.begin(); iter != temp_set2.end(); iter++)
                temp_set.insert(*iter);
            if (temp_set.size() == k + 1){
                    // sets s1-s2 differ for exactly one element -> neighbors
                    newAdjList[s1].push_back(make_pair(edgesColorsG1[setsMap[temp_set]],s2));
                    newAdjList[s2].push_back(make_pair(edgesColorsG1[setsMap[temp_set]],s1));
            }
        }
      }

      // -- free memory
      edgesColorsG1.clear();

      WL_Graph newG1(setsColorsG1, newAdjList);

      setsColorsG1.clear();
      newAdjList.clear();
      newAdjList.resize(oldNumSets);
      // Create new graph G2'
      for (unsigned s1=0; s1 < oldNumSets; s1++){
        for(unsigned s2=s1+1; s2 < oldNumSets; s2++){
            // merge s1-s2
            temp_set= setsVect[s1];
            temp_set2= setsVect[s2];
            for(auto iter= temp_set2.begin(); iter != temp_set2.end(); iter++)
                temp_set.insert(*iter);
            if (temp_set.size() == k + 1){
                    // sets s1-s2 differ for exactly one element -> neighbors
                    newAdjList[s1].push_back(make_pair(edgesColorsG2[setsMap[temp_set]],s2));
                    newAdjList[s2].push_back(make_pair(edgesColorsG2[setsMap[temp_set]],s1));
            }
        }
      }

      // -- free memory
      setsMap.clear(); setsVect.clear();
      edgesColorsG2.clear();

      WL_Graph newG2(setsColorsG2, newAdjList);
      // -- free memory
      setsColorsG2.clear();
      // -- run CR on new graphs
      maybeIsomorphic= isoTest.WLK_Tuples(newG1, newG2, 1);


    } // end if (maybeIsomorphic)

    return (maybeIsomorphic);
};

