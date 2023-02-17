/******************************************************************************
				A graph isomorphism test based on the K-Weisfeiler Lehman
        paritioning/coloring algorithm
			 
        Francesco Nascimben 2023
			 
-- Input (from command-line):
--- 1. first graph instance in standard DIMACS format (G1), 
--- 2. second graph instance in standard DIMACS format (G2), 
--- 3. size of examined tuples (k >= 1)
--- 4. "sets"/"tuples", to run set/tuple-based WL-k

-- Output: false if WL-k declares G1-G2 to be not isomorphic, true otherwise.
*******************************************************************************/

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
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <math.h>	/* ceil */

#include<cstdio>
#include<ctime>

#include "WL_Graph.hh"
#include "WL_Isomorph.hh"


using namespace std;

int main(int argc, char *argv[])
{
  WL_Graph g1= WL_Graph(argv[1]);
  WL_Graph g2= WL_Graph(argv[2]);
  unsigned k= stoi(argv[3]);
  
  WL_IsoTester tester;
 
  
  bool maybeIsomorphic= false;
  const string WL_type = argv[4];
  if (WL_type == "sets") maybeIsomorphic= tester.WLK_Sets(g1,g2,k);
  else if (WL_type == "tuples") maybeIsomorphic= tester.WLK_Tuples(g1,g2,k);
       else {
        cerr << "Error: fourth parameter must be string \"tuples\" or string \"sets\"" << endl;
        exit(1);
       }

  if (maybeIsomorphic) cout << "Maybe isomorphic";
  else                 cout << "NOT isomorphic";
  /*
  // Uncomment this section to run both versions of WL-k on the same instance graph-pair
  string s1= argv[1], s2= argv[2];
  bool maybeIsomorphic1, maybeIsomorphic2= false;
  if (s1 == s2)
    cout << "SAME FILE: ISOMORPHIC" << endl;
  else{
    if (s1 < s2){ // for batch-testing only, avoids test repetitions
    maybeIsomorphic1= tester.WLK_Sets(g1,g2,k);
    maybeIsomorphic2= tester.WLK_Tuples(g1,g2,k);
    if (maybeIsomorphic1 == maybeIsomorphic2)
      cout << "SAME RESULT: " << ((maybeIsomorphic1) ? "maybe iso" : "not iso") << endl;
    else
      cout << "DIFFERENT RESULT: sets " << ((maybeIsomorphic1) ? "maybe iso" : "not iso") << 
        ", tuples " << ((maybeIsomorphic2) ? "maybe iso" : "not iso") << endl;
    }
  }
  */

    
}   
