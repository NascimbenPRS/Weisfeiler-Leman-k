#include "WL_Graph.hh"
#ifndef WL_ISO_HH
#define WL_ISO_HH


// WL isomorphism test wrapper class
class WL_IsoTester{
public:
    WL_IsoTester(){}; // empty Tester constructor
    bool WLK_Tuples(const WL_Graph& g1, const WL_Graph& g2, const unsigned k); // true iff WL-K with tuples distinguishes g1-g2
    bool WLK_Sets(const WL_Graph& g1, const WL_Graph& g2, const unsigned k); // true iff WL-K with sets ""
};


#endif // WL_ISO_HH
