#ifndef ETDRK4B_H
#define ETDRK4B_H
// Author: Matt Pulver <matt@blue-science.org>
// Date: 2014 July

#include "ExpTimeDiff.h"

class ETDRK4b : public ExpTimeDiff {
  protected:
    vector<vector<vector<Real> > > phi, phih;

    ETDRK4b( const char* config_filename, int nreacts );

    void phi_init_thread( int start, int end );

    // Must be called in derived class since it depends on dcoeff.
    void phi_init();

    void calc_k( int ki, int start, int end );
};

#endif // #ifndef ETDRK4B_H
