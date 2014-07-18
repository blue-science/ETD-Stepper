#ifndef MODEL_G_H
#define MODEL_G_H
// Author: Matt Pulver <matt@blue-science.org>
// Date: 2014 July

#include "etdrk4b.h"

class ModelG : public ETDRK4b {
  protected:
    enum { G=0, X, Y }; // reactants. Use as indices for readability as needed.

    Real a, b, dx, dy, p, q, g, s, u, w, G0, X0, Y0;

    void read_config();

    // Called by algorithm constructor since phi_init() depends on it.
    void init_dcoeff();

    // Set t=0 and rv
    void set_initial_conditions(); // Called by ExtTimeDiff::run().

    // Calculate reaction, and FFT-normalize by dividing by N.
    void reaction_normalize_n( int start, int end );

  public:

    ModelG(const char* filename);

    ~ModelG();
};

#endif // #ifndef MODEL_G_H
