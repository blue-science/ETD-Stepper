// ExpTimeDiff.h
// Author: Matt Pulver <matt@blue-science.org>
// Date: 2014 July
#ifndef EXPTIMEDIFF_H
#define EXPTIMEDIFF_H

#include <complex>
#include <cstdio> // rename()
#include <ctime> // clock()
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <signal.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread> // C++11
#include <vector>

#include <boost/filesystem.hpp>
// Compress large output data
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

using namespace std;

typedef double Real;
typedef complex<Real> Complex;

class ExpTimeDiff {
  protected:
    int nthreads;
    ifstream config;
    int nreacts; // Number of reactants
    int nstages; // Number of R-K stages

    double T, // Total time, starting at t=0
           dt, // Computational unit of time
           t; // Current time.
    int fput; // Frames per unit time. Each frame is saved to a file under frames/.

    vector<Real> L; // Length of each spatial dimension
    int *dims; // Number of computational cells along each dimension
    int rank, // Number of spatial dimensions
        N, // Number of spatial computation points = product of dims
        Nh, // Nnumber of elements in Hermitian-reduced matrices Nh < N
        cn; // Index of center point

    // For the following vectors/arrays:
    //  * 0 <= i  <= nstages. i=0 is for the FFT-transformed reactants. 1..nstages are for the stages.
    //  * 0 <= r  <  nreacts
    //  * 0 <= n  <  N
    //  * 0 <= nh <  Nh
    vector<Real> dcoeff; // Diffusion coefficients dcoeff[r]
    vector<Real*> rv; // Vector of reactant arrays. rv[r][n] where r is reactant and 0 <= n < N.
    vector<Real*> rv_output; // Copy of rv for multi-threaded output.
    vector<vector<Complex*> > k; // k[i][r][nh] where i=0 is transformed reactant. 0<i are the R-K values.
    vector<vector<fftw_plan> > planc2r, planr2c; // planc2r[i][r]

    int litmus; // Reaction index for litmus testing

    static bool save_and_exit; // Static since sighandler must be static

    ExpTimeDiff( const char* config_filename, int nreacts, int nstages );
    ~ExpTimeDiff();

    void read_config();

    // dims[0] * dims[1] * ... * dims[rank-2] * (dims[rank-1]/2+1)
    int getNh() const;

    int center() const;

    // Given an index n in a hermitian-reduced matrix where 0 <= n < Nh,
    // return the index it "would have been" in a non-reduced matrix with 0 <= n < N.
    int hermitianMap( int n ) const;

    // Has numerical errors as z -> 0.
    template <typename T> // T is Real or Complex
    static T phi_with_instability( int n, T z )
    {
        T sum = 1.0;
        for( Real k = n-1 ; 0 < k ; --k ) sum = 1.0 + z*sum/k;
        return n == 0 ? exp(z) : (exp(z) - sum) / pow(z,n);
    }

    // Good for |z| << 1, otherwise inaccurate.
    static Real phi_maclaurin( int n, Real z );

    // Good for all z but slow.
    static Real phi_contour( int n, Real z );

    // Return phi_n(z). phi_0(z) = exp(z).
    static Real phi_calc( int n, Real z );

    void normalize_rv();

    static void sighandler( int signum );

    // All double values:
    // t
    // Gval0 Xval0 Yval0
    // Gval1 Xval1 Yval1
    // ...
    // Gval(N-1) Xval(N-1) Yval(N-1)
    void dump_state();

    // Sets t and rv[][] if successful
    bool import_state();

    // t is double-precision. All other numbers are single-precision floating point in binary format.
    // t
    // rank L0 dim0 L1 dim1 ... L(rank-1) dim(rank-1)
    // Gval0 Xval0 Yval0
    // Gval1 Xval1 Yval1
    // ...
    // Gval(N-1) Xval(N-1) Yval(N-1)
    // 
    // where N = product of dims
    void output_frame( double t );

    static void write_float( boost::iostreams::filtering_ostream *pout, float f );

    // Set by algorithm class.
    virtual void calc_k( int ki, int start, int end ) = 0;

    // Set by model class.
    virtual void set_initial_conditions() = 0;

    // Set by model class.
    // Calculate reaction, and FFT-normalize rv by dividing by N.
    // For multi-threading, specify starting index of rv and length.
    virtual void reaction_normalize_n( int start, int n ) = 0;

    // Set by model class.
    virtual void init_dcoeff() = 0;

  public:

    void run();
};

#endif // #ifndef EXPTIMEDIFF_H
