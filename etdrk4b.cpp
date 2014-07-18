// etdrk4b.cpp
// Author: Matt Pulver <matt@blue-science.org>
// Date: 2014 July

#include "etdrk4b.h"

ETDRK4b::ETDRK4b( const char* config_filename, int nreacts ):
    ExpTimeDiff(config_filename,nreacts,4) { } // 4-stage Runge Kutta

void ETDRK4b::phi_init_thread( int start, int end )
{
    int n, m, i, r;
    div_t q;
    Real nom2dt;
    clock_t next_update = clock();
    //for( n = 0 ; n < Nh ; ++n )
    for( n = start ; n != end ; ++n )
    {
        //if( (n&0x7fff) == 0 ) cerr << "In phi_init(). n = " << n << endl;
        if( start == 0 && next_update <= clock() )
        {
            next_update += nthreads * CLOCKS_PER_SEC;
            cerr << "phi_init(): " << 100.0*n/end << '%' << endl;
        }
        q.quot = hermitianMap(n);
        nom2dt = 0;
        for( i = rank-1 ; 0 <= i ; --i )
        {
            m = dims[i];
            q = div( q.quot, m );
            nom2dt += pow( ((q.rem+m/2)%m-m/2)/L[i], 2 );
        }
        nom2dt *= -dt * (4.0*M_PI*M_PI);
        for( i = 0 ; i < 4 ; ++i ) for( r = 0 ; r < nreacts ; ++r )
        {
            phi[i][r][n] = phi_calc(i,nom2dt*dcoeff[r]);
        }
        for( i = 0 ; i < 3 ; ++i ) for( r = 0 ; r < nreacts ; ++r )
        {
            phih[i][r][n] = phi_calc(i,nom2dt*dcoeff[r]/2);
        }
    }
}

// Initialize phi and phih for use in calc_k.
// These are just cached values to speed up runtime calculation.
void ETDRK4b::phi_init()
{
    int total, n;
    div_t q = { 0, 0 };
    thread threads[nthreads];
    phi.resize(4,vector<vector<Real> >(nreacts,vector<Real>(Nh)));
    phih.resize(3,vector<vector<Real> >(nreacts,vector<Real>(Nh)));

    for( total = 0, n = 0 ; n < nthreads ; total += q.quot, ++n )
    {
        q = div( Nh + q.rem, nthreads );
        threads[n] = thread( &ETDRK4b::phi_init_thread, this, total, total+q.quot );
    }
    for( thread &th : threads ) th.join();
}

void ETDRK4b::calc_k( int ki, int start, int end )
{
    int r, n;
    switch(ki)
    {
      case 0: for( r = 0 ; r < nreacts ; ++r ) for( n = start ; n != end ; ++n )
        k[0][r][n] = phi[0][r][n]*k[0][r][n] + dt * (
            ( phi[1][r][n] - 3*phi[2][r][n] + 4*phi[3][r][n] ) * k[1][r][n] +
            ( 2*phi[2][r][n] - 4*phi[3][r][n] ) * ( k[2][r][n] + k[3][r][n] ) +
            (  -phi[2][r][n] + 4*phi[3][r][n] ) * k[4][r][n] );
        break;
      case 1: for( r = 0 ; r < nreacts ; ++r ) for( n = start ; n != end ; ++n )
        k[1][r][n] = k[0][r][n];
        break;
      case 2: for( r = 0 ; r < nreacts ; ++r ) for( n = start ; n != end ; ++n )
        k[2][r][n] = phih[0][r][n]*k[0][r][n] + dt*phih[1][r][n]/2*k[1][r][n];
        break;
      case 3: for( r = 0 ; r < nreacts ; ++r ) for( n = start ; n != end ; ++n )
        k[3][r][n] = phih[0][r][n]*k[0][r][n] +
            dt*( (phih[1][r][n]/2 - phih[2][r][n])*k[1][r][n] + phih[2][r][n]*k[2][r][n] );
        break;
      case 4: for( r = 0 ; r < nreacts ; ++r ) for( n = start ; n != end ; ++n )
        k[4][r][n] = phi[0][r][n]*k[0][r][n] +
            dt*( (phi[1][r][n] - 2*phi[2][r][n])*k[1][r][n] + 2*phi[2][r][n]*k[3][r][n] );
    }
}
