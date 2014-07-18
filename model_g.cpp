#include "model_g.h"
// Author: Matt Pulver <matt@blue-science.org>
// Date: 2014 July

ModelG::ModelG(const char* config_filename):ETDRK4b(config_filename,3) // 3 active reactants: G, X, Y
{
    read_config();
    config.close(); // opened by ExpTimeDiff()
    init_dcoeff();
    phi_init(); // depends on init_dcoeff() which depends on read_config().
    litmus = Y;
}

ModelG::~ModelG() { }

// Read from config ifstream opened and read by ExpTimeDiff constructor.
void ModelG::read_config()
{
    string line, name;
    stringstream ss;
    int i;

    do { getline( config, line ); } while( line[0] == '#' );
    ss.clear();
    ss.str(line);
    while( !ss.eof() )
    {
        ss >> line;
        name = line.substr( 0, i=line.find('=') );
        if( name == "a" ) a = stod(line.substr(i+1));
        else if( name == "b" ) b = stod(line.substr(i+1));
        else if( name == "dx" ) dx = stod(line.substr(i+1));
        else if( name == "dy" ) dy = stod(line.substr(i+1));
        else if( name == "p" ) p = stod(line.substr(i+1));
        else if( name == "q" ) q = stod(line.substr(i+1));
        else if( name == "g" ) g = stod(line.substr(i+1));
        else if( name == "s" ) s = stod(line.substr(i+1));
        else if( name == "u" ) u = stod(line.substr(i+1));
        else if( name == "w" ) w = stod(line.substr(i+1));
    }
    G0 = (   a + g*w ) / ( q - g*p );
    X0 = ( p*a + q*w ) / ( q - g*p );
    Y0 = ( s*X0*X0 + b )*X0 / ( X0*X0 + u );
}

void ModelG::init_dcoeff()
{
    dcoeff.resize(nreacts);
    dcoeff[G] = 1;
    dcoeff[X] = dx;
    dcoeff[Y] = dy;
}

// Set t=0 and rv
void ModelG::set_initial_conditions()
{
    const Real alpha = 3875.0/4096;
    vector<Real> v(rank); // position vector for plotting initial conditions
//    c0 { 0, 0, 0 };
//    c1 { 0,  2,  2/sqrt(3) }, // center of gaussian 1
//    c2 { 0, -2,  2/sqrt(3) }, // center of gaussian 2
//    c3 { 0,  0, -4/sqrt(3) }; // center of gaussian 3
    const double len = 12; // Length of segment along z-axis
    double d2, // distance squared
           dz; // z-distance from end points, 0 if within enpoints
    div_t q;
    int i, j, n;
    t = 0;
    for( i = 0 ; i < N ; ++i )
    {
        q.quot = i;
        for( j = rank-1 ; 0 <= j ; --j )
        {
            n = dims[j];
            q = div( q.quot, n );
            v[j] = ( q.rem - n/2 ) * L[j] / n; // -L[j]/2 <= v[j] < L[j]/2
        }
        //dp2 = v[2]*v[2] + pow(r+r0,2);
        //dn2 = v[2]*v[2] + pow(r-r0,2);
        if( rank == 1 )
        {
            rv[G][i] = -0.161 * alpha * exp(-pow(v[0]/0.363,2)/2);
            rv[X][i] = -8.37  * alpha * exp(-pow(v[0]/0.272,2)/2);
            rv[Y][i] =  0.93  * alpha * exp(-pow(v[0]/0.302,2)/2);
        } else if( rank == 3 ) {
            // Segment
            dz = abs(v[2]) <= len/2 ? 0 : len/2 < v[2] ? v[2]-len/2 : v[2]+len/2;
            d2 = v[0]*v[0] + v[1]*v[1] + dz*dz;
            rv[G][i] = -0.411 * alpha * exp(-d2/pow(1.14,2)/2);
            rv[X][i] = -14.60 * alpha * exp(-d2/pow(1.04,2)/2);
            rv[Y][i] =   1.70 * alpha * exp(-d2/pow(1.07,2)/2);
            /* Torus
            rv[G][i] = -0.308 * alpha * (
                exp(-dp2/pow(0.739,2)/2) +
                exp(-dn2/pow(0.739,2)/2) );
            rv[X][i] = -13.6  * alpha * (
                exp(-dp2/pow(0.634,2)/2) +
                exp(-dn2/pow(0.634,2)/2) );
            rv[Y][i] =  1.5   * alpha * (
                exp(-dp2/pow(0.665,2)/2) +
                exp(-dn2/pow(0.665,2)/2) );
            */
            /* 3 particles
            rv[G][i] = -0.411 * alpha * (
                exp(-(pow(v[0]-c1[0],2)+pow(v[1]-c1[1],2)+pow(v[2]-c1[2],2))/pow(1.14,2)/2) +
                exp(-(pow(v[0]-c2[0],2)+pow(v[1]-c2[1],2)+pow(v[2]-c2[2],2))/pow(1.14,2)/2) +
                exp(-(pow(v[0]-c3[0],2)+pow(v[1]-c3[1],2)+pow(v[2]-c3[2],2))/pow(1.14,2)/2) );
            rv[X][i] = -14.60 * alpha * (
                exp(-(pow(v[0]-c1[0],2)+pow(v[1]-c1[1],2)+pow(v[2]-c1[2],2))/pow(1.04,2)/2) +
                exp(-(pow(v[0]-c2[0],2)+pow(v[1]-c2[1],2)+pow(v[2]-c2[2],2))/pow(1.04,2)/2) +
                exp(-(pow(v[0]-c3[0],2)+pow(v[1]-c3[1],2)+pow(v[2]-c3[2],2))/pow(1.04,2)/2) );
            rv[Y][i] =   1.70 * alpha * (
                exp(-(pow(v[0]-c1[0],2)+pow(v[1]-c1[1],2)+pow(v[2]-c1[2],2))/pow(1.07,2)/2) +
                exp(-(pow(v[0]-c2[0],2)+pow(v[1]-c2[1],2)+pow(v[2]-c2[2],2))/pow(1.07,2)/2) +
                exp(-(pow(v[0]-c3[0],2)+pow(v[1]-c3[1],2)+pow(v[2]-c3[2],2))/pow(1.07,2)/2) );
            */
            /* 2 particles
            rv[G][i] = -0.411 * alpha * (
                exp(-(pow(v[0]-c1[0],2)+pow(v[1]-c1[1],2)+pow(v[2]-c1[2],2))/pow(1.14,2)/2) +
                exp(-(pow(v[0]-c2[0],2)+pow(v[1]-c2[1],2)+pow(v[2]-c2[2],2))/pow(1.14,2)/2) );
            rv[X][i] = -14.60 * alpha * (
                exp(-(pow(v[0]-c1[0],2)+pow(v[1]-c1[1],2)+pow(v[2]-c1[2],2))/pow(1.04,2)/2) +
                exp(-(pow(v[0]-c2[0],2)+pow(v[1]-c2[1],2)+pow(v[2]-c2[2],2))/pow(1.04,2)/2) );
            rv[Y][i] =   1.70 * alpha * (
                exp(-(pow(v[0]-c1[0],2)+pow(v[1]-c1[1],2)+pow(v[2]-c1[2],2))/pow(1.07,2)/2) +
                exp(-(pow(v[0]-c2[0],2)+pow(v[1]-c2[1],2)+pow(v[2]-c2[2],2))/pow(1.07,2)/2) );
            */
        }
    }
}

// Calculate reaction, and FFT-normalize by dividing by N.
void ModelG::reaction_normalize_n( int start, int end )
{
    Real pg,
         *gp = rv[G] + start,
         *xp = rv[X] + start,
         *yp = rv[Y] + start,
         *gend = rv[G] + end;
    for( ; gp != gend ; ++gp, ++xp, ++yp )
    {
        pg = (-q*(*gp) + g*(*xp)) / N;
        *yp = (b*(*xp) - u*(*yp)
               + s*(pow(*xp+X0,3)-pow(X0,3)) - (pow(*xp+X0,2)*(*yp+Y0)-pow(X0,2)*Y0)) / N;
        *xp = (p*(*gp) - *xp) / N - *yp; // shortcut
        *gp = pg;
    }
}


int main()
{
    ModelG model_g("model_g.txt");
    model_g.run();
    return 0;
}
