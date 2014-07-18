// ExpTimeDiff.cpp
// Author: Matt Pulver <matt@blue-science.org>
// Date: 2014 July

#include "ExpTimeDiff.h"

bool ExpTimeDiff::save_and_exit = false;

ExpTimeDiff::ExpTimeDiff( const char* config_filename, int nreacts, int nstages ):
    nthreads(thread::hardware_concurrency()),
    config(config_filename, std::ifstream::in),
    nreacts(nreacts),nstages(nstages),
    dims(nullptr)
{
    int i, r;
    const char* wisdom_file = "wisdom-fftw.txt";
    read_config();
    cerr << "N = " << N << endl;
    Nh = getNh();
    cn = center();

    rv.resize(nreacts);
    rv_output.resize(nreacts);
    k.resize(nstages+1,vector<Complex*>(nreacts));
    planr2c.resize(nstages+1,vector<fftw_plan>(nreacts));
    planc2r.resize(nstages+1,vector<fftw_plan>(nreacts));

    fftw_init_threads();
    fftw_plan_with_nthreads(nthreads);
    fftw_import_wisdom_from_filename(wisdom_file); // returns non-zero on success
    for( i = 0 ; i <= nstages ; ++i ) for( r = 0 ; r < nreacts ; ++r )
    {
        cerr << "Planning FFTW for i = " << i << ", r = " << r << endl;
        if( i == 0 ) rv[r] = fftw_alloc_real(N);
        if( i == 0 ) rv_output[r] = fftw_alloc_real(N); // buffer for multi-threaded output
        k[i][r] = reinterpret_cast<Complex*>(fftw_alloc_complex(Nh));
        // Change FFTW_MEASURE to FFTW_PATIENT or FFTW_EXHAUSTIVE for further optimization
        planr2c[i][r] = fftw_plan_dft_r2c(rank, dims, rv[r], reinterpret_cast<fftw_complex*>(k[i][r]),
            FFTW_MEASURE | FFTW_DESTROY_INPUT );
        planc2r[i][r] = fftw_plan_dft_c2r(rank, dims, reinterpret_cast<fftw_complex*>(k[i][r]), rv[r],
            FFTW_MEASURE ); // FFTW_PRESERVE_INPUT is not available for multi-dimensional c2r
    }
    fftw_export_wisdom_to_filename(wisdom_file); // returns non-zero on success

    signal( SIGINT, sighandler ); // returns previous handler
    //signal( SIGTERM, sighandler );
}

ExpTimeDiff::~ExpTimeDiff()
{
    int i, r;
    delete[] dims;
    for( i = 0 ; i <= nstages ; ++i ) for( r = 0 ; r < nreacts ; ++r )
    {
        fftw_destroy_plan(planc2r[i][r]);
        fftw_destroy_plan(planr2c[i][r]);
        fftw_free(k[i][r]);
        if( i == 0 ) fftw_free(rv_output[r]);
        if( i == 0 ) fftw_free(rv[r]);
    }
    fftw_cleanup_threads();
}

void ExpTimeDiff::read_config()
{
    string line, name;
    stringstream ss;
    int i = 0;
    do getline( config, line ); while( line[0] == '#' );
    T = stod( line );
    do { getline( config, line ); } while( line[0] == '#' );
    dt = 1/stod( line );
    do { getline( config, line ); } while( line[0] == '#' );
    fput = stoi( line );
    do { getline( config, line ); } while( line[0] == '#' );
    ss.str(line);
    while( !ss.eof() )
    {
        L.push_back(0);
        ss >> L.back();
    }
    rank = L.size();
    if( dims != nullptr ) delete[] dims;
    dims = new int[rank];
    do { getline( config, line ); } while( line[0] == '#' );
    ss.clear();
    ss.str(line);
    N = 1;
    while( !ss.eof() )
    {
        if( rank <= i ) throw
            invalid_argument("There are more quantizations of space than the number of lengths.");
        ss >> dims[i];
        N *= dims[i++];
    }
    if( rank != i )
        throw invalid_argument("There are more lengths than the number of quantizations of space.");
}

// dims[0] * dims[1] * ... * dims[rank-2] * (dims[rank-1]/2+1)
int ExpTimeDiff::getNh() const
{
    return N/2 + N/dims[rank-1];
}

int ExpTimeDiff::center() const
{
    int i = rank-2, n = dims[rank-1];
    for( ; 0 <= i ; --i ) n = dims[i] * ( n + 1 );
    return n / 2;
}

// Given an index n in a hermitian-reduced matrix where 0 <= n < Nh,
// return the index it "would have been" in a non-reduced matrix with 0 <= n < N.
int ExpTimeDiff::hermitianMap( int n ) const
{
    div_t q = div( n, dims[rank-1]/2 + 1 );
    return dims[rank-1] * q.quot + q.rem;
}

// Return phi_n(z). phi_0(z) = exp(z).
Real ExpTimeDiff::phi_maclaurin( int n, Real z ) // static
{
    Real i = n + 8, // depends on precision of Real
         retval = 1.0;
    while( n < i ) retval = 1.0 + z*retval/i--;
    while( 1 < i ) retval /= i--;
    return retval;
}

// Return phi_n(z) using complex contour method.  phi_0(z) = exp(z).
// Only upper half is necessary due to conjugate symmetry for real z.
Real ExpTimeDiff::phi_contour( int n, Real z ) // static
{
    const int Nroots = 32; // should be a power of 2
    const Real radius = abs(z) < 2 ? 1+abs(z) : 1; // Contour must enclose z but should stay away from 0
    Real k, sum = 0;
    if( n == 0 ) return exp(z);
    for( k = 0.5/Nroots ; k < 1 ; k += Real(1)/Nroots )
        sum += phi_with_instability( n, z + polar(radius,k*M_PI) ).real();
    return sum / Nroots;
}

// Better handling of small z values
Real ExpTimeDiff::phi_calc( int n, Real z ) // static
{
    if( n == 0 ) return exp(z);
    if( abs(z) < 1.0/8 ) return phi_maclaurin( n, z );
    if( abs(z) < 8 ) return phi_contour( n, z ); // takes the longest
    return phi_with_instability( n, z );
}

void ExpTimeDiff::normalize_rv()
{
    int n;
    for( Real* &rvec : rv ) for( n = 0 ; n < N ; ++n ) rvec[n] /= N;
}

void ExpTimeDiff::sighandler( int signum ) // static
{
    if( save_and_exit ) exit(1);
    save_and_exit = true;
    cerr << "\nSignal " << signum << " caught. Program will exit at next output.\n"
         << "Another signal will exit program immediately." << endl;
}

// All double values:
// t
// Gval0 Xval0 Yval0
// Gval1 Xval1 Yval1
// ...
// Gval(N-1) Xval(N-1) Yval(N-1)
void ExpTimeDiff::dump_state()
{
    int i;
    const char *filename = "state.bin";
    boost::iostreams::filtering_ostream gzout;
    //gzout.push(boost::iostreams::gzip_compressor());
    //snprintf( filename, sizeof(filename), "state.bin", t );
    gzout.push( boost::iostreams::file_descriptor_sink(filename) );

    gzout.write( (char*)&t, sizeof(double) );
    for( i = 0 ; i < N ; ++i ) for( Real* &rvec : rv ) 
        gzout.write( (char*)&rvec[i], sizeof(Real) );
}

// Sets t and rv[][] if successful
bool ExpTimeDiff::import_state()
{
    int i;
    const char *filename = "state.bin";
    char newname[20];
    ifstream in( filename, ios::in | ios::binary );
    if( !in.good() ) return false;

    in.read( (char*)&t, sizeof(double) );
    for( i = 0 ; i < N ; ++i ) for( Real* &rvec : rv ) 
        in.read( (char*)&rvec[i], sizeof(Real) );
    if( in.good() )
    {
        in.close();
        snprintf( newname, sizeof(newname), "state%010f.bin", t );
        rename( filename, newname );
        return true;
    } else return false;
}

// t is double-precision. All other numbers are single-precision floating point in binary format.
// t
// rank L0 dim0 L1 dim1 ... L(rank-1) dim(rank-1)
// Gval0 Xval0 Yval0
// Gval1 Xval1 Yval1
// ...
// Gval(N-1) Xval(N-1) Yval(N-1)
// 
// where N = product of dims
void ExpTimeDiff::output_frame( double t )
{
    int i;
    char filename[26];
    boost::iostreams::filtering_ostream gzout;
    gzout.push(boost::iostreams::gzip_compressor());
    boost::filesystem::create_directory("frames"); // ok if already exists
    snprintf( filename, sizeof(filename), "frames/frame%010f.gz", t );
    gzout.push(boost::iostreams::file_descriptor_sink(filename));

    gzout.write( (char*)&t, sizeof(double) );
    write_float( &gzout, rank );
    for( i = 0 ; i < rank ; ++i )
    {
        write_float( &gzout, L[i] );
        write_float( &gzout, dims[i] );
    }
    for( i = 0 ; i < N ; ++i ) for( Real* &rvec : rv_output ) 
        write_float( &gzout, rvec[i] );
}

void ExpTimeDiff::write_float( boost::iostreams::filtering_ostream *pout, float f ) // static
{
    pout->write( (char*)&f, sizeof(float) );
}

void ExpTimeDiff::run()
{
    int n, r, stage, total;
    div_t q = { 0, 0 };
    clock_t next_update = clock();
    double next_frame;
    thread output_frame_thread,
           threads[nthreads];

    if( !import_state() ) set_initial_conditions();
    next_frame = floor(t*fput)/fput;

    normalize_rv();
    for( fftw_plan &plan : planr2c[0] ) fftw_execute(plan); // work in the frequency domain

    cerr << "Starting main loop. t = " << t << " to " << T << endl;
    for( ; t < T ; t += dt )
    {
        for( stage = 1 ; stage <= nstages ; ++stage ) // Each iteration is a stage
        {
            for( total = 0, n = 0 ; n < nthreads ; total += q.quot, ++n )
            {
                q = div( Nh + q.rem, nthreads );
                threads[n] = thread( &ExpTimeDiff::calc_k, this, stage, total, total+q.quot );
            }
            for( thread &th : threads ) th.join();
            for( fftw_plan &plan : planc2r[stage] ) fftw_execute(plan); // destroys k[stage]
            if( stage == 1 )
            {
                if( next_update <= clock() )
                {
                    next_update += CLOCKS_PER_SEC;
                    cerr << "t = " << t << ", rv["<<litmus<<"]["<<cn<<"] = " << rv[litmus][cn] << endl;
                    if( std::isnan( rv[litmus][cn] ) ) exit(1); // doesn't work with -ffast-math
                }
                if( next_frame <= t )
                {
                    next_frame = floor(t*fput+1)/fput; // Assumes fput is a positive integer
                    if( output_frame_thread.joinable() ) output_frame_thread.join();
                    if( save_and_exit )
                    {
                        cerr << "Saving and exiting." << endl;
                        dump_state();
                        return;
                    }
                    for( r = 0 ; r < nreacts ; ++r ) for( n = 0 ; n < N ; ++n ) rv_output[r][n] = rv[r][n];
                    output_frame_thread = thread( &ExpTimeDiff::output_frame, this, t );
                }
            }
            for( total = 0, n = 0 ; n < nthreads ; total += q.quot, ++n )
            {
                q = div( N + q.rem, nthreads );
                threads[n] = thread( &ExpTimeDiff::reaction_normalize_n, this, total, total+q.quot );
            }
            for( thread &th : threads ) th.join();
            for( fftw_plan &plan : planr2c[stage] ) fftw_execute(plan);
        }
        for( total = 0, n = 0 ; n < nthreads ; total += q.quot, ++n )
        {
            q = div( Nh + q.rem, nthreads );
            threads[n] = thread( &ExpTimeDiff::calc_k, this, 0, total, total+q.quot );
        }
        for( thread &th : threads ) th.join();
    }
    for( fftw_plan &plan : planc2r[0] ) fftw_execute(plan); // destroys k[0]
    if( output_frame_thread.joinable() ) output_frame_thread.join();
    for( r = 0 ; r < nreacts ; ++r ) for( n = 0 ; n < N ; ++n ) rv_output[r][n] = rv[r][n];
    output_frame(t);
}
