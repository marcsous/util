/*
 * Page-wise matrix SVD decomposition based on LAPACK and OMP threads.
 * 
 * Syntax:
 *          S = pagesvd(X)
 *          [U S V] = pagesvd(X)
 *          [__] = pagesvd(X,'econ')
 *          [__] = pagesvd(__,outputForm)
 *
 *          outputForm can be 'vector' or 'matrix' [or 'trans' to return V' not V]
 *
 * Compile:
 *          (MATLAB) mex pagesvd.cpp -lmwlapack -R2018a
 *          (Octave) mex pagesvd.cpp -lopenblas -llapacke -R2018a
 */
#include "mex.h"
#include <omp.h>
#include <cmath>
#include <string>
#include <algorithm>

/* lapacke.h allows a custom complex type */
#define lapack_complex_float  mxComplexSingle
#define lapack_complex_double mxComplexDouble
#include "lapacke.h"

/* only works with interleaved complex */
#if !MX_HAS_INTERLEAVED_COMPLEX
#error "This MEX-file must be compiled with the -R2018a flag."
#endif

template <typename T>
inline void transpose(T *A, int m, int n);

template <typename T>
inline void shiftrows(T *S, int m, int n);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* check arguments */
    if (nrhs > 3) mexErrMsgTxt("Too many input arguments.");    
    if (nlhs > 3) mexErrMsgTxt("Too many output arguments.");          
    if (nrhs < 1) mexErrMsgTxt("Not enough input arguments.");
    if (!mxIsDouble(prhs[0]) && !mxIsSingle(prhs[0])) mexErrMsgTxt("First argument must be numeric array.");
    
    char jobz       = (nlhs > 1) ? 'A' : 'N'; // A(ll), S(mall) or N(o U or V)
    char outputForm = (nlhs > 1) ? 'M' : 'V'; // M(atrix), V(ector) or T(rans)
    bool trans = false; // return transpose V' instead of V [outputForm = 'T']

    for (int i = 1; i < nrhs; i++)
    {   
        char *c = mxArrayToString(prhs[i]); // no need to mxFree
        std::string str = c ? c : "";
        
        if (str.compare("econ")==0)
            jobz = (nlhs > 1) ? 'S' : 'N';
        else if (str.compare("vector")==0)
            outputForm = 'V';
        else if (str.compare("matrix")==0)
            outputForm = 'M';      
        else if (str.compare("trans")==0)
            trans = outputForm = 'M';
        else
            mexErrMsgTxt("Arguments can be 'econ', 'vector' or 'matrix'.");
    }

    /* input array: a */  
    const mxArray *a = prhs[0];
    const mwSize *adims = mxGetDimensions(a);    
    mwSize ndim = mxGetNumberOfDimensions(a);
    
    lapack_int m = ndim>0 ? adims[0] : 1;
    lapack_int n = ndim>1 ? adims[1] : 1;
    lapack_int p = ndim>2 ? adims[2] : 1;
    for (int i = 3; i < ndim; i++) p *= adims[i]; // stack all higher dimensions
    
    lapack_int mn = std::min(m, n);
    lapack_int mx = std::max(m, n);
    
    /* output arrays: s, u, v */   
    lapack_int sdims[3] = {m,n,p}; 
    if (jobz=='S' || jobz=='N') sdims[0] = sdims[1] = mn;
    if (outputForm=='V') sdims[1] = 1;

    lapack_int udims[3] = {m,m,p};
    if (jobz=='S') udims[1] = mn;
    if (jobz=='N') udims[0] = udims[1] = 1;

    lapack_int vdims[3] = {n,n,p};
    if (jobz=='S') vdims[0] = mn;
    if (jobz=='N') vdims[0] = vdims[1] = 1;

#if !MATLAB_MEX_FILE
/* workaround for octave mxCreateNumericArray bug #64687 */
#include "version.h"
if (std::stod(OCTAVE_VERSION)<=8.3 && mxIsComplex(a)) udims[2] = vdims[2] = 2*p;
#endif
    
    mwSize xdims[std::max<mwSize>(ndim,3)]; // for casting lapack_int[] into mwSize[]
    #define cast(arr) std::copy_n(arr,3,xdims)-3 // returns pointer to start of xdims
    
    mxArray *s = mxCreateNumericArray(3, cast(sdims), mxGetClassID(a), mxREAL);    
    mxArray *u = mxCreateNumericArray(3, cast(udims), mxGetClassID(a), mxIsComplex(a) ? mxCOMPLEX : mxREAL);    
    mxArray *v = mxCreateNumericArray(3, cast(vdims), mxGetClassID(a), mxIsComplex(a) ? mxCOMPLEX : mxREAL);
    if (!s || !u || !v) mexErrMsgTxt("Insufficent memory (s, u, v, info).");
   
    /* get no. threads from matlab (maxNumCompThreads) */
    mxArray *mex_plhs[1], *mex_prhs[0];
    mexCallMATLAB(1, mex_plhs, 0, mex_prhs, "maxNumCompThreads");
    int nthreads = (int)mxGetScalar(mex_plhs[0]);
    if (nthreads==1) mexWarnMsgTxt("pagesvd threads equals 1. Try increasing maxNumCompThreads.");

    /* catch errors in omp block (mexErrMsgTxt crashes). on error, info stores bad matrix index */
    volatile int info = 0;

    
/* run in parallel on single threads */
#pragma omp parallel num_threads(nthreads)
if (m*n*p)
{ 
    /* copy of i-th matrix of a (lapack overwrites) */
    void *a_i = LAPACKE_malloc( m * n * mxGetElementSize(a) ); 
    if (!a_i) info = -1; // cannot exit out of the omp block
    
    /* svd and transpose v */   
    #pragma omp for schedule(static,1)
    for (int i = 0; i < p; i++)
    {
        /* point to the i-th matrix (use char* for bytes) */
        void *s_i = (char*)mxGetData(s) + i * sdims[0] * sdims[1] * mxGetElementSize(s);
        void *u_i = (char*)mxGetData(u) + i * udims[0] * udims[1] * mxGetElementSize(u);
        void *v_i = (char*)mxGetData(v) + i * vdims[0] * vdims[1] * mxGetElementSize(v);

        /* on error skip */
        if (info != 0)
            ;
        /* real float */
        else if(!mxIsComplex(a) && !mxIsDouble(a))
        {
            std::copy_n((float*)mxGetData(a) + i * m * n, m * n, (float*)a_i);
            if (LAPACKE_sgesdd(LAPACK_COL_MAJOR, jobz, m, n, (float*)a_i, adims[0], (float*)s_i, (float*)u_i, udims[0], (float*)v_i, vdims[0])) info = i+1;
            if (trans==false) transpose((float*)v_i, vdims[0], vdims[1]);
            if (outputForm=='M') shiftrows((float*)s_i, sdims[0], sdims[1]);
        }
        /* real double */
        else if(!mxIsComplex(a) &&  mxIsDouble(a))
        {
            std::copy_n((double*)mxGetData(a) + i * m * n, m * n, (double*)a_i);   
            if (LAPACKE_dgesdd(LAPACK_COL_MAJOR, jobz, m, n, (double*)a_i, adims[0], (double*)s_i, (double*)u_i, udims[0], (double*)v_i, vdims[0])) info = i+1;
            if (trans==false) transpose((double*)v_i, vdims[0], vdims[1]);
            if (outputForm=='M') shiftrows((double*)s_i, sdims[0], sdims[1]);
        }   
        /* complex float */
        else if( mxIsComplex(a) && !mxIsDouble(a))
        {
            std::copy_n((mxComplexSingle*)mxGetData(a) + i * m * n, m * n, (mxComplexSingle*)a_i);   
            if (LAPACKE_cgesdd(LAPACK_COL_MAJOR, jobz, m, n, (mxComplexSingle*)a_i, adims[0], (float*)s_i, (mxComplexSingle*)u_i, udims[0], (mxComplexSingle*)v_i, vdims[0])) info = i+1;
            if (trans==false) transpose((mxComplexSingle*)v_i, vdims[0], vdims[1]);
            if (outputForm=='M') shiftrows((float*)s_i, sdims[0], sdims[1]);      
        }
        /* complex double */
        else if( mxIsComplex(a) &&  mxIsDouble(a))
        {
            std::copy_n((mxComplexDouble*)mxGetData(a) + i * m * n, m * n, (mxComplexDouble*)a_i);   
            if (LAPACKE_zgesdd(LAPACK_COL_MAJOR, jobz, m, n, (mxComplexDouble*)a_i, adims[0], (double*)s_i, (mxComplexDouble*)u_i, udims[0], (mxComplexDouble*)v_i, vdims[0])) info = i+1;
            if (trans==false) transpose((mxComplexDouble*)v_i, vdims[0], vdims[1]);
            if (outputForm=='M') shiftrows((double*)s_i, sdims[0], sdims[1]);
        }
      
        /* lapack doesn't catch Inf but returns NaN singular values */
        if (mxIsDouble(a) ? std::isnan(*(double*)s_i) : std::isnan(*(float*)s_i)) info = i+1; 
        
    } /* end of pragma omp for loop */

    if (a_i) LAPACKE_free(a_i);
    
} /* end of pragma omp parallel block */


    /* now can throw error */
    if (info)
    {
        std::string str = "LAPACKE_ gesdd() failed ";
        if(!mxIsComplex(a) && !mxIsDouble(a)) str.replace(8,1,"s");
        if(!mxIsComplex(a) &&  mxIsDouble(a)) str.replace(8,1,"d");
        if( mxIsComplex(a) && !mxIsDouble(a)) str.replace(8,1,"c");
        if( mxIsComplex(a) &&  mxIsDouble(a)) str.replace(8,1,"z");
        if (info == -1) str += "due to insufficient memory (a_i).";
        else str += "at matrix " + std::to_string(info) + ".";
        mexErrMsgTxt(str.c_str());
    }
    
    /* reshape to match input */
    if (trans==false) std::swap(vdims[0], vdims[1]);
    
    std::copy_n(adims, ndim, xdims);
    xdims[0] = sdims[0]; xdims[1] = sdims[1]; mxSetDimensions(s, xdims, ndim);    
    xdims[0] = udims[0]; xdims[1] = udims[1]; mxSetDimensions(u, xdims, ndim);
    xdims[0] = vdims[0]; xdims[1] = vdims[1]; mxSetDimensions(v, xdims, ndim);

    if (nlhs > 2)
    {
        plhs[0] = u;
        plhs[1] = s;
        plhs[2] = v;
    }
    else if (nlhs > 1)
    {
        plhs[0] = u;
        plhs[1] = s;
        mxDestroyArray(v);
    }
    else if (nlhs >= 0)
    {
        plhs[0] = s;
        mxDestroyArray(u);
        mxDestroyArray(v);
    }
}
      

// In-place matrix shift
template <typename T>
inline void shiftrows(T *S, int m, int n)
{
    for (int i = 1; i < std::min(m, n); i++)
        std::swap(S[i], S[i+m*i]);
}


// In-place complex conjugate
template<typename T>
inline void conjugate(T *A, int m, int n)
{
    for (int i = 0; i < m*n; i++)
        A[i].imag = -A[i].imag;
}
template<> void conjugate(float *, int, int) {}
template<> void conjugate(double*, int, int) {}


// In-place matrix transpose
template<typename T>
inline void transpose(T *A, int m, int n)
{
    conjugate(A, m, n);
    
    if (m < n) std::swap(m, n);

    while (m > 1 && n > 1)
    {
        for (int i = 1; i < m; i++)
            std::rotate(A+i, A+i*n, A+i*n+1);
        
        A += m;
        n -= 1;
    }
}
