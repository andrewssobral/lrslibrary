/* Stephen Becker, 11/12/08
 * Re-writing the fortran dbdqr.f
 * so that I can compile it on computers without
 * a fortran compiler
 *
 * This function calls dlartg, which is a LAPACK routine
 * that compuets plane rotations.  From the documentation,
 * "this is a slower, more accurate version of the BLAS1
 * routine DROTG" (with some syntax differenes).   */

#include <stdio.h>

#ifdef WINDOWS
extern void dlartg(double *F, double *G, double *CS, double *SN, double *R);
#else
extern void dlartg_(double *F, double *G, double *CS, double *SN, double *R);
#endif

void dbdqr(int *n, double *d, double *e, double *c1, double *c2) {
    int i;
    double cs, sn, r;

    if ( *n < 2 ) {
        return;
    }
    /* in the fortran version, the loop over i is
     * from 1 to n-1.  Fortran is 1-based, so I want
     * to loop from 0 to n-2 in C */
    for ( i=0 ; i < *n-1 ; i++ ) {
#ifdef WINDOWS
        dlartg( d + i, e + i, &cs, &sn, &r );
#else
        dlartg_( d + i, e + i, &cs, &sn, &r );
#endif
        d[i] = r;
        e[i] = sn*d[i+1];
        d[i+1] = cs*d[i+1]; /* need i < n-1 for this */
    }

    /* take care of i = n-1 case */
    i = *n-1;
#ifdef WINDOWS
    dlartg( d+i, e+i, &cs, &sn, &r );
#else
    dlartg_( d+i, e+i, &cs, &sn, &r );
#endif
    d[i] = r;
    e[i] = 0.0;
    *c1 = sn;
    *c2 = cs;
}
                
