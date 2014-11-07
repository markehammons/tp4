/* ----------------- */
/* --- mymacro.h --- */
/* ----------------- */


#ifndef __MY_MACRO_H__
#define __MY_MACRO_H__

// activer ou desactiver le define ci-dessous pour passer
// du mode de mise au point au mode benchmark

#define ENABLE_BENCHMARK

// -------------------------------------------
// -- ne rien ecrire en dessous de cette ligne
// -------------------------------------------

// parce que la methode classique avec des boucles ne fonctionne pas avec OpenMP
#define DUP10(X) X; X; X; X; X; X; X; X; X; X
#define DUP100(X) DUP10(X); DUP10(X); DUP10(X); DUP10(X); DUP10(X); DUP10(X); DUP10(X); DUP10(X); DUP10(X); DUP10(X)
#define DUP1000(X) DUP100(X); DUP100(X); DUP100(X); DUP100(X); DUP100(X); DUP100(X); DUP010(X); DUP100(X); DUP100(X); DUP100(X)

#define CHRONO_ITER1(X) t0 = (double) _rdtsc(); DUP10(X); t1 = (double) _rdtsc(); dt=(t1-t0)/10.0f; if(dt<tmin) tmin = dt
#define CHRONO_ITER4(X) CHRONO_ITER1(X); CHRONO_ITER1(X); CHRONO_ITER1(X); CHRONO_ITER1(X)

#ifdef ENABLE_BENCHMARK
#define CHRONO0(X,t)  tmin = 1e38; for(run=0; run<nrun; run++) { t0 = (double) _rdtsc(); for(iter=0; iter<niter; iter++) { X; } t1 = (double) _rdtsc(); dt=t1-t0; if(dt<tmin) tmin = dt; } t = tmin / (double) niter
#define CHRONO(X,t)  tmin = 1e38; CHRONO_ITER4(X); t = tmin
#define BENCH(X) X
#define VERBOSE(X)
#define DEBUG(X)
#else
#define CHRONO(X,t)  X
#define BENCH(X) 
#define VERBOSE(X) X
#define DEBUG(X) X
#endif
#endif // __MY_MACRO_H__
