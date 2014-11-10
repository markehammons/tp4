/* ------------ */
/* --- pi.c --- */
/* ------------ */

#include <math.h>
#include <float.h>

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h> /* isdigit */
#include <string.h> /* memcpy */

#ifdef OPENMP
#include <omp.h>
#endif

#ifndef TYPE
#define TYPE static
#define BSIZE 20
#endif


#include "def.h"
#include "nrutil.h"

#include "vdef.h"
#include "vnrutil.h"

#include "mutil.h"
#include "simd_macro.h"
#include "mymacro.h"

#include <x86intrin.h> // si compilateur Intel

//#define NB_THREADS 16
#define sqr(x) x*x

double pid = 3.1415926535897932384626433832795028841971693993751058;

typedef long long int64;
/* ----------------- */
void display_math(void)
/* ----------------- */
{
    printf("FLT_MAX      = %e\n", FLT_MAX);
    printf("DBL_MAX      = %e\n", DBL_MAX);
    printf("LDBL_MAX     = %Le\n", LDBL_MAX);
    printf("FLT_EPSILON  = %.40f\n", FLT_EPSILON);
    printf("DBL_EPSILON  = %.40f\n", DBL_EPSILON);
    printf("LDBL_EPSILON = %.40Lf\n", LDBL_EPSILON);
    
    printf("max iter: %.10e\n", 1.0/sqrtf(FLT_EPSILON));
    printf("max iter: %.10e\n", 1.0/sqrt(DBL_EPSILON));
    printf("max iter: %.10Le\n", 1.0/sqrtl(LDBL_EPSILON));
}

/* ------------------ */
double integrale(int64 n)
/* ------------------ */
{
	double p = 0;
	
	//printf("n = %d", n);

	double s = 0.0d;
	double inv = 1/((double) n);	

	#pragma omp parallel for schedule(TYPE,BSIZE) reduction(+:s)		
	
	for(int64 k = 0; k < n; ++k) {
		double kn = k/((double) n);
		s += inv/(1+sqr(kn));
		
	}
	p = s * 4;
    // COMPLETER ICI
	    
	return p;
}
/* ---------------- */
double arctan1(int64 n)
/* ---------------- */
{
	 double sum = 1.0d;
	 
	 //printf("n = %d", n);
	 double addVal;
	 int chunk = n / 16;
	 
		 
#pragma omp parallel for schedule(TYPE,BSIZE) reduction(+:sum)		
	for(int k = 1; k < n+1; ++k)
	{
		double flip = (k & 1)*-2 + 1;
		//printf("s = %f",s);
		addVal = 1 + k*2;
		sum += 1/(flip*addVal);
	}	 
    
    // COMPLETER ICI
    
    return sum;
}
/* ------------------------- */
double arctan(double x, int64 n)
/* ------------------------- */
{
    double s = 0; // somme
   
    // COMPLETER ICI
    
    return s;
}
/* ----------------- */
double arctan_1(int64 n)
/* ----------------- */
{    
    // pi/4 = atan(1/1)
    
    double a1;
    double pia;
    
    a1 = arctan1(n);
    pia = 4.0 * a1;
    
    VERBOSE(printf("%10s\n", "1"));
    VERBOSE(printf("%10s  %.40f\n", "approx", a1));
    VERBOSE(printf("%10s  %.40f\n", "atan", atan(x1)));
    VERBOSE(printf("%10s  %.40f\n", "erreur", fabs(a1-atan(x1))));
    VERBOSE(puts(""));
    VERBOSE(printf("%10s  %.40f\n", "approx", pia));
    VERBOSE(printf("%10s  %.40f\n", "pi", pid));
    VERBOSE(printf("%10s  %.40f\n", "erreur", fabs(pia-pid)));
    VERBOSE(printf("%10s\n", "----------"));
    //printf("%10s  %s\n", "pi", str_pi);
    
    return pia;
}
/* ------------------- */
double arctan_2_3(int64 n)
/* ------------------- */
{
    // pi/4 = 4.atan(1/2)+atan(1/3)
    
    double a2, x2 = 1.0/2.0;
    double a3, x3 = 1.0/3.0;
    double pia;
    
    a2 = arctan(x2, n);
    a3 = arctan(x3, n);
    pia = 4.0 * (a2 + a3);
    
    VERBOSE(printf("%10s\n", "2 3"));
    VERBOSE(printf("%10s  %.40f\n", "approx", a2));
    VERBOSE(printf("%10s  %.40f\n", "atan", atan(x2)));
    VERBOSE(printf("%10s  %.40f\n", "erreur", fabs(a2-atan(x2))));
    VERBOSE(puts(""));
    VERBOSE(printf("%10s  %.40f\n", "approx", a3));
    VERBOSE(printf("%10s  %.40f\n", "atan", atan(x3)));
    VERBOSE(printf("%10s  %.40f\n", "erreur", fabs(a3-atan(x3))));
    VERBOSE(puts(""));
    VERBOSE(printf("%10s  %.40f\n", "approx", pia));
    VERBOSE(printf("%10s  %.40f\n", "pi", pid));
    VERBOSE(printf("%10s  %.40f\n", "erreur", fabs(pia-pid)));
    VERBOSE(printf("%10s\n", "----------"));
    //printf("%10s  %s\n", "pi", str_pi);
    
    return pia;
}
/* --------------------- */
double arctan_5_239(int64 n)
/* --------------------- */
{
    // pi/4 = 4.atan(1/5)+atan(1/239)
    
    double a5, x5     = 1.0/5.0;
    double a239, x239 = 1.0/239.0;
    double pia;
    // atan(1/5)   = 0.19739555985
    // atan(1/239) = 0.00418410041841
    
    a5   = arctan(x5,   n);
    a239 = arctan(x239, n);
    pia = 4.0 * (4.0 * a5 - a239);
    
    VERBOSE(printf("%10s\n", "5 239"));
    VERBOSE(printf("%10s  %.40f\n", "approx", a5));
    VERBOSE(printf("%10s  %.40f\n", "atan", atan(x5)));
    VERBOSE(printf("%10s  %.40f\n", "erreur", fabs(a5-atan(x5))));
    VERBOSE(puts(""));
    VERBOSE(printf("%10s  %.40f\n", "approx", a239));
    VERBOSE(printf("%10s  %.40f\n", "atan", atan(x239)));
    VERBOSE(printf("%10s  %.40f\n", "erreur", fabs(a239-atan(x239))));
    VERBOSE(puts(""));
    VERBOSE(printf("%10s  %.40f\n", "approx", pia));
    VERBOSE(printf("%10s  %.40f\n", "pi", pid));
    VERBOSE(printf("%10s  %.40f\n", "erreur", fabs(pia-pid)));
    VERBOSE(printf("%10s\n", "----------"));
    //printf("%10s  %s\n", "pi", str_pi);
    
    return pia;
}
/* ----------- */
void space(int n)
/* ----------- */
{
    int i;
    for(i=0; i<n; i++) {
        printf(" ");
    }
}
// -------------------------------------------------
void disp(int64 n, char *str, double pia, double dt)
// -------------------------------------------------
{
    // affichage sous la forme
    
    // n
    // nom de la fonction
    // valeur
    // erreur absolue
    // nombre de bits correct
    // nombre de cycle/iteration
    
    char *f0 = "%.0f";
    char *f1 = "%.1f";
    char *f2 = "%.2f";
    char *f3 = "%.3f";
    char *f4 = "%.4f";
    char *f5 = "%.5f";
    char *f6 = "%.6f";
    char *f7 = "%.7f";
    char *f8 = "%.8f";
    char *f9 = "%.9f";
    char *f10 = "%.10f";
    char *f11 = "%.11f";
    char *f12 = "%.12f";
    char *f13 = "%.13f";
    char *f14 = "%.14f";
    char *f15 = "%.15f";
    char *f16 = "%.16f";
    char *f17 = "%.17f";
    
    char *format;
    
    double error;
    double digit;
    int idigit;
    
    error = fabs(pia-pid);;
    digit = -log(error)/log(10);
    idigit = (int) digit;
    
    /*
     printf("e = %f\n", e);
     printf("d = %f\n", d);
     printf("n = %d\n", n);
     */
    switch(idigit) {
        case 0: format = f0; break;
        case 1: format = f1; break;
        case 2: format = f2; break;
        case 3: format = f3; break;
        case 4: format = f4; break;
        case 5: format = f5; break;
        case 6: format = f6; break;
        case 7: format = f7; break;
        case 8: format = f8; break;
        case 9: format = f9; break;
        case 10: format = f10; break;
        case 11: format = f11; break;
        case 12: format = f12; break;
        case 13: format = f13; break;
        case 14: format = f14; break;
        case 15: format = f15; break;
        case 16: format = f16; break;
        case 17: format = f17; break;
    }
    
    //printf("n = %12ld", n);
    printf("n = %10.3e", (double)n);
    printf("%10s  ", str);
    if(n==0) {
        printf(format, pia); space(21);
    } else {
        printf(format, pia); space(20-idigit);
    }
    
    printf(" %.5e %2d", error, idigit);
    printf(" %15.3f cycle/iteration\n", dt/(double)n);
}
/* ---------------------- */
void routine_arctan(int64 n)
/* ---------------------- */
{
    double a1;
    double t0, t1, dt;
    //double a2;
    //double a5;
    
    //printf("--------------------\n");
    //printf("--- ArcTangente n = %8ld ---\n", n);
    //printf("--------------------\n");
    
    t0 = (double) _rdtsc();
    a1 = arctan_1(n);
    t1 = (double) _rdtsc();
    dt = t1-t0;
    //a2 = arctan_2_3(n);
    //a5 = arctan_5_239(n);
    
    disp(n, "arctan 1",   a1, dt);
    //disp("atan 2 3", a2);
    //disp("atan 5 239", a5);
}
/* ---------------- */
void main_arctan(void)
/* ---------------- */
{
    // exemples d'appels
    
    //routine_arctan(1);
    //routine_arctan(2);
    //routine_arctan(5);
    routine_arctan(10);
    routine_arctan(100);
    routine_arctan(1000);
    routine_arctan(10000);
    routine_arctan(100000);
    routine_arctan(1000000);
    routine_arctan(10000000);
    routine_arctan(100000000);
    routine_arctan(1000000000);
#ifdef OPENMP
    routine_arctan(10000000000);
#endif
    puts("");
}
/* ------------------------- */
void routine_integrale(int64 n)
/* ------------------------- */
{
    double i;
    double t0, t1, dt;
        
    t0 = (double) _rdtsc();
    i = integrale(n);
    t1 = (double) _rdtsc();
    dt = t1-t0;
    
    disp(n, "integrale", i, dt);
}
/* ------------------- */
void main_integrale(void)
/* ------------------- */
{
    
    // exemples d'appels
    routine_integrale(10);
    routine_integrale(100);
    routine_integrale(1000);
    routine_integrale(10000);
    routine_integrale(100000);
    routine_integrale(1000000);
    routine_integrale(10000000);
    routine_integrale(100000000);
    routine_integrale(1000000000);
#ifdef OPENMP
    routine_integrale(10000000000);
#endif
    puts("");
}
// ========================================
int main_pi (int argc, const char * argv[])
// ========================================
{
    printf("sizeof(int) = %d\n",         (int) sizeof(int));
    printf("sizeof(long) = %d\n",        (int) sizeof(long));
    
    
    printf("sizeof(float) = %d\n",       (int) sizeof(float));
    printf("sizeof(double) = %d\n",      (int) sizeof(double));
    printf("sizeof(long double) = %d\n", (int) sizeof(long double));
    puts("");
    
    main_arctan(); // OK en SIMD sans OpenMP
    main_integrale();
    return 0;
}