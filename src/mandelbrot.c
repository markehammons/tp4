/* -------------------- */
/* --- mandelbrot.c --- */
/* -------------------- */

#include <math.h>

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h> /* isdigit */
#include <string.h> /* memcpy */

#ifdef OPENMP
#include <omp.h>
#endif

#include "def.h"
#include "nrutil.h"

#include "vdef.h"
#include "vnrutil.h"

#include "mutil.h"
#include "simd_macro.h"
#include "mymacro.h"

#include <x86intrin.h>  // si compilateur Intel
#define lte(x,y) _mm_or_si128(_mm_cmpeq_epi32(x,y),_mm_cmplt_epi32(x,y))

// --------------------------------------------------
int mandelbrot_scalar(float a, float b, int max_iter)
// --------------------------------------------------
{// conseil: afficher le contenu des variables dans la boucles *et* valider via excel

// COMPARER LES DEUX METHODES POUR VOIR LA PLUS EFFICACE
// WE HAVE TO COMPARE TO SEE THE MOST EFFICIENT APPROACH /* 
 /*   
    int iter = 0;
	     
    // COMPLETER ICI
	 
	 float x,y,z,yNew;
	 x = y = z = 0;
	 
	 for(iter = 0;iter < max_iter && x*x + y*y <= 4 ; ++iter) {
		 yNew = 2*x*y+b;
		 x = x*x - y*y + a;
		 y = yNew;
	 }
        
    return iter;
*/
 float aa =a; 
    float bb = b; 
    int iter = 0;
    float xn = 0.0; 
    float yn = 0.0;
    float xn1, yn1, zn;
	do{
		xn1 = xn*xn - yn*yn + aa;
		yn1 = 2*xn*yn + bb;
		zn = (xn1*xn1 + yn1*yn1);
		xn = xn1; 
		yn = yn1;
		iter++;		
	}while((iter < max_iter) && (zn<=4.0));

   	return iter;
}
// ------------------------------
void test_mandelbrot_scalar(void)
// ------------------------------
{
    int iter, max_iter = 20;
    float x0, y0;
    float x1, y1;
    
    puts("------------------------------");
    puts("--- test_mandelbrot_scalar ---");
    puts("------------------------------");
    // tests unitaire pour valider
    
    x0 = -0.8; y0 = 0.3; iter = mandelbrot_scalar(x0, y0, max_iter); printf("(%4.2f %4.2f) -> %3d\n", x0, y0, iter);
    x0 = -0.7; y0 = 0.2; iter = mandelbrot_scalar(x0, y0, max_iter); printf("(%4.2f %4.2f) -> %3d\n", x0, y0, iter);
}
// --------------------------------------------------------------
vuint32 mandelbrot_SIMD_F32(vfloat32 a, vfloat32 b, int max_iter)
// --------------------------------------------------------------
{
    // version avec test de sortie en float
    
    vuint32 iter = _mm_set1_epi32(0);
	 vfloat32 inc = _mm_set1_ps(1);
	 vfloat32 deux = _mm_set1_ps(2);
	 vfloat32 lim = _mm_set1_ps(4);
	 vuint32 ones = _mm_set1_epi32(1);
	 
	 vfloat32 x,y,yNew,z;
	 x = y = yNew = z = _mm_set1_ps(0);
	 uint16 incomplete = 0xffff;
		 
	 //if incomplete == 0, then we're no longer incrementing iter, and can halt looping early
	 for(int i = 0; i < max_iter && incomplete != 0; ++i) {
		 yNew = _mm_add_ps(_mm_mul_ps(deux, _mm_mul_ps(x,y)),b);
		 x = _mm_add_ps(_mm_sub_ps(_mm_mul_ps(x,x), _mm_mul_ps(y,y)),a);
		 
		 y = yNew;
		 
		 z = _mm_add_ps(_mm_mul_ps(x,x),_mm_mul_ps(y,y));
		 
		 //checks if values in z are <= values in lim. If a val is greater, the result in the resulting vec is 0.
		 //else, it's 0xfffff... something. the _mm_and_ps sets the increment for iter to either 1 or 0. then it converts the 
		 //vfloat32 vector into a vuint32 vec.
		 
		 vuint32 add = _mm_cvtps_epi32(_mm_and_ps(_mm_cmple_ps(z,lim), inc));
		 iter = _mm_add_epi32(add, iter);
		 
		 incomplete = _mm_movemask_epi8((vuint32)_mm_cmpeq_epi32(add, ones)); // extract results of comparison
	 }
    // COMPLETER ICI
    
    return iter;
}
// --------------------------------------------------------------
vuint32 mandelbrot_SIMD_I32(vfloat32 a, vfloat32 b, int max_iter)
// --------------------------------------------------------------
{
    // version avec test de sortie en int
    
    vuint32 iter = _mm_set1_epi32(0);
	 
	 vfloat32 shifter = _mm_set1_ps(pow(2,29));
	 vuint32 inc = _mm_set1_epi32(1);
	 vuint32 deux = _mm_set1_epi32(2);
	 
	 
	 
	 vuint32 ai,bi;
	 vuint32 lim = _mm_cvtps_epi32(_mm_mul_ps(_mm_set1_ps(4),shifter));
	 vuint32 ones = _mm_set1_epi32(1);
	 vuint32 add;
	 
	 
	 ai = _mm_cvtps_epi32(_mm_mul_ps(a,shifter));
	 bi = _mm_cvtps_epi32(_mm_mul_ps(b,shifter));
	 
	 uint16 incomplete = 0xffff;
	 vuint32 x,y,yNew,z;
	 x = y = yNew = z = _mm_set1_epi32(0);
	 
	 for(int i = 0; i < max_iter && incomplete != 0; ++i) {
		 yNew = _mm_add_epi32(_mm_mul_epi32(deux, _mm_mul_epi32(x,y)),bi);
		 x = _mm_add_epi32(_mm_sub_epi32(_mm_mul_epi32(x,x), _mm_mul_epi32(y,y)),ai);
		 
		 y = yNew;
		 
		 z = _mm_add_epi32(_mm_mul_epi32(x,x),_mm_mul_epi32(y,y));
		 
		 //checks if values in z are <= values in lim. If a val is greater, the result in the resulting vec is 0.
		 //else, it's 0xfffff... something. the _mm_and_ps sets the increment for iter to either 1 or 0. then it converts the 
		 //vfloat32 vector into a vuint32 vec.
		 add = _mm_and_si128(lte(z,lim), inc);
		 iter = _mm_add_epi32(add, iter);
		 
		 incomplete = _mm_movemask_epi8((vuint32)_mm_cmpeq_epi32(add, ones)); // extract results of comparison
	 }
	 
	 
    
    return iter;
}
// ----------------------------
void test_mandelbrot_SIMD(void)
// ----------------------------
{
    int iter = 20;
    vuint32 i;
    vfloat32 a, b;
    double t0, t1, dt;
    
    puts("----------------------------");
    puts("--- test_mandelbrot_SIMD ---");
    puts("----------------------------");
    
    a = _mm_setr_ps(-0.8, -0.7, -0.8, -0.7);
	 b = _mm_setr_ps(+0.3, +0.2, -0.3, -0.2);
        
    iter = 15;
    puts("mandelbrot_SIMD_F32");
    a = _mm_setr_ps(-1.00, -0.90, -0.80, -0.70);
    b = _mm_setr_ps(+0.40, +0.30, +0.30, +0.10);
    i = mandelbrot_SIMD_F32(a, b, iter);
    display_vuint32(i, "%10d ", "iter"); puts("");
    
    puts("mandelbrot_SIMD_I32");
    a = _mm_setr_ps(-1.00, -0.90, -0.80, -0.70);
    b = _mm_setr_ps(+0.40, +0.30, +0.30, +0.10);
    i = mandelbrot_SIMD_I32(a, b, iter);
    display_vuint32(i, "%10d ", "iter"); puts("");
}

// --------------------------------------------------------------------------------------------------------
void calc_mandelbrot_scalar(uint32 **M, int h, int w, float a0, float a1, float b0, float b1, int max_iter)
// --------------------------------------------------------------------------------------------------------
{
    // intervale de valeurs: [a0:a1]x[b0:b1]
    
    // la seule chose a modifier dans cette fonction est la ligne de pragma OpenMP
    
    int i, j;
    
    float da, db;
    float32 a, b;
    uint32 iter;
    
    da = (a1-a0)/w;
    db = (b1-b0)/h;
    
// OPENMP 
#ifdef OPENMP
#pragma omp parallel for private(a,b,iter) schedule(static, 20)// schedule(dynamique, 5)
#endif    
    
    for(i=0; i<h; i++) {
        for(j=0; j<w; j++) {
            
            // conversion (i,j) -> (x,y)
            a = a0 + i * da;
            b = b0 + j * db;
            
            iter = mandelbrot_scalar(a, b, max_iter);
            
            M[i][j] = iter;
        }
    }
}
// -----------------------------------------------------------------------------------------------------------
void calc_mandelbrot_SIMD_F32(vuint32 **M, int h, int w, float a0, float a1, float b0, float b1, int max_iter)
// -----------------------------------------------------------------------------------------------------------
{
    int i, j;
    
    float da, db;
    float sa, sb;
    vfloat32 a, b;
    vuint32 iter;
    
    da = (a1-a0)/w;
    db = (b1-b0)/h;
    

#ifdef OPENMP
	 #pragma omp parallel for schedule(static) private(a,b,iter) shared(M)	
#endif   

    for(i=0; i<h; i++) {
        for(j=0; j<w/4; j++) {
            
            // conversion (i,j) -> (x,y)
            sa = a0 + i * da;
            sb = b0 + j * db * 4;
            
            a = _mm_setr_ps(sa, sa+da, sa+2*da, sa+3*da);
            b = _mm_set1_ps(sb);
            
            iter = mandelbrot_SIMD_F32(a, b, max_iter);
            M[i][j] = iter;
        }
    }
}
// -----------------------------------------------------------------------------------------------------------
void calc_mandelbrot_SIMD_I32(vuint32 **M, int h, int w, float a0, float a1, float b0, float b1, int max_iter)
// -----------------------------------------------------------------------------------------------------------
{
    int i, j;
    
    float da, db;
    float sa, sb;
    vfloat32 a, b;
    vuint32 iter;
    
    da = (a1-a0)/w;
    db = (b1-b0)/h;

#ifdef OPENMP
	 #pragma omp parallel for schedule(static) private(a,b,iter) shared(M)	
#endif   

    for(i=0; i<h; i++) {
        for(j=0; j<w/4; j++) {
            
            // conversion (i,j) -> (x,y)
            sa = a0 + i * da;
            sb = b0 + j * db * 4;
            
            a = _mm_setr_ps(sa, sa+da, sa+2*da, sa+3*da);
            b = _mm_set1_ps(sb);
            
            iter = mandelbrot_SIMD_I32(a, b, max_iter);
            M[i][j] = iter;
        }
    }
}
// -----------------------------------------------------------------------------------
convert_ui32matrix_ui8matrix(uint32 **m32, int i0, int i1, int j0, int j1, uint8 **m8)
// -----------------------------------------------------------------------------------
{
    int i, j;
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            m8[i][j] = (uint8) m32[i][j];
        }
    }
}
// ----------------------------------------------
void bench_mandelbrot_scalar(int n, int max_iter)
// ----------------------------------------------
{
    // ne rien modifier dans cette fonction
    
    int h, w;
    int i0, i1, j0, j1;
    float a0, a1, b0, b1;
    uint32 **M32;
    uint8 **M8;
    
    
    // chronometrie
    int iter, niter = 4;
    int run, nrun = 5;
    double t0, t1, dt, tmin, t;
    double cycles;
    
    //puts("-----------------------------");
    //puts("-- bench_mandelbrot_scalar --");
    //puts("-----------------------------");
    
    h = w = n;
    M32 = ui32matrix(0, h-1, 0, w-1);
    M8  = ui8matrix(0, h-1, 0, w-1);
    
    a0 = -1.5; a1 = +0.5;
    b0 = -1.0; b1 = +1.0;
    
    CHRONO(calc_mandelbrot_scalar(M32, h, w, a0, a1, b0, b1, max_iter), cycles);  printf("scalar:   %10.2f\n", cycles/(n*n));
    
    DEBUG(convert_ui32matrix_ui8matrix(M32, 0, h-1, 0, w-1, M8));
    DEBUG(SavePGM_ui8matrix(M8, 0, h-1, 0, w-1, "M_scalar.pgm"));
    
    free_ui32matrix(M32, 0, h-1, 0, w-1);
    free_ui8matrix(M8, 0, h-1, 0, w-1);
}
// --------------------------------------------
void bench_mandelbrot_SIMD(int n, int max_iter)
// --------------------------------------------
{   
    // ne rien modifier dans cette fonction
    
    int h, w;
    float a0, a1, b0, b1;
    vuint32 **M32;
    uint32 **wM32;
    uint8 **M8;
    
    // chronometrie
    int iter, niter = 4;
    int run, nrun = 5;
    double t0, t1, dt, tmin, t;
    double cycles;
    
    //puts("---------------------------");
    //puts("-- bench_mandelbrot_SIMD --");
    //puts("---------------------------");
    h = w = n;
    
    M32 = vui32matrix(0, h-1, 0, w/4-1);
    M8  = ui8matrix(0, h-1, 0, w-1);
    wM32 = (uint32**) M32;
    
    // ne pas changer
    a0 = -1.5; a1 = +0.5;
    b0 = -1.0; b1 = +1.0;
    
    CHRONO(calc_mandelbrot_SIMD_F32(M32, h, w, a0, a1, b0, b1, max_iter), cycles);
    printf("SIMD F32: %10.2f\n", cycles/(n*n));
    
    CHRONO(calc_mandelbrot_SIMD_I32(M32, h, w, a0, a1, b0, b1, max_iter), cycles); // facultatif
    printf("SIMD I32: %10.2f\n\n", cycles/(n*n));
    
    DEBUG(convert_ui32matrix_ui8matrix(wM32, 0, h-1, 0, w-1, M8));
    DEBUG(SavePGM_ui8matrix(M8, 0, h-1, 0, w-1, "M_v.pgm"));
    
    free_vui32matrix(M32, 0, h-1, 0, w/4-1);
    free_ui8matrix(M8, 0, h-1, 0, w-1);
}

// =========================================
int main_mandelbrot(int argc, char * argv[])
// =========================================
{
    int n, max_iter; // pour avoir les meme param en scalar et SIMD ...
    
    test_mandelbrot_scalar();
    test_mandelbrot_SIMD();
    
    n = 1000; max_iter = 256;
    printf("n = %4d max_iter = %d\n", n, max_iter);
    bench_mandelbrot_scalar(n, max_iter);
    bench_mandelbrot_SIMD(n, max_iter);
    
    return 0;
}
