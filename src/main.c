/* -------------- */
/* --- main.c --- */
/* -------------- */

#include <stdio.h>
#include <stdlib.h>


#ifdef OPENMP
#include <omp.h>
#endif

#include "def.h"
#include "nrutil.h"

#include "vdef.h"
#include "vnrutil.h"

#include "mutil.h"

#include "mandelbrot.h"
#include "pi.h"

// ------------
void info(int cores)
// ------------
{
    int p;
	 

#ifdef ENABLE_BENCHMARK
    puts("mode Benchmark ON");
    puts("DEBUG OFF");
#else
    puts("mode Benchmark OFF");
    puts("DEBUG ON");
#endif
    
    #ifdef OPENMP
    puts("OpenMP ON");
    printf("En mode %d coeurs!", cores);
	 puts("adaptez p a votre machine !");
    p = 1;
    //p = 8;
    p = 16;
    omp_set_num_threads(cores);
#endif

}
// -----------------------------
int main(int argc, char *argv[])
// -----------------------------
{
#ifdef OPENMP
	for(int cores = 8; cores < 32; cores*=2) {
    info(cores);
#else
	 info(8);
#endif
    main_mandelbrot(argc, argv);
    main_pi(argc, argv);
#ifdef OPENMP
	}
#endif

    return 0;   
}