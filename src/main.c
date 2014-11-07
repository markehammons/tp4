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
void info(void)
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
    puts("adaptez p a votre machine !");
    p = 1;
    //p = 8;
    p = 16;
    omp_set_num_threads(p);
#endif

}
// -----------------------------
int main(int argc, char *argv[])
// -----------------------------
{
    info();
    main_mandelbrot(argc, argv);
    main_pi(argc, argv);
    return 0;   
}