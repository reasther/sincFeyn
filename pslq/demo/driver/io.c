/* io.c -- screen/log output and progress lines. */
#include "sq.h"
#include <stdarg.h>

FILE *LOGF, *RESF;
double PROG = 30, T0, TLAST;                 /* progress interval (s), start time, last report */
const char *CURNAME = "";                    /* graph being evaluated */
int CURI, CURN;                              /* its index in this run, number of graphs */

double wtime(void){
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return (double)clock()/CLOCKS_PER_SEC;
#endif
}

void say(const char *fmt, ...){
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap); fflush(stdout);
    if(LOGF){ va_start(ap, fmt); vfprintf(LOGF, fmt, ap); va_end(ap); fflush(LOGF); }
}

/* called from inside the parallel loops; whichever thread is first after the interval reports */
void progress(const char *stage, long done, long total){
    if(PROG <= 0) return;
    double t = wtime(); if(t - TLAST < PROG) return;
    #pragma omp critical(sincqed_progress)
    if(t - TLAST >= PROG){
        TLAST = t;
        time_t now = time(NULL); char tb[32]; strftime(tb, sizeof tb, "%H:%M:%S", localtime(&now));
        say("  [%s  %6.0fs]  graph %d/%d %-8s %-24s %ld/%ld (%.1f%%)\n", tb, t - T0, CURI, CURN, CURNAME, stage,
            done, total, total > 0 ? 100.0*done/total : 0.0);
    }
}
