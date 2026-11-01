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

static double TSTAGE; static long DSTAGE;    /* start of the current stage and its done count then (progress_start), for the ETA */
void progress_start(long done0){ TSTAGE = wtime(); DSTAGE = done0; }
static void eta(char *b, size_t n, double s){
    if(!(s >= 0) || s > 1e7){ snprintf(b, n, "?"); return; }
    long x = (long)(s + 0.5);
    if(x >= 3600) snprintf(b, n, "%ldh%02ldm", x/3600, x/60%60); else if(x >= 60) snprintf(b, n, "%ldm%02lds", x/60, x%60); else snprintf(b, n, "%lds", x);
}
/* called from inside the parallel loops; whichever thread is first after the interval reports.
 * A stage with a fixed total (the refinement's cells) also shows the time left at the stage's rate. */
void progress(const char *stage, long done, long total){
    if(PROG <= 0) return;
    double t = wtime(); if(t - TLAST < PROG) return;
    #pragma omp critical(sincqed_progress)
    if(t - TLAST >= PROG){
        TLAST = t;
        time_t now = time(NULL); char tb[32]; strftime(tb, sizeof tb, "%H:%M:%S", localtime(&now));
        char et[24] = ""; if(TSTAGE > 0 && done > DSTAGE && total > done){ char e[16]; eta(e, sizeof e, (t - TSTAGE)*(double)(total - done)/(done - DSTAGE)); snprintf(et, sizeof et, "  ETA %s", e); }
        say("  [%s  %6.0fs]  graph %d/%d %-8s %-24s %ld/%ld (%.1f%%)%s\n", tb, t - T0, CURI, CURN, CURNAME, stage,
            done, total, total > 0 ? 100.0*done/total : 0.0, et);
    }
}
/* the same line with free text after the stage (the grow, whose total is not known in advance) */
void progress_text(const char *stage, const char *text){
    if(PROG <= 0) return;
    double t = wtime(); if(t - TLAST < PROG) return;
    #pragma omp critical(sincqed_progress)
    if(t - TLAST >= PROG){
        TLAST = t;
        time_t now = time(NULL); char tb[32]; strftime(tb, sizeof tb, "%H:%M:%S", localtime(&now));
        say("  [%s  %6.0fs]  graph %d/%d %-8s %-24s %s\n", tb, t - T0, CURI, CURN, CURNAME, stage, text);
    }
}
