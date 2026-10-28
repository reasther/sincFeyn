/* graphs.c -- the graphs linked into this executable, selectors and reference values. */
#include "sq.h"

const graph_t **ALL;
int NALL;
static int CAPALL;

void add_graph(const graph_t *g){
    for(int i = 0; i < NALL; i++) if(!strcmp(ALL[i]->name, g->name)) return;     /* first one wins */
    if(NALL == CAPALL){ CAPALL = CAPALL ? 2*CAPALL : 512; ALL = (const graph_t**)realloc(ALL, CAPALL*sizeof *ALL); }
    ALL[NALL++] = g;
}

/* selectors: a graph's name or Volkov string, or a group */
static const char *GROUPS[] = {"all2", "all3", "all4", "q3", "vp3", "vp4", "lbl4", "pi6", "lbl", NULL};
int isgroup(const char *s){ for(int i = 0; GROUPS[i]; i++) if(!strcmp(s, GROUPS[i])) return 1; return 0; }
int selected(const graph_t *g, const char *s){
    if(!strcmp(s, g->name) || !strcmp(s, g->volkov)) return 1;
    if(!strcmp(s, "all2")) return g->L == 2;
    if(!strcmp(s, "all3")) return g->L == 3 && strncmp(g->name, "w3_", 3);   /* w3_* are checks, not graphs */
    if(!strcmp(s, "all4")) return g->L == 4;
    if(!strcmp(s, "q3"))   return !strncmp(g->name, "g3_", 3) && strcmp(g->name, "g3_lbl");
    if(!strcmp(s, "vp3"))  return !strncmp(g->name, "v3_", 3);
    if(!strcmp(s, "vp4"))  return !strncmp(g->name, "v4_", 3) || !strncmp(g->name, "w4_", 3);
    if(!strcmp(s, "lbl4")) return !strncmp(g->name, "l4_", 3);
    if(!strcmp(s, "pi6"))  return !strncmp(g->name, "w4_p6", 5);
    if(!strcmp(s, "lbl"))  return !strcmp(g->name, "g3_lbl");
    return 0;
}
const graph_t *find(const char *s){
    for(int i = 0; i < NALL; i++) if(!strcmp(ALL[i]->name, s) || !strcmp(ALL[i]->volkov, s)) return ALL[i];
    return NULL;
}

/* references (values include the mirror graph) */
static const ref_t REFS[] = {
    {"3; 1-4, 2-5", -0.46764544609361365, 0, "exact 1/6+13pi^2/36-5pi^2ln2/6+5z3/4 (crossed)"},
    {"3; 1-5, 2-4",  0.77747802228274215, 0, "exact 11/48+pi^2/18 (ladder)"},
    {"2; 1-4, 3-5", -0.064020941344682816, 0, "exact -55/24+pi^2/18+pi^2ln2/3-z3/2 (corner, PSLQ)"},
    {"2; 1-5, 3-4", -0.58997802228274215, 0, "exact -1/24-pi^2/18 (SE, PSLQ)"},
    {"3; 1-6, 2-5, 4-7", -0.334695, 5e-7, "analytic, 6 digits"},
    {"VP: 1-loop vertex + Pi2", 0.015687421859103, 0, "exact 119/36-pi^2/3"},
    {"VP: 2; 1-3 +Pi2(1-3)+Pi2(1-3)+Pi2(1-3)", 0.000876865858890, 0, "1D integral (mpmath)"},
    {"LbL (6 orderings)", 0.3710052921, 1e-10, "analytic (electron loop)"},
    {NULL, 0, 0, NULL}
};
const ref_t *findref(const char *v){ for(int i = 0; REFS[i].volkov; i++) if(!strcmp(REFS[i].volkov, v)) return &REFS[i]; return NULL; }
