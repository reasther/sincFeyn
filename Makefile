# sincqed build.  No Python needed: the kernels are pre-generated C, in kernels/<class>/ (in git:
# g2 g3 v3 v4 w4 pi6) and kernels-ext/<class>/ (not in git: g4 and lbl4, from the release archive
# or generator/regen.sh; KEXT=DIR to keep them elsewhere).
#
#   make                      all classes that have kernels (g2 g3 v3 g4 v4 w4 pi6 lbl4)
#   make CLASSES="g2 g3"      only these classes (a smaller, faster-linking executable)
#   make DD=0                 __float128 (x86-64) or long double as the high-precision fallback
#                             instead of double-double (the default; needed on Apple silicon)
#   make test                 quick regression runs against known values (tests/run_tests.sh),
#                             with the executable as last built (it does not rebuild)
#   later plain `make` calls keep the CLASSES of the last build (build/classes)
#   make -j8                  parallel build; each class becomes build/libsq_<class>.a
#   make GPU=1 -j16           CUDA build bin/sincqed-gpu (needs nvcc and Python 3; see docs/GPU.md)
#
# On macOS the newest Homebrew GCC on the PATH (gcc-16 .. gcc-12) is the default, with its matching
# g++; make CC=gcc-14 picks one explicitly (the matching g++-14 follows automatically).
UNAME := $(shell uname -s)
CC  ?= gcc
CXX ?= g++
ifeq ($(origin CC),default)
ifeq ($(UNAME),Darwin)
CC := $(or $(firstword $(foreach v,16 15 14 13 12,$(shell command -v gcc-$(v) 2>/dev/null))),gcc)
else
CC := gcc
endif
endif
# Homebrew GCC targets the macOS version it was built for, while Apple's assembler and linker default
# to the running system, which gives "overriding deployment version" warnings (harmless).  Using the
# running system's version for both silences them.
ifeq ($(UNAME),Darwin)
export MACOSX_DEPLOYMENT_TARGET ?= $(shell sw_vers -productVersion 2>/dev/null | cut -d. -f1).0
endif
# The default double-double build compiles as C++.  If only CC is given (e.g. make CC=gcc-14),
# use the matching C++ compiler (g++-14), so OpenMP works on macOS without also passing CXX.
ifeq ($(origin CXX),default)
CXX := $(if $(findstring gcc,$(notdir $(CC))),$(if $(findstring /,$(CC)),$(dir $(CC)),)$(subst gcc,g++,$(notdir $(CC))),g++)
endif
OPT    ?= -O3 -march=native
OPT_K4 ?= -O2 -march=native      # four-loop kernels: MBs of straight-line code, -O2 compiles much faster at little cost
OMP    ?= -fopenmp
# The high-precision (QREAL) half of each kernel is compiled at this level via #pragma GCC optimize.
# Its runtime hardly depends on it (the fallback is <0.1% of points), but the compile time does:
# l4_vp (0.8 MB) takes 600 s as double-double C++ at -O2, 157 s with QOPT=O1, 44 s with QOPT=O0.
# QOPT= (empty) compiles the whole file at OPT/OPT_K4.
QOPT   ?= O1
DD     ?= 1
ALLCLASSES := g2 g3 v3 g4 v4 w4 pi6 lbl4 g5
# classes whose kernels are not in git: kept in $(KEXT)/<class>/ (kernels/<class>/, their place before
# 2026-10-25, is still read when $(KEXT)/<class>/ is empty)
EXTCLASSES := g4 lbl4 g5
KEXT ?= kernels-ext
kdir = $(if $(filter $(1),$(EXTCLASSES)),$(if $(wildcard $(KEXT)/$(1)/*.c),$(KEXT),$(if $(wildcard kernels/$(1)/*.c),kernels,$(KEXT))),kernels)/$(1)
$(foreach c,$(EXTCLASSES),$(if $(wildcard $(KEXT)/$(c)/*.c),,$(if $(wildcard kernels/$(c)/*.c),$(warning *** using kernels/$(c): move it to $(KEXT)/$(c) (mv kernels/$(c) $(KEXT)/)))))
# default: the classes of the last build (build/classes), else every class that has kernels
CLASSES ?= $(or $(strip $(shell cat build/classes 2>/dev/null)),$(foreach c,$(ALLCLASSES),$(if $(wildcard $(call kdir,$(c))/*.c),$(c),)))
K4CLASSES := g4 v4 w4 pi6 lbl4 g5

# The four-loop g4 kernels (release archive, or generator/regen.sh with FORM installed) are made with
# FORM.  Warn when a build uses g4 kernels from the Python path (gen.py without FORM): they run 1.5-2x slower.
ifneq ($(filter g4,$(CLASSES)),)
G4FILES := $(wildcard $(call kdir,g4)/*.c)
ifeq ($(G4FILES),)
$(warning *** CLASSES includes g4 but $(KEXT)/g4 is empty: unpack the release archive or run generator/regen.sh with FORM installed)
else
G4NONFORM := $(shell grep -L 'Dirac algebra: FORM' $(G4FILES) | wc -l | tr -d ' ')
ifneq ($(G4NONFORM),0)
$(warning *** $(G4NONFORM) of $(words $(G4FILES)) kernels in $(call kdir,g4) were not generated with FORM and run 1.5-2x slower)
endif
endif
endif

ifeq ($(DD),1)
XCC   := $(CXX)
XLANG := -x c++ -DQREAL_DD -D_Alignof=alignof   # macOS 15 SDK headers use the C11 keyword _Alignof, which g++ rejects
else
XCC   := $(CC)
XLANG := -x c
endif
LDLIBS ?= -lm
ifneq ($(DD),1)
ifneq ($(shell $(CC) -dM -E - </dev/null 2>/dev/null | grep -c __SIZEOF_FLOAT128__),0)
LDLIBS += -lquadmath
endif
endif
AR ?= ar

B := build
DRV := main io graphs lattice eval sum checkpoint
DRVOBJ := $(DRV:%=$(B)/%.o)
# keep main() in the ordinary text section, next to the driver's helpers (the arm64 branch range).
# GCC on macOS still puts main in __text_startup after all the code, so the driver objects are
# linked LAST, next to it (the kernels, reached through function pointers, can be anywhere).
MAINFLAGS ?= -fno-reorder-functions -fno-reorder-blocks-and-partition
HDRS := driver/sq.h driver/common.h driver/dd.h
# lattice points carry MAXD coordinates (driver/sq.h, default 12); five-loop graphs need D = 13, so a
# build with g5 uses 14, at 2 more bytes per stored point
MAXDFLAG := $(if $(filter g5,$(CLASSES)),-DMAXD=14)
# link order: on macOS everything main() calls directly (driver + registry) goes last, next to main
# (ld64 resolves archive members whatever the order); GNU ld needs the registry before the archives
ifeq ($(UNAME),Darwin)
LINKOBJ = $(CLASSES:%=$(B)/libsq_%.a) $(B)/registry.o $(DRVOBJ)
else
LINKOBJ = $(B)/registry.o $(CLASSES:%=$(B)/libsq_%.a) $(DRVOBJ)
endif

ifeq ($(GPU),1)
all: bin/sincqed-gpu
else
all: bin/sincqed
endif

bin/sincqed: $(B)/registry.o $(DRVOBJ) $(CLASSES:%=$(B)/libsq_%.a) | bin
	$(XCC) $(OPT) $(OMP) -o $@ $(LINKOBJ) $(LDLIBS)
	@echo "$(CLASSES)" > $(B)/classes
	@echo "built $@ with classes: $(CLASSES) (fallback: $(if $(filter 1,$(DD)),double-double,__float128/long double))"

bin $(B):
	mkdir -p $@

$(B)/main.o: driver/main.c $(HDRS) $(B)/maxd.want | $(B)
	$(XCC) $(OPT) $(OMP) $(MAINFLAGS) $(XLANG) $(MAXDFLAG) -Idriver -c -o $@ $<
$(B)/%.o: driver/%.c $(HDRS) $(B)/maxd.want | $(B)
	$(XCC) $(OPT) $(OMP) $(XLANG) $(MAXDFLAG) -Idriver -c -o $@ $<
# rewritten only when MAXDFLAG changes, so the driver is rebuilt when g5 enters or leaves CLASSES
$(B)/maxd.want: FORCE | $(B)
	@echo "$(MAXDFLAG)" | cmp -s - $@ || echo "$(MAXDFLAG)" > $@

# registry of the graphs in the chosen classes (the driver's list of graph_t descriptors)
# rewritten only when CLASSES changes, so a build with fewer classes regenerates the registry
$(B)/classes.want: FORCE | $(B)
	@echo "$(CLASSES)" | cmp -s - $@ || echo "$(CLASSES)" > $@
.PHONY: FORCE
FORCE:

# per-class listing of the kernel files (name, size, date), rewritten only when it changes: kernels
# unpacked with tar keep their archive dates, possibly older than the objects built from the files they
# replace, so the objects and the registry follow this listing rather than the kernel dates
$(B)/%.kfiles: FORCE | $(B)
	@ls $(call kdir,$*)/*.c >/dev/null 2>&1 || { echo "*** no kernels for class $* in $(call kdir,$*)/: unpack its release archive there (e.g. tar xJf kernels-$*_<date>.tar.xz) or drop it from CLASSES" >&2; exit 1; }
	@{ echo $(call kdir,$*); ls -ln $(call kdir,$*)/ 2>/dev/null | grep '\.c$$'; } | cmp -s - $@ || \
	  { echo $(call kdir,$*); ls -ln $(call kdir,$*)/ 2>/dev/null | grep '\.c$$'; } > $@

$(B)/registry.c: $(foreach c,$(CLASSES),$(wildcard $(call kdir,$(c))/*.c) $(B)/$(c).kfiles) Makefile $(B)/classes.want | $(B)
	@{ echo '/* generated by make: graphs in classes $(CLASSES) */'; echo '#include "sq.h"'; \
	  for d in $(foreach c,$(CLASSES),$(call kdir,$(c))); do grep -ho '^const graph_t G_[A-Za-z0-9_]*' $$d/*.c | sed 's/^const graph_t /extern const graph_t /; s/$$/;/'; done; \
	  echo 'void register_graphs(void){'; \
	  for d in $(foreach c,$(CLASSES),$(call kdir,$(c))); do grep -ho '^const graph_t G_[A-Za-z0-9_]*' $$d/*.c | sed 's/^const graph_t \(.*\)/    add_graph(\&\1);/'; done; \
	  echo '}'; } > $@
$(B)/registry.o: $(B)/registry.c $(HDRS)
	$(XCC) $(OPT) $(XLANG) -Idriver -c -o $@ $<

# kernels: one object per graph file, one static library per class.  In C++ a namespace-scope
# const has internal linkage, so the graph descriptor is declared extern while compiling; and the
# QREAL half is preceded by the QOPT pragma.  (awk rather than sed -i/a, which differ on macOS.)
# In the double half, each accumulation of a |.| companion sum ("; n_..._Wa[..]+=") is made
# conditional on SQ_AB, so the driver can skip that half of the work where no fallback can be needed.
# ($$$$ because KPREP is expanded twice, once by $(eval) of CLASS_RULES and once in the recipe)
KPREP = awk '/^\#define REAL QREAL/ {q = 1} /^\#define REAL QREAL/ && "$(QOPT)" != "" {print "\#pragma GCC optimize (\"$(QOPT)\")"} \
  !q && index($$$$0, "_Wa[") { o = ""; s = $$$$0; while (match(s, /; n_[A-Za-z0-9_]+_Wa(\[[^]]*\])+\+=/)) { o = o substr(s, 1, RSTART + 1) "if(SQ_AB) " substr(s, RSTART + 2, RLENGTH - 2); s = substr(s, RSTART + RLENGTH) } $$$$0 = o s } \
  {if ("$(DD)" == "1") sub(/^const graph_t G_/, "extern const graph_t G_"); print}'
define CLASS_RULES
$(B)/$(1):
	mkdir -p $$@
$(B)/$(1)/%.o: $(call kdir,$(1))/%.c $(HDRS) $(B)/$(1).kfiles | $(B)/$(1)
	$(KPREP) $$< | $(XCC) $(if $(filter $(1),$(K4CLASSES)),$(OPT_K4),$(OPT)) $(XLANG) -Idriver -c -o $$@ -
$(B)/libsq_$(1).a: $(patsubst $(call kdir,$(1))/%.c,$(B)/$(1)/%.o,$(wildcard $(call kdir,$(1))/*.c))
	rm -f $$@; $(AR) rcs $$@ $$^
endef
$(foreach c,$(ALLCLASSES),$(eval $(call CLASS_RULES,$(c))))

test:   # runs the executable as built; does not rebuild it
	@[ -x bin/sincqed ] || { echo "build first: make CLASSES=..."; exit 1; }
	sh tests/run_tests.sh bin/sincqed

clean:
	rm -rf $(B) bin

# ---- CUDA back end: make GPU=1 -> bin/sincqed-gpu (bin/sincqed stays the CPU build); see docs/GPU.md ----
# The double half of each kernel in GPUCLASSES runs on the GPU; the double-double fallback stays on the
# CPU, so the CPU kernel libraries are linked too.  build/gpu/<graph>.cu is made by scripts/gpugen.py.
# Graphs of other classes (v4 w4 pi6 lbl4: not yet tested on the GPU) run on the CPU as before.
#   make GPU=1 -j16          GPUARCH=sm_80 (A100), sm_90 (H100), ... (default: native = the GPU of this machine)
CUDA    ?= /usr/local/cuda
NVCC    ?= $(CUDA)/bin/nvcc
GPUARCH ?= native
# -fmad=false: nvcc's default fused multiply-add breaks exact cancellations in 11 g4 kernels
# (g4_191 came out at -1.46e7 instead of -0.0516); without fusing all 269 agree with the CPU.  Costs 5-10%.
NVOPT   ?= -O3 -arch=$(GPUARCH) -fmad=false
PYTHON  ?= python3
GB := $(B)/gpu
GPUCLASSES ?= $(filter g2 g3 v3 g4,$(CLASSES))
GPUNAMES := $(foreach c,$(GPUCLASSES),$(patsubst $(call kdir,$(c))/%.c,%,$(wildcard $(call kdir,$(c))/*.c)))
GPUDRV := $(DRV:%=$(GB)/drv_%.o)
# The double-double fallback (points whose terms cancel beyond double precision) stays on the CPU, but
# vectorised: scripts/ddvgen.py makes a SIMD version of each kernel's double-double half that runs 8
# points at once with AVX-512 (driver/ddv.h).  Kernels that branch on computed values (v3_p4*, w3_*)
# get an empty stub and keep the scalar fallback.  DDVGRAPHS: the graphs that get one; DDVOPT: its
# compiler flags (it is compiled like the QREAL half, cf. QOPT).
DDVGRAPHS ?= $(GPUNAMES)
DDVOPT    ?= -O1 -march=native
DB := $(B)/ddv
DDVNAMES := $(filter $(DDVGRAPHS),$(GPUNAMES))
DDVOBJ := $(DDVNAMES:%=$(DB)/%.o) $(DB)/registry.o
GPUOBJ := $(GPUNAMES:%=$(GB)/k_%.o) $(GB)/gpu.o $(GB)/registry.o
.SECONDARY: $(GPUNAMES:%=$(DB)/%.cc) $(GPUNAMES:%=$(GB)/%.cu)
$(GB) $(DB):
	mkdir -p $@
define GPU_RULES
$(GB)/%.cu: $(call kdir,$(1))/%.c scripts/gpugen.py | $(GB)
	$(PYTHON) scripts/gpugen.py $$< $$@
$(DB)/%.cc: $(call kdir,$(1))/%.c scripts/ddvgen.py scripts/gpugen.py | $(DB)
	$(PYTHON) scripts/ddvgen.py $$< $$@
endef
$(foreach c,$(GPUCLASSES),$(eval $(call GPU_RULES,$(c))))
$(DB)/%.o: $(DB)/%.cc driver/ddv.h driver/ddvk.h
	$(CXX) $(DDVOPT) -Idriver -c -o $@ $<
$(DB)/registry.c: $(B)/classes.want Makefile | $(DB)
	@{ echo '/* generated by make: SIMD double-double fallbacks (weak: a kernel that is not vectorised has none) */'; echo '#include "ddvk.h"'; \
	  for f in $(DDVNAMES); do echo "extern const sq_ddv_t GV_$$f __attribute__((weak));"; done; \
	  echo 'const sq_ddv_t *DDV_ALL[] = {'; for f in $(DDVNAMES); do echo "    &GV_$$f,"; done; echo '    0 };'; \
	  echo 'const int DDV_N = $(words $(DDVNAMES));'; } > $@
$(DB)/registry.o: $(DB)/registry.c driver/ddvk.h
	$(CC) -O2 -Idriver -c -o $@ $<
$(GB)/k_%.o: $(GB)/%.cu driver/gpu_kernel.h driver/gpu.h
	$(NVCC) $(NVOPT) -Idriver -c -o $@ $<
$(GB)/gpu.o: driver/gpu.cu driver/gpu.h $(B)/maxd.want | $(GB)
	$(NVCC) $(NVOPT) $(MAXDFLAG) -Xcompiler -fopenmp -Idriver -c -o $@ $<
$(GB)/registry.c: $(B)/classes.want Makefile | $(GB)
	@{ echo '/* generated by make: graphs with a GPU kernel */'; echo '#include "gpu.h"'; \
	  for f in $(GPUNAMES); do echo "extern const sq_gpu_t GG_$$f;"; done; \
	  echo 'const sq_gpu_t *GPU_ALL[] = {'; for f in $(GPUNAMES); do echo "    &GG_$$f,"; done; echo '    0 };'; } > $@
$(GB)/registry.o: $(GB)/registry.c driver/gpu.h
	$(CC) -O2 -Idriver -c -o $@ $<
$(GB)/drv_main.o: driver/main.c $(HDRS) driver/gpu.h driver/ddvk.h $(B)/maxd.want | $(GB)
	$(XCC) $(OPT) $(OMP) $(MAINFLAGS) $(XLANG) $(MAXDFLAG) -DSQ_GPU -DSQ_DDV -Idriver -c -o $@ $<
$(GB)/drv_%.o: driver/%.c $(HDRS) driver/gpu.h driver/ddvk.h $(B)/maxd.want | $(GB)
	$(XCC) $(OPT) $(OMP) $(XLANG) $(MAXDFLAG) -DSQ_GPU -DSQ_DDV -Idriver -c -o $@ $<
bin/sincqed-gpu: $(B)/registry.o $(GPUDRV) $(CLASSES:%=$(B)/libsq_%.a) $(GPUOBJ) $(DDVOBJ) | bin
	$(XCC) $(OPT) $(OMP) -o $@ $(B)/registry.o $(CLASSES:%=$(B)/libsq_%.a) $(GPUOBJ) $(DDVOBJ) $(GPUDRV) $(LDLIBS) -L$(CUDA)/lib64 -lcudart -Wl,-rpath,$(CUDA)/lib64
	@echo "$(CLASSES)" > $(B)/classes
	@echo "built $@ with classes: $(CLASSES); GPU kernels for: $(GPUCLASSES) ($(words $(GPUNAMES)) graphs, $(GPUARCH))"

.PHONY: all test clean
