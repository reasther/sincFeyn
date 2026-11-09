# Extended-precision Sinc sums of the reduced two-loop graphs.  No Python needed.
#
#   make            portable quad-double builds (~64 digits; plain doubles + fma, works on Apple silicon):
#                     se_qd      self-energy g2_01, 3-D   (the PSLQ self-energy value)
#                     corner_qd  corner g2_00, 4-D        (slow; the PSLQ corner value came from ../demo)
#                     crossed_qd crossed g2_02, 4-D
#   make quadmath   the original x86-64 builds with libquadmath: se_qq (double-quad), corner_q, crossed_q (__float128)
#
# On macOS the newest Homebrew g++ is used (Apple clang has no OpenMP); make CXX=g++-14 to pick one.
UNAME := $(shell uname -s)
ifeq ($(origin CXX),default)
ifeq ($(UNAME),Darwin)
CXX := $(or $(firstword $(foreach v,16 15 14 13 12,$(shell command -v g++-$(v) 2>/dev/null))),g++)
else
CXX := g++
endif
endif
# -ffp-contract=off: the error-free transformations in qd.h must not be fused into fma by the compiler
CXXFLAGS ?= -O2
QDFLAGS  := -ffp-contract=off -DQREAL_QD -D_Alignof=alignof
OMP ?= -fopenmp
QD = se_qd corner_qd crossed_qd
QM = se_qq corner_q crossed_q
all: $(QD)
quadmath: $(QM)
se_qd: q2.cpp k_g2_01.cpp g2_01.c qd.h qio.h common.h
	$(CXX) $(CXXFLAGS) $(QDFLAGS) $(OMP) -I. q2.cpp k_g2_01.cpp -o $@
corner_qd: qn.cpp k_g2_00.cpp g2_00p.c qd.h qio.h common.h
	$(CXX) $(CXXFLAGS) $(QDFLAGS) $(OMP) -I. -DGRAPH=G_g2_00 qn.cpp k_g2_00.cpp -o $@
crossed_qd: qn.cpp k_g2_02.cpp g2_02p.c qd.h qio.h common.h
	$(CXX) $(CXXFLAGS) $(QDFLAGS) $(OMP) -I. -DGRAPH=G_g2_02 qn.cpp k_g2_02.cpp -o $@
se_qq: q2.cpp k_g2_01.cpp g2_01.c qq.h qio.h common.h
	$(CXX) $(CXXFLAGS) $(OMP) -DQREAL_QQ -I. q2.cpp k_g2_01.cpp -lquadmath -o $@
corner_q: qn.cpp k_g2_00.cpp g2_00p.c qio.h common.h
	$(CXX) $(CXXFLAGS) $(OMP) -I. -DGRAPH=G_g2_00 qn.cpp k_g2_00.cpp -lquadmath -o $@
crossed_q: qn.cpp k_g2_02.cpp g2_02p.c qio.h common.h
	$(CXX) $(CXXFLAGS) $(OMP) -I. -DGRAPH=G_g2_02 qn.cpp k_g2_02.cpp -lquadmath -o $@
clean:
	rm -f $(QD) $(QM)
.PHONY: all quadmath clean
