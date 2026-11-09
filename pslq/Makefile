# High-precision two-loop Sinc runs behind the PSLQ closed forms.  No Python needed to build or run.
#
#   make -j8            demo/bin/sincqed (unreduced sums, double + double-double) and quad/se_qd,
#                       quad/corner_qd, quad/crossed_qd (quad-double, ~64 digits).  Linux or macOS
#                       (Apple silicon too: no __float128 is needed).  On macOS the newest Homebrew gcc/g++.
#   make test           QUICK=1 run_pslq.sh and a check against the closed forms (a few minutes)
#   sh run_pslq.sh      the runs used for the paper (README.md)
#   make quadmath       x86-64 only: the original double-quad / __float128 builds (quad/se_qq, corner_q, crossed_q)
all:
	$(MAKE) -C demo CLASSES=g2
	$(MAKE) -C quad
quadmath:
	$(MAKE) -C quad quadmath
test: all
	sh tests/check.sh
clean:
	$(MAKE) -C demo clean
	$(MAKE) -C quad clean
.PHONY: all quadmath test clean
