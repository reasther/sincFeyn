# sincqed demo (standalone)

The unreduced one- and two-loop Sinc sums from an earlier sincFeyn driver (CPU only),
cut down to what `sincqed demo` needs: the driver and the two-loop kernels (`kernels/g2`: the one-loop
vertex g1_00, the four two-loop graphs g2_00..g2_03 and the Pi2 insertion v2_00). No Python is needed to
build or run.

    make -j8                       # bin/sincqed (on macOS the newest Homebrew gcc; make CC=gcc-14 to pick one)
    make test                      # two-loop and demo checks against 1/2 and A2 (a minute or two)
    bin/sincqed demo 1.0,0.8 1e-14 # one- and two-loop sums at two steps, eps = 1e-14
    HS="1.0,0.8,0.6,0.5" EPSS="1e-16 1e-18" THREADS=24 sh scripts/demo_scan.sh

Converged h ladder (h = 1.6 to 0.6 in steps of 0.02, eps 1e-18, lower below h = 0.76; totals in demo_conv/sums.tsv):

    THREADS=24 sh scripts/demo_converge.sh

`demo` sums over every free Schwinger parameter, s_i = e^{z_i}, z_i = k_i h (2 at one loop, 5 at two
loops; the Pi2 graph keeps both photon segments), with no series classes, external-vertex split, coarse
levels or predictor.  The sum is a flood fill from z = 0 that visits a node's neighbours while its
density |f prod s_i| exceeds eps; eps is the only truncation besides the box |z_i| <= SINCQED_ZMAX
(default 300, a safety limit: points are evaluated after a projective rescaling, so large |z| cannot overflow, and the fill stops on its own).  Results: stdout, sincqed.log and sincqed_demo.txt (SINCQED_DEMO_RESULTS); progress on
stderr every SINCQED_PROGRESS seconds (default 30).  Threads: OMP_NUM_THREADS.

Memory: the visited set is stored per column (a list of row intervals along the last variable), about
h^-4: 1.1 GB at h = 0.6, eps = 1e-16 (3.4e8 evaluations of g2_02, 27 s on 4 cores).  Time goes as h^-5.

Precision: a point is recomputed in double-double when its rounding bound exceeds G2_ABSERR (default
1e-15 per point).  For results below ~1e-12 set G2_ABSERR=1e-18 as well.

The full code (three and four loops, the generator) is the sincFeyn repository.
