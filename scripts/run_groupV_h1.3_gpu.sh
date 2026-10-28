#!/bin/sh
# Group V (all 269 g4 graphs = 518 diagrams) at h=1.3 on the A100: the second h point next to the Mac's
# h=1.5 N=57 run, so that each graph gets its own h error.
#   h=1.3, N=65 (window |z| <= 41.6, as close as N allows to h=1.5's 42), eps=1e-11 for every graph
#   (pruning loss per eps decade is larger at smaller h), predictor on, SINCQED_DDSLACK=100 (as on the Mac).
# Run from the top of a checkout of main with the FORM g4 kernels in kernels-ext/g4
# (sincqed_2026-10-10_form.tar.xz), after
#     sh scripts/check_g4.sh                         # the right 269 FORM kernels, nothing else
#     make GPU=1 -j16 CLASSES="g2 g3 g4"             # bin/sincqed-gpu (g2, g3 for the regression tests)
#     sh tests/run_tests.sh bin/sincqed-gpu
#     sh scripts/check_g4.sh bin/sincqed-gpu         # four quick g4 graphs against the Mac Pro values
# then
#     PAR=2 sh scripts/run_groupV_h1.3_gpu.sh [OUTDIR]        (default OUTDIR = groupV_h1.3)
# The kernel check runs again at the start and stops the run if kernels-ext/g4 is wrong (CHECK=0 skips it).
# PAR graphs run at the same time (separate processes sharing the GPU; 2 is what the GPU team used for g4).
# Graphs are taken longest first (by their h=1.5 times) and dealt round-robin to the PAR workers.
# Output: OUTDIR/res_<graph>.txt and log_<graph>.log per graph; at the end OUTDIR/results.txt and
# OUTDIR/sincqed.log (all graphs).  Rerunning the same command resumes: graphs with a result are skipped,
# and a graph that was killed restarts from its checkpoint (OUTDIR/ckpt<w>/).
set -u
OUT=${1:-groupV_h1.3}
BIN=${BIN:-bin/sincqed-gpu}
PAR=${PAR:-2}
H=${H:-1.3} N=${N:-65} EPS=${EPS:-1e-11}
export SINCQED_PREDICT=${SINCQED_PREDICT:-1}
export SINCQED_DDSLACK=${SINCQED_DDSLACK:-100}
GRAPHS=${GRAPHS:-"
g4_155 g4_219 g4_118 g4_229 g4_226 g4_153 g4_148 g4_213 g4_216 g4_22 g4_152 g4_119 g4_03 g4_121 g4_174
g4_107 g4_208 g4_222 g4_193 g4_137 g4_150 g4_211 g4_230 g4_151 g4_189 g4_198 g4_268 g4_170 g4_181 g4_238
g4_162 g4_166 g4_195 g4_225 g4_257 g4_168 g4_227 g4_183 g4_106 g4_256 g4_234 g4_109 g4_185 g4_176 g4_249
g4_149 g4_215 g4_196 g4_163 g4_154 g4_126 g4_199 g4_134 g4_214 g4_128 g4_74 g4_161 g4_175 g4_92 g4_96
g4_156 g4_179 g4_94 g4_143 g4_261 g4_102 g4_250 g4_251 g4_212 g4_158 g4_217 g4_244 g4_84 g4_224 g4_210
g4_264 g4_223 g4_235 g4_14 g4_160 g4_220 g4_113 g4_35 g4_101 g4_00 g4_48 g4_263 g4_10 g4_205 g4_52
g4_116 g4_31 g4_86 g4_28 g4_88 g4_182 g4_141 g4_200 g4_201 g4_254 g4_39 g4_125 g4_253 g4_172 g4_105
g4_242 g4_25 g4_54 g4_79 g4_42 g4_18 g4_13 g4_186 g4_05 g4_37 g4_190 g4_262 g4_233 g4_34 g4_76
g4_108 g4_243 g4_247 g4_77 g4_15 g4_140 g4_202 g4_02 g4_95 g4_187 g4_248 g4_122 g4_146 g4_104 g4_173
g4_110 g4_111 g4_06 g4_120 g4_98 g4_252 g4_27 g4_23 g4_100 g4_87 g4_265 g4_12 g4_36 g4_228 g4_51
g4_267 g4_04 g4_89 g4_78 g4_30 g4_266 g4_99 g4_178 g4_245 g4_49 g4_184 g4_55 g4_45 g4_133 g4_75
g4_177 g4_24 g4_192 g4_142 g4_26 g4_53 g4_204 g4_01 g4_239 g4_91 g4_80 g4_115 g4_11 g4_43 g4_20
g4_241 g4_117 g4_191 g4_33 g4_206 g4_258 g4_218 g4_16 g4_236 g4_159 g4_40 g4_171 g4_93 g4_188 g4_144
g4_138 g4_180 g4_124 g4_165 g4_231 g4_97 g4_114 g4_41 g4_127 g4_112 g4_237 g4_90 g4_164 g4_17 g4_123
g4_38 g4_167 g4_197 g4_203 g4_129 g4_19 g4_232 g4_136 g4_207 g4_259 g4_85 g4_32 g4_221 g4_194 g4_169
g4_255 g4_81 g4_07 g4_240 g4_56 g4_08 g4_57 g4_145 g4_63 g4_82 g4_47 g4_131 g4_67 g4_157 g4_44
g4_103 g4_135 g4_21 g4_46 g4_130 g4_69 g4_246 g4_209 g4_260 g4_139 g4_66 g4_29 g4_147 g4_09 g4_83
g4_68 g4_132 g4_58 g4_50 g4_60 g4_64 g4_70 g4_62 g4_59 g4_72 g4_71 g4_65 g4_61 g4_73"}
[ -x "$BIN" ] || { echo "no binary $BIN (make GPU=1 CLASSES=\"g2 g3 g4\" first, or set BIN)"; exit 1; }
if [ "${CHECK:-1}" != 0 ]; then
    c=$(sh "$(dirname "$0")/check_g4.sh") || { echo "$c" | tail -20; exit 1; }
    echo "$c" | tail -1
fi
mkdir -p "$OUT"
worker(){
    w=$1; i=0
    mkdir -p "$OUT/ckpt$w"
    for g in $GRAPHS; do
        if [ $((i % PAR)) -eq "$w" ]; then
            if ! grep -q "^$g " "$OUT/res_$g.txt" 2>/dev/null; then
                echo "$(date '+%F %T') worker $w: $g"
                SINCQED_RESULTS="$OUT/res_$g.txt" SINCQED_LOG="$OUT/log_$g.log" SINCQED_CKPT_DIR="$OUT/ckpt$w" \
                    "$BIN" run "$H" "$N" "$EPS" "$g" > "$OUT/out_$g.txt" 2>&1
            fi
        fi
        i=$((i + 1))
    done
}
w=0; while [ $w -lt "$PAR" ]; do worker $w & w=$((w + 1)); done
wait
cat "$OUT"/res_g4_*.txt | grep -v '^#' > "$OUT/results.txt"
cat "$OUT"/log_g4_*.log > "$OUT/sincqed.log"
echo "$(grep -c '^g4_' "$OUT/results.txt") of 269 graphs done; results in $OUT/results.txt"
