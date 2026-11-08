#!/bin/sh
# Group V: all 269 g4 graphs (518 four-loop diagrams without lepton loops), as in the paper.
#
#   sh scripts/run_groupV.sh                      h = 1.5, N = 57 (window |z| <= 42)
#   HS="1.5 1.7 1.9 2.1 2.2" sh scripts/run_groupV.sh     the Fig. 4 h scan, eps matched to h
#   python3 scripts/sum_groupV.py sincqed_results_groupV_h1.5.txt sincqed_groupV_h1.5.log
#
# Settings: predictor on (SINCQED_PREDICT=1), fallback slack 100, eps = 1e-10 for 199 graphs and
# 1e-11 for the 70 whose pruning loss at 1e-10 was 1e-3 or more, both scaled by (h/1.5)^P (P=10, the
# matched eps of Fig. 4; P=0 keeps them fixed).  N = 2 round(42/h) + 1 unless N is set.  Checkpoints
# are off (CKPT=3600 turns on hourly ones).  Run from the top of the tree after
# "bin/sincqed list | grep -c '^ *g4_'" shows 269 (scripts/check_g4.sh checks the kernel tree).
# Results and log go to sincqed_results_groupV_h<h>.txt and sincqed_groupV_h<h>.log; rerunning the same
# command skips every graph already in the results file.
# About 14 hours on a 24-thread Apple M2 Ultra at h = 1.5; 4.5, 1.6, 0.6 and 0.4 hours at h = 1.7 ... 2.2.
# keep a Mac awake for the whole run (caffeinate -i), unless already under it
if [ -z "$SQ_AWAKE" ] && command -v caffeinate >/dev/null; then export SQ_AWAKE=1; exec caffeinate -i sh "$0" "$@"; fi
HS=${HS:-${H:-1.5}} P=${P:-10} BIN=${BIN:-bin/sincqed}
export SINCQED_PREDICT=1 SINCQED_DDSLACK=${DDSLACK:-100} SINCQED_CKPT=${CKPT:-0}
for H in $HS; do
NH=${N:-$(awk -v h=$H 'BEGIN{printf "%d", 2*int(42/h + 0.5) + 1}')}
E10=$(awk -v h=$H -v p=$P 'BEGIN{printf "%.2g", 1e-10*(h/1.5)^p}'); E11=$(awk -v h=$H -v p=$P 'BEGIN{printf "%.2g", 1e-11*(h/1.5)^p}')
export SINCQED_RESULTS=${RESULTS:-sincqed_results_groupV_h$H.txt} SINCQED_LOG=${LOG:-sincqed_groupV_h$H.log}
echo "== h=$H N=$NH eps $E10 / $E11 -> $SINCQED_RESULTS"
todo(){ for g in "$@"; do grep -q "^$g " "$SINCQED_RESULTS" 2>/dev/null || printf '%s ' $g; done; }
L=$(todo g4_01 g4_02 g4_04 g4_05 g4_06 g4_07 g4_08 g4_09 g4_10 g4_11 g4_12 g4_13 g4_14 g4_15 g4_16 g4_17 g4_18 g4_19 g4_20 g4_21 g4_23 g4_24 g4_25 g4_26 g4_27 g4_28 g4_29 g4_30 g4_31 g4_32 g4_33 g4_35 g4_36 g4_37 g4_38 g4_39 g4_40 g4_41 g4_42 g4_43 g4_44 g4_45 g4_46 g4_47 g4_48 g4_49 g4_50 g4_51 g4_52 g4_53 g4_54 g4_55 g4_56 g4_57 g4_58 g4_59 g4_60 g4_61 g4_62 g4_63 g4_64 g4_65 g4_66 g4_67 g4_68 g4_69 g4_70 g4_71 g4_72 g4_73 g4_81 g4_83 g4_85 g4_86 g4_87 g4_88 g4_89 g4_90 g4_91 g4_92 g4_93 g4_96 g4_97 g4_98 g4_99 g4_100 g4_101 g4_102 g4_103 g4_104 g4_105 g4_108 g4_109 g4_110 g4_111 g4_112 g4_113 g4_114 g4_115 g4_116 g4_117 g4_120 g4_122 g4_123 g4_124 g4_125 g4_126 g4_127 g4_128 g4_129 g4_130 g4_131 g4_132 g4_133 g4_135 g4_136 g4_138 g4_139 g4_140 g4_141 g4_142 g4_143 g4_144 g4_145 g4_147 g4_157 g4_158 g4_159 g4_160 g4_161 g4_163 g4_164 g4_165 g4_167 g4_169 g4_171 g4_172 g4_173 g4_175 g4_176 g4_177 g4_178 g4_179 g4_180 g4_182 g4_184 g4_185 g4_186 g4_187 g4_188 g4_190 g4_191 g4_192 g4_194 g4_197 g4_200 g4_201 g4_202 g4_203 g4_204 g4_205 g4_206 g4_207 g4_209 g4_215 g4_217 g4_218 g4_220 g4_221 g4_224 g4_228 g4_231 g4_232 g4_233 g4_235 g4_236 g4_237 g4_239 g4_240 g4_241 g4_242 g4_243 g4_244 g4_245 g4_246 g4_247 g4_248 g4_249 g4_250 g4_251 g4_252 g4_253 g4_254 g4_255 g4_258 g4_259 g4_260 g4_263 g4_266)
[ -n "$L" ] && $BIN run $H $NH $E10 $L
L=$(todo g4_00 g4_03 g4_22 g4_34 g4_74 g4_75 g4_76 g4_77 g4_78 g4_79 g4_80 g4_82 g4_84 g4_94 g4_95 g4_106 g4_107 g4_118 g4_119 g4_121 g4_134 g4_137 g4_146 g4_148 g4_149 g4_150 g4_151 g4_152 g4_153 g4_154 g4_155 g4_156 g4_162 g4_166 g4_168 g4_170 g4_174 g4_181 g4_183 g4_189 g4_193 g4_195 g4_196 g4_198 g4_199 g4_208 g4_210 g4_211 g4_212 g4_213 g4_214 g4_216 g4_219 g4_222 g4_223 g4_225 g4_226 g4_227 g4_229 g4_230 g4_234 g4_238 g4_256 g4_257 g4_261 g4_262 g4_264 g4_265 g4_267 g4_268)
[ -n "$L" ] && $BIN run $H $NH $E11 $L
done
exit 0
