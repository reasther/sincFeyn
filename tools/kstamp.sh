#!/bin/sh
# kstamp.sh BUILD CLASS DIR: per-kernel stamps for make.  BUILD/CLASS/NAME.ks holds the kernel's directory
# and its ls -ln line (size, date) and is rewritten only when that changes, so each kernel object depends on
# its own file alone: adding or replacing one kernel recompiles that kernel, not the class.  (Kernels
# unpacked with tar keep their archive dates, possibly older than the objects built from the files they
# replace, hence the listing rather than the kernel dates.)  BUILD/CLASS.kfiles, the whole listing, is
# rewritten when any line changes (the registry follows it).  A missing stamp for an unchanged line (the
# first build with this script) takes the old listing's date, so existing objects are not rebuilt.
b=$1/$2; mkdir -p "$b"
{ echo "$3"; ls -ln "$3"/ 2>/dev/null | grep '\.c$'; } | awk -v b="$b" -v lst="$b.kfiles" '
BEGIN { n = 0; while ((getline l < lst) > 0) { if (n++ == 0) od = l; else { k = l; sub(/\.c$/, "", k); sub(/.* /, "", k); old[k] = l } } close(lst) }
NR == 1 { dir = $0; all = $0; next }
{ all = all "\n" $0; k = $0; sub(/\.c$/, "", k); sub(/.* /, "", k); s = b "/" k ".ks"
  if (dir != od || old[k] != $0) { close(s); print dir > s; print $0 > s; close(s) }
  else if ((getline x < s) <= 0) { close(s); print dir > s; print $0 > s; close(s); system("touch -r \"" lst "\" \"" s "\"") }
  else close(s) }
END { cur = ""; m = 0; while ((getline l < lst) > 0) cur = (m++ ? cur "\n" : "") l; close(lst)
      if (cur != all) print all > lst }'
