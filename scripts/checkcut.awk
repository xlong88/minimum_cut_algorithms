#!/usr/bin/awk -f
/capacity/ {cap = $5}
/^n/ {n[$2]=1}
/^a/ { if (n[$2] != n[$3]) ccap+=$4}
END {if (cap != ccap) print "MISTAKE:", cap, ccap; else print "OK:", cap, ccap}
