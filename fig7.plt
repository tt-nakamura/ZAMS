reset
#set terminal aqua
set terminal postscript eps enhanced 24
set output "fig7.eps"
set xlabel "T_s = surface temperature  / K"
set ylabel "L = luminosity  / erg/s"
set xrange [1.e5:2.e3]
set yrange [1e30:1e39]
set label "X=0.71" at graph 0.84,0.95
set label "Y=0.27" at graph 0.84,0.88
set key bottom left Left reverse font "Helvetica,20"
set format xy "10^{%T}"
set logscale xy
plot \
	"near_stars_hr.txt" t "nearest stars" w p,\
	"fig7-9.txt" u 4:3 t "ZAMS" w l lt 1