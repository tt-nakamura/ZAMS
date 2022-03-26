reset
#set terminal aqua
set terminal postscript eps enhanced 24
set output "fig6.eps"
set xlabel "T = temperature  / K"
set ylabel "degree of ionization"
set xrange [5577.14:1.30474e+07]
set yrange [0:1]
set label "M=M_{sun}" at graph 0.02,0.74
set label "X=0.71" at graph 0.02,0.67
set label "Y=0.27" at graph 0.02,0.6
set key left Left reverse
set format x "10^{%T}"
set logscale x
plot \
	"fig6.txt" u 2:3 t 'H^+' w l lt 1,\
	"fig6.txt" u 2:4 t 'He^+' w l lt 2,\
	"fig6.txt" u 2:5 t 'He^{++}' w l lt 4
