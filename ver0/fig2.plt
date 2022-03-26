reset
#set terminal aqua
set terminal postscript eps enhanced 24
set output "fig2.eps"
set xlabel "T = temperature  / K"
set ylabel "{/Symbol k} = opacity  / cm^2/g"
set logscale xy
set xrange [1e2:3e7]
set yrange [2e-1:1e13]
set format xy "10^{%T}"
set label "X=0.71, Y=0.27" at graph 0.02,0.05
set label "dyne/cm^2" at graph 0.81,0.72
set key Left reverse
plot \
	"fig2.txt" u 1:2 t 'P = 10^{16}' w l lt 1,\
	"fig2.txt" u 1:3 t 'P = 10^{12}' w l lt 2,\
	"fig2.txt" u 1:4 t 'P = 10^8' w l lt 4,\
	"fig2.txt" u 1:5 t 'P = 10^4' w l lt 5
