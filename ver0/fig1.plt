reset
#set terminal aqua
set terminal postscript eps enhanced 24
set output "fig1.eps"
set xlabel "T = temperature  / K"
set ylabel "{/Symbol e} = nuclear energy  / erg/g/s"
set logscale xy
set xrange [1e6:1e8]
set yrange [1e-14:1e13]
set format xy "10^{%T}"
set key left Left reverse
set label "X=0.71, Y=0.27" at graph 0.68,0.05
set label "{/Symbol r} = 80 g/cm^3" at graph 0.68,0.12
plot \
	"fig1.txt" u 1:2 t "pp chain" w l lt 1,\
	"fig1.txt" u 1:3 t "CNO cycle" w l lt 2
