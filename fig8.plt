reset
#set terminal aqua
set terminal postscript eps enhanced 24
set output "fig8.eps"
set xlabel "M = mass  / g"
set ylabel "L = luminosity  / erg/s"
set label "X=0.71" at graph 0.84,0.12
set label "Y=0.27" at graph 0.84,0.05
set key left Left reverse font "Helvetica,20"
set format xy "10^{%T}"
set logscale xy
M_o = 1.98892e33
L_o = 3.839e33
set xrange [0.1*M_o:30*M_o]
set yrange [7e30:1e39]
plot \
	"popper_table2.txt" u (M_o*$2):(L_o*10**$5) t "binary stars" w p pt 1,\
	"popper_table8.txt" u (M_o*$2):(L_o*10**$4) not w p pt 1,\
	"fig7-9.txt" u 1:3 t "ZAMS" w l lt 1
