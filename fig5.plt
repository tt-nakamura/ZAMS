reset
#set terminal aqua
set terminal postscript eps enhanced 24
set output "fig5.eps"
set xlabel "T = temperature  / K"
set ylabel "dlnT / dlnP = temperature gradient"
set xrange [5577.14:1.30474e+07]
#set xrange [1.27388e+17:221191]
set yrange [0:1]
set label "M=M_{sun}" at graph 0.84,0.12
set label "X=0.71, Y=0.27" at graph 0.69,0.05
set key left Left reverse
set format x "10^{%T}"
set logscale x
plot \
	"fig5.txt" u 2:5 t '{/Symbol \321}' w l lt 1,\
	"fig5.txt" u 2:3 t '{/Symbol \321}_{ad}' w l lt 2,\
	"fig5.txt" u 2:4 t '{/Symbol \321}_{rad}' w l lt 4
