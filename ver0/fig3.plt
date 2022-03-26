reset
#set terminal aqua
set terminal postscript eps enhanced 24
set output "fig3.eps"
set xlabel "r = distance from center  / solar radius"
set ylabel "m, l, P, T  in solar units"
set xrange [0:1.105]
set yrange [0:1.1]
set label "X=0.71" at graph 0.84,0.17
set label "Y=0.27" at graph 0.84,0.1
set key left Left reverse
plot \
	"fig3.txt" u 1:2 t 'mass' w l lt 1,\
	"fig3.txt" u 1:3 t 'luminosity' w l lt 2,\
	"fig3.txt" u 1:4 t 'Pressure' w l lt 4,\
	"fig3.txt" u 1:5 t 'Temperature' w l lt 5
