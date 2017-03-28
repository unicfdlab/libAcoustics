set terminal png enhanced font "arial,15" size 500, 500 
set output "exp-calc-kt3000-cp0.png"
#set tics front
set grid nopolar
set grid front   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000
set key right bottom
#set key bottom outside
#set title "Static power coefficient APC SF 10x7" offset -3,0
set xlabel "N, rev/s" offset -1,0
set ylabel "Cp0" offset 0,0 
set xtics 1000
set yrange [0:0.2]
plot 'apcsf_10x7_static_kt0827.txt' using 1:2 title 'APC 10x7 SF exp' with points ls 3,\
'staticCalc.dat' using 1:3 title 'APC 10x7 SF calc' with lines lw 2 lc 'red'

#set terminal win 
#replot

