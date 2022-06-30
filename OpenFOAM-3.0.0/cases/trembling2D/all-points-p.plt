set encoding koi8r
set terminal png enhanced size 900, 500 
set output "trembling_09.png"
#set tics front
set grid nopolar
set grid front linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000
set title "Probes pressure"
set xlabel "Time, s" offset -3,0
set ylabel "Pressure, Pa" offset -1,0 
set key horizontal outside bottom center 
set yrange [-3:3]
set xrange [0:0.04]
#set ytics  0, 0.5
plot 'postProcessing/probes/0/p' u (column(1)):(column(6)-101325) title 'probe A (pi/6)' with lines lw 2,\
'postProcessing/probes/0/p' u (column(1)):(column(7)-101325) title 'probe B (pi/6)' with lines lw 2,\
'postProcessing/probes/0/p' u (column(1)):(column(8)-101325) title 'probe C (pi/4)' with lines lw 2,\
'postProcessing/probes/0/p' u (column(1)):(column(9)-101325) title 'probe D (pi/4)' with lines lw 2,\
'trembling_r_0.9.csv' u 1:2 title 'eq A (0)' with lines lw 2,\
'trembling_r_0.9.csv' u 1:3 title 'eq B (pi/6)' with lines lw 2,\
'trembling_r_0.9.csv' u 1:4 title 'eq C (pi/4)' with lines lw 2,\
'trembling_r_0.9.csv' u 1:5 title 'eq D (2pi/3)' with lines lw 2

set terminal png enhanced size 750, 500 
set output "trembling_150.png"
set yrange [-0.015:0.015]

plot 'acousticData/dipoleCurleTest-time.dat' u 1:7 title 'Curle B (150)' with lines lw 2,\
'acousticData/dipoleCurleTest-time.dat' u 1:8 title 'Curle C (150)' with lines lw 2,\
'acousticData/dipoleFwhTest-time.dat' u 1:3 title 'FWH B (150)' with lines lw 2,\
'acousticData/dipoleFwhTest-time.dat' u 1:5 title 'FWH D (150)' with lines lw 2,\
'trembling_r_150.0.csv' u 1:2 title 'eq A (0)' with lines lw 2,\
'trembling_r_150.0.csv' u 1:3 title 'eq B (pi/6)' with lines lw 2,\
'trembling_r_150.0.csv' u 1:4 title 'eq C (pi/4)' with lines lw 2,\
'trembling_r_150.0.csv' u 1:5 title 'eq D (2pi/3)' with lines lw 2

#set terminal png enhanced size 750, 500 
#set output "near-field-comp.png"

#plot 'postProcessing/probes/0/p' u (column(1)):(column(6)-101325) title 'probe B (pi/6)' with lines lw 2,\
#'trembling_r_0.9.csv' u 1:3 title 'eq B (pi/6)' with lines lw 2,\
#'acousticData/dipoleCurleTest-time.dat' u 1:3 title 'Curle B (0.9)' with lines lw 2,\
#'acousticData/dipoleFwhTest-time.dat' u 1:3 title 'FWH B (0.9)' with lines lw 2
