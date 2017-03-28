set encoding koi8r
#set terminal latex 
#set terminal epslatex
set terminal png enhanced font 'Times-Roman, 12' size 750, 500
#set terminal svg
set output "curle-results.png"
#set tics front
set grid nopolar
set grid front linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000
set title "Probes pressure"
set xlabel "St" offset -3,0
set ylabel "SPL, dB" offset -1,0 
#set yrange [0:0.001]
set xrange [0:0.6]
#set logscale x
#set ytics  0, 0.5
plot 'acousticData/fft-CurleAnalogy1-mic.dat' u (column(1)*0.019/68):(column(3)) title 'mic A' with lines lw 2
#'acousticData/fft-CurleAnalogy1-microphone-B.dat' u (column(1)):(column(3)) title 'mic B' with lines lw 2
