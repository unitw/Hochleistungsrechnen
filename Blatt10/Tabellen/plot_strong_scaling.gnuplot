set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
#set ylabel "Speedup"
set ylabel "Zeit [s]"
set xlabel "Anzahl Prozesse"
set title "STRONG SCALING"

#base_time_JA=system("awk 'FNR == 2 {print $4}' STRONG_SCALING_JA.dat")
#base_time_GS=system("awk 'FNR == 2 {print $4}' STRONG_SCALING_GS.dat")

#plot  "STRONG_SCALING_JA.dat" using 1:((base_time_JA/$4)) title 'Jacobi' with linespoints,\
#      "STRONG_SCALING_GS.dat" using 1:((base_time_GS/$4)) title 'Gauss-Seidel' with linespoints 

plot  "STRONG_SCALING_JA.dat" using 1:4 title 'Jacobi' with linespoints,\
      "STRONG_SCALING_GS.dat" using 1:4 title 'Gauss-Seidel' with linespoints 