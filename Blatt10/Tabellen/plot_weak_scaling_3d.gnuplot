set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set ztic auto    
set dgrid3d 20,20
set xlabel "Anzahl Prozesse"
set ylabel "Interlines"
set zlabel "Zeit [s]"
set title "WEAK SCALING"

#base_time_JA=system("awk 'FNR == 2 {print $4}' WEAK_SCALING_JA.dat")
#base_time_GS=system("awk 'FNR == 2 {print $4}' WEAK_SCALING_GS.dat")

splot  "WEAK_SCALING_JA.dat" using 1:3:4 title 'Jacobi' with lines,\
      "WEAK_SCALING_GS.dat" using 1:3:4 title 'Gauss-Seidel' with lines 
