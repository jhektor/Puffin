set terminal png
set output 'all.png'
set title 'All Postprocessors'
set xlabel 'time'
set ylabel 'values'
plot 'crysp_substep_out.dat' using 1:2 title 'e_zz' with linespoints, \
 'crysp_substep_out.dat' using 1:3 title 'fp_zz' with linespoints, \
 'crysp_substep_out.dat' using 1:4 title 'gss' with linespoints, \
 'crysp_substep_out.dat' using 1:5 title 'stress_zz' with linespoints

set output 'e_zz.png'
set ylabel 'e_zz'
plot 'crysp_substep_out.dat' using 1:2 title 'e_zz' with linespoints

set output 'fp_zz.png'
set ylabel 'fp_zz'
plot 'crysp_substep_out.dat' using 1:3 title 'fp_zz' with linespoints

set output 'gss.png'
set ylabel 'gss'
plot 'crysp_substep_out.dat' using 1:4 title 'gss' with linespoints

set output 'stress_zz.png'
set ylabel 'stress_zz'
plot 'crysp_substep_out.dat' using 1:5 title 'stress_zz' with linespoints

