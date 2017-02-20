# Normalization factor for the guess
norm = 15.9/.9394

set grid

set terminal qt 1 size 1820, 980

set multiplot layout 2,1

set title  "Fit result of the S-wave decay cumulative phase in D^+ -> {/Symbol p}^+ {/Symbol p}^- {/Symbol p}^+"
set xlabel "Mass-bin low edge [GeV/c^2]"
set ylabel "Cumulative phase [deg]"

set autoscale
plot "par_guess.txt" using 1:6    with linespoints title "Guess cumulative phase", \
     "par_fit.txt"   using 1:6:7  with yerrorbars  title "Best fit parameters (w. error bars)", \
     "par_fit.txt"   using 1:6    with lines       notitle

set title  "Fit result of the S-wave decay phase difference in D^+ -> {/Symbol p}^+ {/Symbol p}^- {/Symbol p}^+"
set xlabel "Mass-bin low edge [GeV/c^2]"
set ylabel "Phase difference [deg]"

set autoscale
plot "par_guess.txt" using 1:4    with linespoints title "Guess phase difference", \
     "par_fit.txt"   using 1:4:5  with yerrorbars  title "Best fit parameters (w. error bars)", \
     "par_fit.txt"   using 1:4    with lines       notitle

unset multiplot

set terminal qt 2

set title  "Fit result of the S-wave decay intensity in D^+ -> {/Symbol p}^+ {/Symbol p}^- {/Symbol p}^+"
set xlabel "Mass-bin low edge [GeV/c^2]"
set ylabel "Decay intensity"

set autoscale
plot "par_fit.txt"   using 1:2           with steps       notitle, \
     "par_fit.txt"   using 1:2:3         with yerrorbars  title "Best fit parameters (w. error bars)", \
     "par_guess.txt" using 1:2           with linespoints title "Guess parameters",

     #"par_guess.txt" using 1:(norm * $2) with linespoints title "(15.9/.9394) * guess par"

set terminal qt 3

set title "Fit result of the S-wave decay intensity in D^+ -> {/Symbol p}^+ {/Symbol p}^- {/Symbol p}^+"
set xlabel "Re(bin amplitude)"
set ylabel "Im(bin amplitude)"

set autoscale
# I can't yet manage to unset this!
set size square
plot "par_guess.txt" using ($2 * cos(acos(-1) * $6 / 180)):($2 * sin(acos(-1) * $6 / 180)) with linespoints title "Guess parameters", \
     "par_fit.txt" using ($2 * cos(acos(-1) * $6 / 180)):($2 * sin(acos(-1) * $6 / 180)) with linespoints title "Fit parameters", \

unset size 
