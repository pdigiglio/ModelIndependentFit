# Normalization factor for the guess
norm = 15.9/.9394

set grid

set title  "Fit result of the S-wave decay intensity in D^+ -> {/Symbol p}^+ {/Symbol p}^- {/Symbol p}^+"
set xlabel "Mass-bin low edge [GeV/c^2]"
set ylabel "Decay intensity"

plot "par_fit.txt"   using 1:2:3         with yerrorbars  title "Best fit parameters (w. error bars)", \
     "par_fit.txt"   using 1:2           with lines       notitle, \
     "par_guess.txt" using 1:2           with linespoints title "Guess parameters", \
     "par_guess.txt" using 1:(norm * $2) with linespoints title "(15.9/.9394) * guess par"
