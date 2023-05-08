plot "percErrAnalysis.csv" every ::1 u 1:2:(sqrt(($3 - $2)**2 + ($4 - $2)**2 + ($5 - $2)**2 + ($6 - $2)**2)) with yerrorbars title "retinal beam diameter"

