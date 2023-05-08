plot "percErrAnalysis.csv" every ::1 u 1:($3 - $2)/$2 with linespoints, "" u 1:($4 - $2)/$2 with linespoints, "" u 1:($5 - $2)/$2 with linespoints, "" u 1:($6 - $2)/$2 with linespoints 

