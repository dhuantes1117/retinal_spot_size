set title "retinal beam diameter with 1% variation in various parameters"
set key autotitle columnhead
set ylabel "retinal beam diameter (um)"
set xlabel "wavelength (nm)"
set xrange [400:1400]
set yrange [0:140]
plot "percErrAnalysis.csv" every ::1, "" u 1:3, "" u 1:4, "" u 1:5, "" u 1:6
