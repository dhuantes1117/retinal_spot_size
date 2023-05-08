sterm(la, Ba, Ca) = la**2 * Ba / (la**2 - Ca)

sellmeier(l, B1, B2, B3, C1, C2, C3) = 1 + sterm(l, B1, C1) + sterm(l, B2, C2) + sterm(l, B3, C3)

mag_adj_sellmeier(l_, b1, b2, b3, c1, c2, c3) = sellmeier(l_, b1*10**(-1), b2*10**(-3), b3*10, c1*10**(-2), c2*10**(-2), c3*10**3)

B1_ =  7.516e-1
B2_ = -4.484e-3
B3_ = -1.503*10
C1_ =  1.641e-2
C2_ =  8.596e-2
C3_ = -1.028e3

set xrange [400:1400]
set samples 200

#fit mag_adj_sellmeier(x, B1_, B2_, B3_, C1_, C2_, C3_) "Vincelette-RefractiveError.csv" u ($1 * 0.001):((1 / ) - (1 / $2)) via B1_, B2_, B3_, C1_, C2_, C3_

unset key
set xlabel "wavelength (nm)"
set ylabel "n"
set title "Refractive Index of Reduced Eye Model"

set term png
set output "sellmeierEquation.png"

plot sqrt(sellmeier(x / 1000, B1_, B2_, B3_, C1_, C2_, C3_)) lt -1
