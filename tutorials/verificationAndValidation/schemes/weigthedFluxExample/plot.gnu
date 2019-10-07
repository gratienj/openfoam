# plot temperature and heat flux:
# - with wrong and correct schemes

set title "temperature profile" font "Helvetica,12"

set xlabel "z" font "Helvetica,12"
set ylabel "T" font "Helvetica,12"

set key bottom right
set logscale y
set format y "10^{%L}"
set key font ",12"

s1=1
s2=1e6

fo(x)=s1/(s1+s2)+s2*x/(s1+s2)
fu(x)=s1/(s1+s2)+s1*x/(s1+s2)

plot (x<=0)?fu(x):fo(x) w l lc 0 t 'analytical solution', \
    "postProcessing/singleGraph/1/line_T.xy" u 1:2 w p pt 2 ps 2 lt rgb "red" t 'linear interpolation', \
    "postProcessing/singleGraph/2/line_T.xy" u 1:2 w p pt 6 ps 2 lt rgb "green" t 'harmonic interpolation'

pause -1

set key top right
set title "heat flux"

plot "postProcessing/singleGraph/1/line_flux.xy" u 1:4 w p pt 6 ps 2 lt rgb "red" t 'laplacian: linear, grad: linear', \
     "postProcessing/singleGraph/2/line_flux.xy" u 1:4 w p pt 2 ps 2 lt rgb "blue" t 'laplacian: harmonic, grad: linear',\
     "postProcessing/singleGraph/3/line_flux.xy" u 1:4 w l lt rgb "green" t 'laplacian: harmonic, grad: weightedFlux'
pause -1
quit
