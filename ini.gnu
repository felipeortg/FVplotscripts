#set out "modify.pdf"
#set term pdfcairo noenhanced color solid font "Times,14"
LW = 1.4
PS = 0.5

################################
set size 1.0, 1.0
set origin 0.0, 0.0
set multiplot
################################


################
## DATA
set size 1, 0.85
set origin 0, 0.15

unset key
set xrange[0.5:*]
#set yrange[-10.5455:-7.35996]
set logscale y 2.71828
set xzeroaxis lt 1 lc 0
set xtics nomirror offset 0,0.5
set ytics nomirror offset 0.5,0
set grid
set border 3  # lb
set style fill transparent solid 0.3 noborder
set pointintervalbox 0
set errorbars 1.4
plot "modify" index 0 using 1:(abs($2)):3 with yerr pt 4 ps PS lw LW lc rgb "#000000"
unset label
unset logscale y
unset grid

################
## CORR
set size 1, 0.2
set origin 0, 0

unset key
set xrange[0.5:*]
#set yrange[-1:1]
set style fill transparent solid 0.3 noborder
set xzeroaxis lt 1 lc 0
set xtics axis offset 0,0.6
set ytics -1,1,1 nomirror offset 0.8,0
set border 2  # l
plot \
     "modify" index 1 using 1:3 with lines             lw LW lc rgb "#C0272D",\
     "modify" index 1 using 1:3 with filledcurves y1=0       lc rgb "#C0272D",\
     "<echo '1 1'"               with points pt 7 ps PS       lc rgb "#C0272D",\
     "modify" index 1 using 1:8 with lines             lw LW lc rgb "#F96600",\
     "modify" index 1 using 1:8 with filledcurves y1=0       lc rgb "#F96600",\
     "<echo '6 1'"               with points pt 7 ps PS       lc rgb "#F96600",\
     "modify" index 1 using 1:13 with lines             lw LW lc rgb "#2F7A79",\
     "modify" index 1 using 1:13 with filledcurves y1=0       lc rgb "#2F7A79",\
     "<echo '11 1'"               with points pt 7 ps PS       lc rgb "#2F7A79",\
     "modify" index 1 using 1:18 with lines             lw LW lc rgb "#417D0A",\
     "modify" index 1 using 1:18 with filledcurves y1=0       lc rgb "#417D0A",\
     "<echo '16 1'"               with points pt 7 ps PS       lc rgb "#417D0A",\
     "modify" index 1 using 1:23 with lines             lw LW lc rgb "#800080",\
     "modify" index 1 using 1:23 with filledcurves y1=0       lc rgb "#800080",\
     "<echo '21 1'"               with points pt 7 ps PS       lc rgb "#800080",\
     "modify" index 1 using 1:28 with lines             lw LW lc rgb "#D9BB76",\
     "modify" index 1 using 1:28 with filledcurves y1=0       lc rgb "#D9BB76",\
     "<echo '26 1'"               with points pt 7 ps PS       lc rgb "#D9BB76"

################################
unset multiplot
################################