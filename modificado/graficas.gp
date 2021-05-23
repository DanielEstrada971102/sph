lim = 4000
retardo = 0.0001
f = .6

nPart = 30

parts = nPart*nPart - 1
bBorder = (parts + 1) + 2*nPart
rBorder = (bBorder + 1) + 2*nPart -2
tBorder = (rBorder + 1) + 2*nPart
lBorder = (tBorder + 1) + 2*nPart -2
Lx = 1
Ly = 1
 
set size square
set xrange [-5e-5/2.0:Lx+5e-5/2.0]
set yrange [-5e-5/2.0:Ly+5e-5/2.0]


# do for [i=0:lim]{
#   titulo = sprintf("paso = %.4d - tiempo = %.5f",i,i*5e-5)  
#   set title titulo
#   file = sprintf("./output/state_%.4d",i)
#     plot file every ::0::parts u 2:3 w p ps 1 pt 7 lc rgb "blue" not,\
# 	 "" every ::parts+1::bBorder u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
# 	 "" every ::bBorder+1::rBorder u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
# 	 "" every ::rBorder+1::tBorder u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
# 	 "" every ::tBorder+1::lBorder u 2:3 w lp ps 1 pt 7 lc rgb "black" not
# 	 #"" every ::1920::1960 u 2:3 w lp ps 1 pt 7 lc rgb "black" not

#   pause retardo
# }

do for [i=0:lim]{
  titulo = sprintf("paso = %.4d - tiempo = %.5f \t #particulas=%d", i, i*5e-5, nPart*nPart)  
  set title titulo
  file = sprintf("./output/state_%.4d",i)
    plot file every ::0::parts u 2:3:($4*f):($5*f) w vectors lc rgb "blue" not,\
   "" every ::parts+1::bBorder u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
     "" every ::bBorder+1::rBorder u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
     "" every ::rBorder+1::tBorder u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
     "" every ::tBorder+1::lBorder u 2:3 w lp ps 1 pt 7 lc rgb "black" not
	 #"" every ::1920::1960 u 2:3 w lp ps 1 pt 7 lc rgb "black" not

  pause retardo
}

