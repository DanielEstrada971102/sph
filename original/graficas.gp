Skip to content
Search or jump to…

Pull requests
Issues
Marketplace
Explore
 
@DanielEstrada971102 
Learn Git and GitHub without any code!
Using the Hello World guide, you’ll start a branch, write comments, and open a pull request.


DanielEstrada971102
/
sph
1
00
Code
Issues
Pull requests
Actions
Projects
Wiki
Security
Insights
Settings
sph/graficas.gp
@DanielEstrada971102
DanielEstrada971102 Se borran todos los archivos producidos por el codigo para dejarlo co…
…
Latest commit 14f6d15 2 hours ago
 History
 1 contributor
37 lines (30 sloc)  1.24 KB
  


lim = 4000
retardo = 0.00001
f = 0.01

set size square
set xrange [-5e-5/2.0:1e-3+5e-5/2.0]
set yrange [-5e-5/2.0:1e-3+5e-5/2.0]

#do for [i=0:lim]{
#  titulo = sprintf("paso = %.4d - tiempo = %.5f",i,i*5e-5)  
#  set title titulo
#  file = sprintf("./output/state_%.4d",i)
#    plot file every ::0::1599 u 2:3 w p ps 1 pt 7 lc rgb "blue" not,\
#    "" every ::1600::1680 u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
#    "" every ::1681::1759 u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
#    "" every ::1760::1840 u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
#    "" every ::1841::1919 u 2:3 w lp ps 1 pt 7 lc rgb "black" not
#    #"" every ::1920::1960 u 2:3 w lp ps 1 pt 7 lc rgb "black" not

#  pause retardo
#}

do for [i=0:lim]{
  titulo = sprintf("paso = %.4d - tiempo = %.5f",i,i*5e-5)  
  set title titulo
  file = sprintf("./output/state_%.4d",i)
    plot file every ::0::1599 u 2:3:($4*f):($5*f) w vectors lc rgb "blue" not,\
     "" every ::1600::1680 u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
     "" every ::1681::1759 u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
     "" every ::1760::1840 u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
     "" every ::1841::1919 u 2:3 w lp ps 1 pt 7 lc rgb "black" not
     #"" every ::1920::1960 u 2:3 w lp ps 1 pt 7 lc rgb "black" not

  pause retardo
}