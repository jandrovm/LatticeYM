Nt = 8
Ns = 4
dimSp = 3 #Spatial dimensions

set title sprintf("%d",Ns).'^'.sprintf("%d",dimSp).' x '.sprintf("%d",Nt).' lattice' font "Times-Roman, 14"

set xlabel 'time (MC steps)' font "Times-Roman,13"
set ylabel '<{/Symbol e}>' font "Times-Roman,14"

set key font "Times-Roman, 12" bottom right
plot 'hot_start_action.dat' w lines lw 1.2 t 'Hot Start', 'warm_start_action.dat' w lines lw 1.2 t 'Warm Start', 'cold_start_action.dat' w lines lw 1.2 t 'Cold Start'

pause -1
