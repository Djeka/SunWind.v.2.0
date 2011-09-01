 set term gif
 set autoscale x
 set autoscale y
 set output "test/res/1D/HD_N_x.gif"
 set title "t=0.0686075 [sec], Bx [Gs]"
 set xlabel "x [Rz]"
 plot 'test/res/1D/LF_HD_N_x.data' with linespoints pt 1, 'test/res/1D/Roe_HD_N_x.data' with linespoints pt 2
