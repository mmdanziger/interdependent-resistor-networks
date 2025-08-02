from __future__ import division,print_function
import figureManager2 as fm2
import better_fit as bf
reload(bf)
L=1500
p,S,N = fm2.calculate_critical_exponent(q=1,r=1,L=L,depType=3,min_above_pc=1.3e-3,max_above_pc=0.03,avg_pc=True)#force_pc=0.8665)
#p,S,N = fm2.calculate_critical_exponent(q=1,r=0,L=L,depType=3,min_above_pc=1.3e-4,max_above_pc=0.03,avg_pc=True)
#p,S,N = fm2.calculate_critical_exponent(q=1,r=2,L=L,depType=3,min_above_pc=0.001,max_above_pc=0.02,avg_pc=True)
f = lambda x, A, beta : A * x ** beta
fit = bf.BetterFit(p,N)
fit.set_fit_function(f)
fit.minimal_range = 0.5
fit.check_ranges()
fit.plot_errors()
#fit.filter_fits("deviation", -0.002, 0.002)
#fit.filter_fits("error",0,0.008)
#fit.filter_fits("params_err",0,0.02,getter=lambda x: x[1])
