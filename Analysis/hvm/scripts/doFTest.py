from ROOT import *

def do_f_test(sumsqi, sumsqip, parami, paramip, nbins):
  f = (double(sumsqi-sumsqip))/double(paramip-parami)*(double(nbins-paramip)/double(sumsqip))
  print "Value of F: ",f 
  print "p-value: ", Math.fdistribution_cdf_c(f,paramip-parami,nbins-paramip)


print "Cheby1 -> Cheby2, mumu:"
do_f_test(51550.7120715,19992.994121,2,3,42)

print "Cheby2 -> Cheby3, mumu:"
do_f_test(19992.994121,16053.254482,3,4,42)

print "Cheby3 -> Cheby4, mumu:"
do_f_test(16053.254482,15398.0883467,4,5,42)

print "Cheby4 -> Cheby5, mumu:"
do_f_test(15398.0883467,15282.1161932,5,6,42)

print "Cheby1 -> Cheby2, ee:"
do_f_test(13490.2147623,13227.8361783,2,3,42)

print "Cheby2 -> Cheby3, ee:"
do_f_test(13227.8361783,11016.2400591,3,4,42)

print "Cheby3 -> Cheby4, ee:"
do_f_test(11016.2400591,10502.6234122,4,5,42)

print "Cheby4 -> Cheby5, ee:"
do_f_test(10502.6234122,10643.7430154,5,6,42)

