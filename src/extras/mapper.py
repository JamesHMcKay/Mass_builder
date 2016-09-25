#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special


xp = "o";yp = "o";zp = "o";up = "o";vp = "o"


# map from V(x,y,z,u) to -TVI(u,x,z,y)

x = "a"
y = "b"
z = "c"
u = "d"
v = "e"

print "for the FeynCalc function V(a, b, c ,d ) "


# now define TSIL primed parameters

ap = y
bp = u
cp = z
dp = x
ep = "o"

# so the corresponding TSIL function we seek is
print "the corresponding TSIL function is:"
print "TVI ( " + ap + " , " + bp +" , " + cp +" , " + dp +" )"

xp = "o";yp = "o";zp = "o";up = "o";vp = "o"
# using get function Uzxyv

zp = ap
xp = bp
yp = cp
vp = dp

#print " "
print "add_eval_obj( " + xp + " , " + yp + " , " + zp + " , " + up + " , " + vp + " ) "
#print "get_function ( " + "Uzxyv )"
#print " "


xp = "o";yp = "o";zp = "o";up = "o";vp = "o"
# using get function Uuyxv

up = ap
yp = bp
xp = cp
vp = dp

#print " "
print "add_eval_obj( " + xp + " , " + yp + " , " + zp + " , " + up + " , " + vp + " ) "
#print "get_function ( " + "Uuyxv )"
#print " "


xp = "o";yp = "o";zp = "o";up = "o";vp = "o"
# using get function Uxzuv

xp = ap
zp = bp
up = cp
vp = dp

#print " "
print "add_eval_obj( " + xp + " , " + yp + " , " + zp + " , " + up + " , " + vp + " ) "
#print "get_function ( " + "Uxzuv )"
#print " "


xp = "o";yp = "o";zp = "o";up = "o";vp = "o"
# using get function Uyuzv

yp = ap
up = bp
zp = cp
vp = dp

#print " "
print "add_eval_obj( " + xp + " , " + yp + " , " + zp + " , " + up + " , " + vp + " ) "
#print "get_function ( " + "Uyuzv )"
#print " "

#############################################################################
#####  now do the same as above with the last two arguments switched
###########################################################################

xp = "o";yp = "o";zp = "o";up = "o";vp = "o"
# using get function Uzxvy

zp = ap
xp = bp
vp = cp
yp = dp

#print " "
print "add_eval_obj( " + xp + " , " + yp + " , " + zp + " , " + up + " , " + vp + " ) "
#print "get_function ( " + "Uzxvy )"
#print " "


xp = "o";yp = "o";zp = "o";up = "o";vp = "o"
# using get function Uuyvx

up = ap
yp = bp
vp = cp
xp = dp

#print " "
print "add_eval_obj( " + xp + " , " + yp + " , " + zp + " , " + up + " , " + vp + " ) "
#print "get_function ( " + "Uuyvx )"
#print " "


xp = "o";yp = "o";zp = "o";up = "o";vp = "o"
# using get function Uxzvu

xp = ap
zp = bp
vp = cp
up = dp

#print " "
print "add_eval_obj( " + xp + " , " + yp + " , " + zp + " , " + up + " , " + vp + " ) "
#print "get_function ( " + "Uxzvu )"
#print " "


xp = "o";yp = "o";zp = "o";up = "o";vp = "o"
# using get function Uyuvz

yp = ap
up = bp
vp = cp
zp = dp

#print " "
print "add_eval_obj( " + xp + " , " + yp + " , " + zp + " , " + up + " , " + vp + " ) "
#print "get_function ( " + "Uyuvz )"
#print " "

















