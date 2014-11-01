# coding=utf-8
from __future__ import division
import sympy
from sympy.abc import s
from matplotlib import pyplot
import sympy.mpmath as mpmath
import numpy
import scipy.signal as signal
import pprint
j = sympy.I

ct = 1000

mpmath.mp.dps = 40;
mpmath.mp.pretty = True

blacks = lambda G_s, H_s: G_s/(1+G_s*H_s)

def eeval(expression, w):
	""" evaluate a sympy expression at omega. return magnitude, phase."""
	num, den = e2nd(expression)
	y = numpy.polyval(num, 1j*w) / numpy.polyval(den, 1j*w)
	phase = numpy.arctan2(y.imag, y.real) * 180.0 / numpy.pi
	mag = abs(y)
	return mag, phase

def bode(expression, n = 10):
	decibels = lambda lin: 20*numpy.log10(numpy.abs(lin))
	num, den = e2nd(expression)
	freqs = signal.findfreqs(num, den, n)
	magnitude = numpy.array([])
	phase = numpy.array([])
	for freq in freqs:
		(m, p) = eeval(expression, freq)
		magnitude = numpy.append(magnitude, m)
		if p >= 0:
			p = -360+p
		phase = numpy.append(phase, p)
	magnitude = numpy.array(map(decibels, magnitude))
	return freqs, magnitude, phase

def e2nd(expression):
	""" basic helper function that accepts a sympy expression, expands it,
	attempts to simplify it, and returns a numerator and denomenator pair for the instantiation of a scipy
	LTI system object. """
	expression = expression.expand()
	expression = expression.cancel()
	#print expression
	n = sympy.Poly(sympy.numer(expression), s).all_coeffs()
	d = sympy.Poly(sympy.denom(expression), s).all_coeffs()
	n = map(float, n)
	d = map(float, d)
	return (n, d)

def pade(t, n):
	""" pade approximation of a time delay, as per mathworks' controls toolkit.
	supports arbitrary precision mathematics via mpmath and sympy"
	for more information, see:
		* http://home.hit.no/~hansha/documents/control/theory/pade_approximation.pdf
		* http://www.mathworks.com/help/control/ref/pade.html
	"""
	# e**(s*t) -> laplace transform of a time delay with 't' duration
	# e**x -> taylor series
	taylor = mpmath.taylor(sympy.exp, 0, n*2)
	(num, den) = mpmath.pade(taylor, n, n)
	num = sum([x*(-t*s)**y for y,x in enumerate(num[::-1])])
	den = sum([x*(-t*s)**y for y,x in enumerate(den[::-1])])
	return num/den

def sys2e(system):
	""" basic helper function that accepts an instance of the scipy.signal.lti class
	and returns a sympy expression given by the ratio of the expanded numerator and denomenator"""
	den = sum([x*s**y for y,x in enumerate(system.den[::-1])])
	num = sum([x*s**y for y,x in enumerate(system.num[::-1])])
	return sympy.factor(num/den)

def findZero(arr):
	return numpy.argmin(numpy.abs(arr))

def phaseMargin(expression):
	w, mag, phase = bode(expression, n=ct)
	crossingPoint = findZero(mag) # point where magnitude is closest to 0dB
	return {"w_c": w[crossingPoint], "p_m": phase[crossingPoint]+180}

def gainMargin(expression):
	w, mag, phase = bode(expression, n=ct)
	unstablePoint = findZero(-phase-180)
	return {"w": w[unstablePoint], "g_m": -mag[unstablePoint]}

def steadyStateError(expression):
	""" use sympy to determine limit at 0. """
	errors = {}
	for order, name in enumerate(['dc', 'step', 'ramp', 'parabola']):
		errors[name] = 1/(1+sympy.limit(s*expression*1/s**order, s, 0))
	return errors #{"sse": 1/(1+fvt)}

def reducedGainCompensate(expression, target):
	""" find the uncompensated system's magnitude at the right phase to give the right margin, normalize against it"""
	w, mag, phase = bode(expression, n=ct)
	phaseTarget = (-180 + target)
	targetIndex = findZero(phaseTarget - phase)
	magnitude = mag[targetIndex]
	magnitude = 10**(magnitude/20)
	return {"K": 1/magnitude}, 1/magnitude

def dominatePoleCompensate(expression, target):
	""" add a pole, call reducedGainCompensate"""
	K = reducedGainCompensate(1/s*expression, target)[0]['K']
	return {"K": K}, K/s

def lagCompensate(expression, target):
	w, mag, phase = bode(expression, n=ct)
	phaseTarget = -180 + (target + 6) # find location of new crossing goal; extra '-6' is from 10/wc rule
	targetIndex = findZero(phaseTarget - phase)
	w_c = w[targetIndex]
	tau = 10/w_c
	# equation for basic lag compensator
	alpha = eeval(expression, w_c)[0]
	G_c = (tau*s+1)/(alpha*tau*s+1)
	return {"alpha": alpha, "tau": tau}, G_c

def leadCompensate(expression, target):
	# alpha = 10 for 55 degrees phase margin
	w, mag, phase = bode(expression, n=ct)
	alpha = 10
	phaseTarget = -(180 + (55 - target))
	targetIndex = findZero(phase-phaseTarget)
	w_c = w[targetIndex]
	tau = 1 / (numpy.sqrt(alpha) * w_c)
	G_c = (alpha*tau*s+1)/(tau*s+1)
	# get magnitude of loop transfer function L(s) at w_c
	K_l = 1/ eeval(G_c*expression, w_c)[0]
	# equation for basic lead compensator
	G_c = K_l * (alpha*tau*s+1)/(tau*s+1)
	return {"K_l": K_l, "alpha": alpha, "tau":tau}, G_c

def drawBode(expression, f1=pyplot.figure(), color="k", labeled=""):
	""" draw bode plot for a given scipy.signal.lti instance. """
	adjustprops = dict(left=0.1, bottom=0.1, right=0.97, top=0.93, hspace=0.2)
	f1.subplots_adjust(**adjustprops)
	magPlot = f1.add_subplot(2, 1, 1)
	phasePlot = f1.add_subplot(2, 1, 2, sharex=magPlot)
	w, mag, phase = bode(expression, n=ct)
	magPlot.semilogx(w, mag, label=labeled, color=color)
	magPlot.set_title("magnitude")
	magPlot.set_ylabel("amplitude (dB)")
	phasePlot.semilogx(w, phase, label=labeled, color=color)
	phasePlot.set_title("phase")
	phasePlot.set_xlabel("radians per second")
	phasePlot.set_ylabel("degrees")
	pyplot.setp(magPlot.get_xticklabels(), visible=False)
	pyplot.legend(loc="best")
	f1.savefig("./figure1.png", dpi = 300)
	#f2 = pyplot.figure()
	#nichols_grid()
	#pyplot.plot(phase,mag)
	#f2.savefig("./figure2.png", dpi = 300)
	return f1#, f2

def rLocus(expression, f1=pyplot.figure(), color="k", labeled=""):
	r, k = control.matlab.rlocus(sys=control.TransferFunction(*e2nd(expression)), klist = numpy.linspace(0, 10, 10000))
	pyplot.plot(k, r, color)
	return f1
