import sympy
import numpy as np

def getMapping(tf):
	syms = tf.atoms(sympy.Symbol)
	return dict(zip(map(str, syms), syms))

def reduceTF(tf, cir):
	for e in cir:
		if not e.is_nonlinear:
			name = e.part_id
			obs = getMapping(tf)
			if name in map(str, tf.atoms()):
				tf = tf.subs(obs[name], float(e.value))
	return tf

def tustin(tf, sRate):
	s =	getMapping(tf)['s']
	T = 1/sRate
	poles = sympy.roots(sympy.denom(tf), s, multiple=True)
	centroid = np.mean(np.abs(poles))
	invz = sympy.Symbol('invz', real=True)
	# normalized center frequency derived from poles of filter SISO transfer function
	w0 = 2*np.pi*centroid/(sRate*2*np.pi)
	# modified bilinear transform w/ frequency warping
	bt = w0/sympy.tan(w0*T/2) * ((1-invz)/(1+invz))
	dt = sympy.simplify(tf.subs({s: bt}))
	b = sympy.Poly(sympy.numer(dt)).all_coeffs()[::-1]
	a = sympy.Poly(sympy.denom(dt)).all_coeffs()[::-1]
	normalize = lambda x: float(x/a[0])
	return (map(normalize, b), map(normalize, a))

def e2nd(expression):
	""" basic helper function that accepts a sympy expression, expands it,
	attempts to simplify it, and returns a numerator and denomenator pair for the instantiation of a scipy
	LTI system object. """
	expression = expression.expand()
	expression = expression.cancel()
	n = sympy.Poly(sympy.numer(expression), s).all_coeffs()
	d = sympy.Poly(sympy.denom(expression), s).all_coeffs()
	n = map(float, n)
	d = map(float, d)
	return (n, d)
