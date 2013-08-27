import numpy as np
import sympy
from ahkab import *
import pickle
import pylab

mycircuit = circuit.circuit(title="Tow-Thomas biquad")

gnd = mycircuit.get_ground_node()

def buildsvf(svf):
	ar = svf.add_resistor
	ac = svf.add_capacitor
	ao = lambda name, p, n: svf.add_vcvs(name, "U"+name[1:]+"o", gnd, p, n, 1e6)
	ar("R1", "in", "U1n", 10e3)
	ar("R3", "U1n", "U1o", 10e3)
	ar("R2", "U3o", "U1n", 10e3)
	ar("R4", "U1o", "U2n", 10e3)
	ar("R5", "U2o", "U3n", 10e3)
	ar("R6", "U3n", "U3o", 10e3)
	ac("C1", "U1n", "U1o", 15e-9)
	ac("C2", "U2n", "U2o", 15e-9)
	ao("E1", gnd, "U1n")
	ao("E2", gnd, "U2n")
	ao("E3", gnd, "U3n")

buildsvf(mycircuit)

printing.print_circuit(mycircuit)

mycircuit.add_vsource(name="V1", ext_n1="in", ext_n2=gnd, vdc=5, vac=1)

subs = symbolic.parse_substitutions(('R2=R1', 'R3=R1', 'C2=C1', 'E2=E1', 'E3=E1', "R3=R1", "R4=R1", "R5=R1", "R6=R1"))
print subs

symbolic_sim = {"type":"symbolic", "ac": True, "source":None, 'subs':subs}

ac_sim = {'type':'ac', 'start':0.1, 'stop':100e6, 'nsteps':1000}

an_list = [ac_sim, symbolic_sim]

try:
	r = pickle.load(open("results-ttb.pk"))
except:
	r = ahkab.process_analysis(an_list, circ=mycircuit, outfile="./.ahkab_data", verbose=2, cli_tran_method=None, guess=True, disable_step_control=False)
	pickle.dump(r, open("results-ttb.pk", "wb"))

# bandpass output
output = r['symbolic'][0].as_symbol('VU1o')
tf = r['symbolic'][0][output]
syms = filter(lambda x: x.is_Symbol, tf.atoms())
locals().update(
	dict(zip(map(str, syms), syms))
)

print "The band-pass output has expression: \n%s = %s"%(output.name, str(tf))
# get all the other symbols we need
w = sympy.Symbol('w', real=True)
tf = sympy.limit(tf, E1, sympy.oo, '+')

v1v, R1v, C1v = (1, 10e3, 15e-9)
print "We set: v1=%d (AC), R1=%g ohm, C1=%g F." % (v1v, R1v, C1v)
tf = tf.subs({R1:R1v, C1:C1v, V1:v1v})
print tf

sRate = 44.1e3
T = 1/sRate

fs = lambda x: sympy.N(abs(tf.subs({s:sympy.I*x})))
ws = np.logspace(1, np.log10(sRate), 5e2)
mags = map(fs,ws)

def tustin(expression):
	poles = sympy.roots(sympy.denom(expression), s, multiple=True)
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

(b, a) = tustin(tf)
print b,a

pylab.semilogx(ws, map(fs, ws), 'v', label="from transfer function")
aclabel = "|" + output.name + "|"
pylab.semilogx(r['ac']['w'].T[::10], r['ac'][aclabel].T[::10], '-', label='from AC simulation')
pylab.vlines(np.abs(sympy.roots(sympy.denom(tf), s, multiple=True)), 0, 1, 'r')

# build white noise input vector, normalized to \pm 1.0
x = list(2*np.random.random(sRate)-1.0)
# allocate output vector for y[n-2] indexing to work
y = [0.0]*len(x)

# Direct Form I
for n in range(len(x)):
	y[n] = b[0]*x[n] + b[1]*x[n-1] + b[2]*x[n-2] - a[1]*y[n-1] - a[2]*y[n-2]

# frequencies to radians
freqs = (np.fft.fftfreq(len(y), T))[0:len(y)/2] * (2*np.pi)
mags = 	(np.abs(np.fft.fft(np.array(y))))[0:len(y)/2]

# normalize magnitude to 1.0
mags = mags/max(mags)
pylab.semilogx(freqs, mags, '.', label='DFT of IIR filtered white noise')
pylab.xlabel('freq (w/s)')
pylab.ylabel('normalized amplitude')
pylab.legend()
pylab.show()
