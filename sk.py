import numpy as np
import ahkab
from ahkab import circuit, printing
import sympy
import pickle
import ahkabHelpers

mycircuit = circuit.Circuit(title="sallen-key 2pole lpf", filename=None)

# https://upload.wikimedia.org/wikipedia/commons/5/5c/Sallen-Key_Lowpass_Example.svg

def buildsk(svf):
	gnd = svf.get_ground_node()
	ar = svf.add_resistor
	ac = svf.add_capacitor
	al = svf.add_inductor
	av = svf.add_vsource
	ao = lambda name, p, n: svf.add_vcvs("E"+name[1:], name+"o", gnd, p, n, 1e6)
	ao("u1", "u1p", "u1o")
	ar("R1", "in", "nR2", 10e3)
	ar("R2", "nR2", "u1p", 10e3)
	ac("C1", "u1o", "nR2", 1e-9)
	ac("C2", "u1p", gnd, 1e-9)
	av("V1", "in", gnd, 2.5, 1.0)

buildsk(mycircuit)

print(mycircuit)

symbolic_sim = ahkab.new_symbolic(ac_enable=True, source='V1')

ac_sim = ahkab.new_ac(start=10, stop=100e6, points=1000, x0=None)

try:
	r = pickle.load(open("results-sallenkey.pk"))
except:
	ahkab.queue(ac_sim, symbolic_sim)
	r = ahkab.run(mycircuit)
	pickle.dump(r, open("results-sallenkey.pk", "wb"))

print r['symbolic'][0]
print "Symbolic transfer functions:"
printing.print_symbolic_transfer_functions(r['symbolic'][1])

# substitute the actual values to the symbols and plot
V1 = r['symbolic'][0].as_symbol('V1')
tf = r['symbolic'][0]['Vu1o']/V1
tf = ahkabHelpers.reduceTF(tf, mycircuit)
evalTF = sympy.lambdify(ahkabHelpers.getMapping(tf)['s'], tf)

out = map(evalTF, 1j*r['ac']['w'][::50])

Vu1o_symb_mag = np.abs(out)
Vu1o_symb_arg = np.angle(out, deg=True)
import pylab

pylab.subplot(211)
pylab.semilogx(r['ac']['w'], np.abs(20*np.log10(r['ac']['Vu1o'])), '-', label='From AC simulation')
pylab.semilogx(r['ac']['w'][::50], 20*np.log10(Vu1o_symb_mag), 'v', label='From SYMB simulation')
pylab.ylabel("|VOUT| (dB)")
pylab.legend()
pylab.subplot(212)
pylab.semilogx(r['ac']['w'], (180.0/np.pi*np.angle(r['ac']['Vu1o'])), '-', label='From AC simulation')
pylab.semilogx(r['ac']['w'][::50], Vu1o_symb_arg, 'v', label='From SYMB simulation')
pylab.ylabel("arg(VOUT)/deg")
pylab.xlabel("$\omega$ (rad/s)")
pylab.legend()
pylab.savefig('sk.png')
