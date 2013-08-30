import numpy as np
from ahkab import *
import sympy
import pickle
import ahkabHelpers

mycircuit = circuit.circuit(title="sallen-key 2pole lpf", filename=None)

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

printing.print_circuit(mycircuit)

symbolic = {"type":"symbolic", "ac": True, "source":"V1"}

ac_sim = {'type':'ac', 'start':10, 'stop':100e6, 'nsteps':1000}

try:
	r = pickle.load(open("results-sallenkey.pk"))
except:
	r = ahkab.process_analysis(an_list=[symbolic, ac_sim], circ=mycircuit, outfile="/tmp/ahkab_data", verbose=2, cli_tran_method=None, guess=True, disable_step_control=False)
	pickle.dump(r, open("results-sallenkey.pk", "wb"))

# substitute the actual values to the symbols and plot
tf = ahkabHelpers.reduceTF(r['symbolic'][0]['Vu1o'], mycircuit)
evalTF = sympy.lambdify(ahkabHelpers.getMapping(tf)['s'], tf)

out = map(evalTF, 1j*numpy.array(r['ac']['w'][::50].tolist()[0]))

Vu1o_symb_mag = np.abs(out)
Vu1o_symb_arg = np.angle(out, deg=True)
import pylab

pylab.subplot(211)
pylab.semilogx(r['ac']['w'].T, 20*np.log10(r['ac']['|Vu1o|'].T), '-', label='From AC simulation')
pylab.semilogx(r['ac']['w'][::50].T, 20*np.log10(Vu1o_symb_mag), 'v', label='From SYMB simulation')
pylab.ylabel("|VOUT| (dB)")
pylab.legend()
pylab.subplot(212)
pylab.semilogx(r['ac']['w'].T, (180.0/np.pi*r['ac']['arg(Vu1o)'].T), '-', label='From AC simulation')
pylab.semilogx(r['ac']['w'][::50].T, Vu1o_symb_arg, 'v', label='From SYMB simulation')
pylab.ylabel("arg(VOUT)/deg")
pylab.xlabel("$\omega$ (rad/s)")
pylab.legend()
print r['symbolic'][0]
pylab.show()
