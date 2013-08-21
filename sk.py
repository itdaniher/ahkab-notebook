from ahkab import *
import pickle

mycircuit = circuit.circuit(title="sallen-key 2pole lpf", filename=None)

# https://upload.wikimedia.org/wikipedia/commons/5/5c/Sallen-Key_Lowpass_Example.svg

def buildsk(svf):
	gnd = svf.get_ground_node()
	ar = svf.add_resistor
	ac = svf.add_capacitor
	al = svf.add_inductor
	av = svf.add_vsource
	ao = lambda name, p, n: svf.add_vcvs(name, name+"o", gnd, p, n, 1e6)
	ao("u1", "u1p", "u1o")
	ar("R1", "in", "R2", 10e3)
	ar("R2", "R2", "u1p", 10e3)
	ac("C1", "u1o", "R2", 1e-9)
	ac("C2", "u1p", gnd, 1e-9)
	av("V1", "in", gnd, 2.5, 1.0)

buildsk(mycircuit)

printing.print_circuit(mycircuit)

symbolic = {"type":"symbolic", "ac": True, "source":"V1"}

ac_sim = {'type':'ac', 'start':10, 'stop':100e6, 'nsteps':1000}


try:
	r = pickle.load(open("results-ttb.pk"))
except:
	r = ahkab.process_analysis(an_list=[symbolic, ac_sim], circ=mycircuit, outfile="/tmp/ahkab_data", verbose=2, cli_tran_method=None, guess=True, disable_step_control=False)
	pickle.dump(r, open("results-ttb.pk", "wb"))

	print r['ac'].keys()

import pylab
import numpy as np

pylab.subplot(211)
pylab.semilogx(r['ac']['w'].T, 20*np.log10(r['ac']['|Vu1o|'].T))
pylab.ylabel("|VOUT| (dB)")
pylab.subplot(212)
pylab.semilogx(r['ac']['w'].T, (180.0/np.pi*r['ac']['arg(Vu1o)'].T))
pylab.ylabel("arg(VOUT)/deg")
pylab.xlabel("$\omega$ (rad/s)")
print r['symbolic'][0]
pylab.show()
