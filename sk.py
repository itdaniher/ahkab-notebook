from ahkab import *
import sympy
import pickle

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

	print r['ac'].keys()

# substitute the actual values to the symbols and plot
Vu1o, s, R1, R2, C1, C2, V1 = r['symbolic'][0].as_symbols("Vu1o s R1 R2 C1 C2 V1")
R1v, R2v, C1v, C2v, V1MAGv = 10e3, 10e3, 1e-9, 1e-9, 1.0
w = sympy.Symbol('w', real=True)
Vu1o = r['symbolic'][0]['Vu1o'].subs({V1:V1MAGv, R1:R1v, C1:C1v, C2:C2v, R2:R2v, s:w*sympy.I})
Vu1o_f = lambda x: Vu1o.subs({w:x})
Vu1o_symb_mag = numpy.array(map(sympy.N, map(sympy.Abs, map(Vu1o_f, r['ac']['w'][::50].tolist()[0]))), dtype='float')
Vu1o_symb_arg = numpy.array(map(sympy.N, map(sympy.arg, map(Vu1o_f, r['ac']['w'][::50].tolist()[0]))), dtype='float')


import pylab
import numpy as np

pylab.subplot(211)
pylab.semilogx(r['ac']['w'].T, 20*np.log10(r['ac']['|Vu1o|'].T), '-', label='From AC simulation')
pylab.semilogx(r['ac']['w'][::50].T, 20*np.log10(Vu1o_symb_mag), 'v', label='From SYMB simulation')
pylab.ylabel("|VOUT| (dB)")
pylab.legend()
pylab.subplot(212)
pylab.semilogx(r['ac']['w'].T, (180.0/np.pi*r['ac']['arg(Vu1o)'].T), '-', label='From AC simulation')
pylab.semilogx(r['ac']['w'][::50].T, 180.0/np.pi*Vu1o_symb_arg, 'v', label='From SYMB simulation')
pylab.ylabel("arg(VOUT)/deg")
pylab.xlabel("$\omega$ (rad/s)")
pylab.legend()
print r['symbolic'][0]
pylab.show()
