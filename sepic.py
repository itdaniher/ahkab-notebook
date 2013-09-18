import numpy as np
from ahkab import *
import sympy
import pickle
import ahkabHelpers
mycircuit = circuit.circuit(title="sepic", filename=None)

# https://en.wikipedia.org/wiki/File:SEPIC_Schematic.gif

#.MODEL 2N7000 NMOS (LEVEL=3 RS=0.205 NSUB=1.0E15 DELTA=0.1 KAPPA=0.0506 TPG=1 CGDO=3.1716E-9 RD=0.239 VTO=1.000 VMAX=1.0E7 ETA=0.0223089 NFS=6.6E10 TOX=1.0E-7 LD=1.698E-9 UO=862.425 XJ=6.4666E-7 THETA=1.0E-5 CGSO=9.09E-9 L=2.5E-6 W=0.8E-2)

def buildsepic(circuit):
	gnd = circuit.get_ground_node()
	ar = circuit.add_resistor
	ac = circuit.add_capacitor
	al = circuit.add_inductor
	av = circuit.add_vsource
#	circuit.add_model("diode", "aDiodeModel", {"name": "aDiodeElement"})
#	circuit.add_diode("D1", "da", "dc", "aDiodeModel")
	av("Vin", "in", gnd, 2.5, 1.0)
	ac("Cin", "in", gnd, 100e-3)
	al("L1", "in", "s1", 47e-6)
	ac("C1", "s1", "da", 10e-6)
	al("L2", "da", gnd, 47e-6)
	ac("C2", "dc", gnd, 100e-3)
	#circuit.add_inductor_coupling("K1", "L1", "L2", 0.99)

buildsepic(mycircuit)

printing.print_circuit(mycircuit)

symbolic_sim = ahkab.new_symbolic()

try:
	r = pickle.load(open("results-sepic.pk"))
except:
	r = ahkab.run(mycircuit, (symbolic_sim,))
	pickle.dump(r, open("results-sepic.pk", "wb"))

print mycircuit
print r['symbolic'][0]

tf = r['symbolic'][0]['I[L1]']

tf = ahkabHelpers.reduceTF(tf, mycircuit)
import compensators
print compensators.drawBode(tf)
from pylab import show
show()
