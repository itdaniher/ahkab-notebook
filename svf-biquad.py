from ahkab import *
import pickle
import sympy

mycircuit = circuit.Circuit(title="state variable filter")

gnd = mycircuit.get_ground_node()

def buildsvf(svf):
	ar = svf.add_resistor
	ac = svf.add_capacitor
	al = svf.add_inductor
	ao = lambda name, p, n: svf.add_vcvs(name, "U"+name[1:]+"o", gnd, p, n, 1e6)
	ar("Rin", gnd, "U1p", 4.7e3)
	ar("R00", "U1p", "U2o", 10e3)
	ar("R01", "in", "U1n", 10e3)
	ar("Rf1", "U1o", "U2n", 10e3)
	ar("R10", "U1o", "U1n", 10e3)
	ar("R11", "U1n", "U3o", 10e3)
	ar("Rf2", "U2o", "U3n", 10e3)
	ar("R02", "U2o", "out", 10e3)
	ac("C10", "U2o", "U2n", 15e-9)
	ac("C11", "U3o", "U3n", 15e-9)
	ao("E1", "U1p", "U1n")
	ao("E2", gnd, "U2n")
	ao("E3",  gnd, "U3n")

buildsvf(mycircuit)
mycircuit.add_vsource(part_id="V1", n1="in", n2=gnd, dc_value=5, ac_value=1)

printing.print_circuit(mycircuit)

subs = symbolic.parse_substitutions(('E2=E1', 'E3=E1', 'R01=R00', 'R02=R00', 'R11=R00', 'R10=R00', 'C11=C10', 'Rf2=Rf1', 'Rin=R00'))

symbolic_sim = ahkab.new_symbolic(ac_enable=True, source="V1", subs=subs)
ac_sim = ahkab.new_ac(start=0.1, stop=100e6, points=1000, x0=None)

try:
	r = pickle.load(open("results-svf.pk"))
except:
	ahkab.queue(symbolic_sim, ac_sim)
	r = ahkab.run(mycircuit)
	pickle.dump(r, open("results-svf.pk", "wb"))

E = r['symbolic'][0].as_symbol('E1')

out = sympy.limit(r['symbolic'][0]['VU1o'], E, sympy.oo, '+')

print "VU1o =", out.simplify()
