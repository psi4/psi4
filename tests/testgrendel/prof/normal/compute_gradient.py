from itertools import combinations_with_replacement
import sys, os

sys.path.append('../../')
from grendel import *


mol = Molecule("""
 H         -0.8001600401       -0.8822498902        0.4394135499
 O         -0.0148376410       -0.6961963579       -0.0579905499
 O          0.0148376410        0.6961963579       -0.0579905499
 H          0.8001600401        0.8822498902        0.4394135499
""", units=Angstroms
)

# Doesn't work right now...
#mol.convert_units(Bohr)

rep = InternalRepresentation(mol, """
    bond 2 1
    bond 3 2
    bend 1 2 3
    bond 3 4
    bend 2 3 4
    tors 1 2 3 4
""",
    units={
        DistanceUnit: Angstrom,
        AngularUnit: Radians
    }
)


data_getter = LegacyXMLResultGetter(rep, 'data.xml')

mol.use_result_getter(data_getter)

print(mol.internal_representation)

dm = DisplacementManager(
    mol.internal_representation,
    Energy,
    deltas = [0.01*Angstrom, 0.01*Angstrom, 0.02*Radians, 0.01*Angstrom, 0.02*Radians, 0.02*Radians]
)

rep = mol.internal_representation
ncoords = len(mol.internal_representation)

ff = FiniteDifferenceForceField(
    max_order=4,
    robustness=2,
    representation=rep,
    property=Energy,
    property_units=Hartrees,
    deltas = [0.01*Angstrom, 0.01*Angstrom, 0.02*Radians, 0.01*Angstrom, 0.02*Radians, 0.02*Radians],
    run_now=True
)

#print(ff.for_order(1))
#print(ff.for_order(2))
#print(ff.for_order(2).formatted_string())
#print(MatrixFormatter.quick_format(ff.for_order(2).view(Tensor) * Hartrees.to(Attojoules), float_digits=3, float_type_string='f'))

funit = (Hartrees/Angstrom**3).to(Joules/(Meter*Angstrom**2)) * 10**5 * 10**3
##FUNIT = 4.359813653 / BohrRadius.in_units(Angstroms)**2
#F = ff.for_order(2)
#ffcart = ff.in_representation(mol.cartesian_representation)
#FX = ffcart.for_order(2).view(Matrix)
#FX = FX*(Angstroms**-2).to(Bohr**-2)
#print(FX.formatted_string(name="FX", float_digits=7, float_type_string='f'))
#sqm = mol.inverse_sqrt_mass_matrix
#FXM = FX.transformed(sqm) * funit
#print(FXM.formatted_string(name="FXM", float_digits=7, float_type_string='f'))
#evals, evecs = FXM.eigensystem()
##evals, evecs = np.real(evals), np.real(evecs)
#print([e for e in evals])
#print(evecs.formatted_string(name="LXM", float_digits=7, float_type_string='f'))
#LX = sqm * evecs
#print(LX.formatted_string(name="LX", float_digits=7, float_type_string='f'))
#print((funit*FX.transformed(LX.T)).formatted_string(name="Lambda", float_digits=7, float_type_string='f'))
ff.for_order(1)[...] = 0
nrep = NormalRepresentation(ff)
print(nrep.frequencies)
ffcart = ff.in_representation(mol.cartesian_representation)
print(chopped(ffcart.for_order(2).view(Matrix).transformed(nrep.b_matrix)).formatted_string())
ffnorm = ff.in_representation(nrep)
f3x = ffnorm.for_order(3)
f4x = ffnorm.for_order(4)
#f3x[...] = f3x*funit
f3x[...] = NormalRepresentation.convert_to_wavenumbers(f3x, Hartrees)
print(f3x.formatted_string(name="F3", float_digits=3, float_type_string='f'))
f4x[...] = NormalRepresentation.convert_to_wavenumbers(f4x, Hartrees)
print(f4x.formatted_string(name="F4", float_digits=3, float_type_string='f'))
#print(nrep.l_matrix.formatted_string(
#    name="L matrix",
#    float_digits=3,
#    float_type_string='f')
#)
#Linv = np.linalg.pinv(nrep.l_matrix)
#print(Linv.formatted_string(
#    name="L inv",
#    float_digits=3,
#    float_type_string='f')
#    )
#G = rep.g_matrix
#print(MatrixFormatter.quick_format(G, float_digits=3, float_type_string='f'))
#
#G_12 = G.sqrt_matrix()
#GFG = F.transformed(G_12)
#print(MatrixFormatter.quick_format(GFG* Hartrees.to(Attojoules), float_digits=3, float_type_string='f'))
#evals, evecs = GFG.eigensystem()
#print(evals*Hartrees.to(Attojoules))
#print([ev.view(Vector).norm() for ev in evecs])
#evecs = G_12 * Matrix(evecs).T
#evals, evecs = zip(*sorted((val, vec) for val, vec in zip(evals, evecs.col_iter)))
##evals = Vector(evals) * Hartrees.to(Joules) / (AMU * mol.cartesian_representation.units**2).to(Kilograms * Meters**2)
##evals = Vector(evals) * 2.37492753340298e+28 / Angstroms.to(Bohr)**2
#evals = Vector(evals) * (Hartrees**(1.0/2.0)/(AMU*Bohr)).to(Joules**(1.0/2.0)/(Kilograms*Meters))
##evals = Vector(evals) * Kilograms.to(AMU)*Hartrees.to(Joules)**(1.0/2.0)*Meters.to(Bohr)
#evecs = Matrix(evecs).T
#frequencies = Vector([math.sqrt(val) if val > 0 else -math.sqrt(-val) for val in evals])
#print(frequencies * Hertz.to(Wavenumbers) / (rep.units[DistanceUnit]).to(Bohr))
#print([ev.norm() for ev in evecs.row_iter])
##print(frequencies * Hertz.to(Wavenumbers))

