import unittest
import sys, os

# Add the directory containing 'grendel' to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir))

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel import *
from grendel.util.units.errors import UnknownUnitError, IncompatibleUnitsError
from grendel.util.units.unit import PrefixedUnit
from grendel_tests import allow_generators

@allow_generators
class UnitPackageTest(unittest.TestCase):

    @classmethod
    def back_and_forth_test_generator(cls):
        def units_back_and_forth(self, unit1, unit2):
            # Check that they have the same superclass
            test_value = 3.14
            if unit1.__mro__[1] == unit2.__mro__[1]:
                # Convert from unit1 to unit 2 and back;  We should get test_value out
                self.assertAlmostEqual(
                    convert_units(convert_units(test_value, unit1, unit2), unit2, unit1),
                    test_value,
                    places=12
                )
            else:
                self.assertRaises(
                    IncompatibleUnitsError,
                    convert_units, test_value, unit2, unit1
                )
                self.assertRaises(
                    IncompatibleUnitsError,
                    convert, test_value, unit1, unit2
                )
        units_to_test = []
        for item in dir(util.units):
            klass = None
            try:
                klass = eval(item)
            except NameError:
                pass

            if isunit(klass) and not isinstance(klass, PrefixedUnit):
                units_to_test.append(klass)

        # Delete duplicates
        units_to_test = list(set(units_to_test))
        for unit1 in units_to_test:
            for unit2 in units_to_test:
                units_back_and_forth.__name__ = "test_" + str(unit1) + "_to_" + str(unit2) + "_and_back"
                yield units_back_and_forth, unit1, unit2

    #--------------------------------------------------------------------------------#

    def test_aliases_1(self):
        self.assertIs(Angstroms, Angstrom)

    def test_aliases_2(self):
        self.assertIs(AtomicUnitsOfCharge, AtomicUnitOfElectricCharge)

    def test_aliases_3(self):
        self.assertIs(MicroGram, Microgram)

    #--------------------------------------------------------------------------------#

    def test_unit_conv_errors_1(self):
        with self.assertRaises(UnknownUnitError):
            convert_units(3.14, Tensor, Angstroms)

    def test_unit_conv_errors_2(self):
        with self.assertRaises(UnknownUnitError):
            convert_units(3.14, Degrees, Atom)

    def test_unit_conv_errors_3(self):
        with self.assertRaises(UnknownUnitError):
            ValueWithUnits(1.5, Angstroms).in_units(float)

    def test_unit_conv_errors_4(self):
        with self.assertRaises(TypeError):
            Angstrom/Bohr in Angstrom/Bohr

    def test_unit_conv_errors_5(self):
        with self.assertRaises(TypeError):
            Angstrom/Bohr in Angstrom

    def test_unit_conv_errors_6(self):
        with self.assertRaises(IncompatibleUnitsError):
            (Angstrom/Attojoule).to(Meter/Attojoule**2)

    def test_unit_conv_errors_7(self):
        with self.assertRaises(UnknownUnitError):
            ValueWithUnits(10.5, str)

    #--------------------------------------------------------------------------------#

    def test_contains_1(self):
        self.assertIn(10.*Angstroms, Angstrom)

    def test_contains_2(self):
        self.assertNotIn(10.*Angstroms, Milliangstroms)

    def test_contains_3(self):
        self.assertNotIn(10.*Kilograms, Micrometers)

    def test_contains_4(self):
        self.assertNotIn(10.*Angstroms/Joule, Meters/Wavenumber)

    def test_contains_5(self):
        self.assertNotIn(-15.*Bohr/Microhartree, Meters)

    def test_contains_6(self):
        self.assertNotIn(-15.*Bohr/Microhartree, Kilograms)

    def test_contains_7(self):
        self.assertNotIn(-15.*Bohr/Microhartree, Bohrs/Year)

    def test_contains_8(self):
        self.assertNotIn(-15.*Bohr/Microhartree, Kilograms/Megahartree)

    def test_contains_9(self):
        self.assertNotIn(-15.*Bohr/Microhartree, Bohr/Hartree)

    def test_contains_10(self):
        self.assertNotIn(10.*Angstroms, Kilograms/Angstrom)

    def test_contains_11(self):
        self.assertIn(10.*Angstroms/Hartree, Angstrom/Hartrees)

    #--------------------------------------------------------------------------------#

    def test_common_sense_1(self):
        self.assertEqual(Angstrom, Angstroms)

    def test_common_sense_2(self):
        self.assertNotEqual(Joules, Millijoules)

    def test_common_sense_3(self):
        self.assertNotEqual(AtomicUnitsOfMass, Millijoules)

    def test_common_sense_4(self):
        self.assertNotEqual(AtomicUnitsOfEnergy, Millihartrees)

    def test_common_sense_5(self):
        self.assertEqual(AtomicUnitsOfEnergy, Hartrees)

    def test_common_sense_6(self):
        self.assertEqual((Bohr**2/Megabohr).to(Microbohr), 1.0)

    def test_common_sense_7(self):
        self.assertEqual((Microbohr/Bohr**2).to(Megabohr**(-1)), 1.0)

    def test_common_sense_8(self):
        self.assertEqual((Microbohr/(Bohr*Kilogram**2)).to(Megagram**(-2)), 1.0)

    def test_common_sense_9(self):
        self.assertEqual((Micrometer/Millennia).genre, CompositeUnit)

    def test_common_sense_10(self):
        self.assertNotEqual(Meter, 1.75)

    def test_composite_unit_1(self):
        self.assertIsInstance((Angstrom*Joule)/(Kilogram*Second), CompositeUnit)

    def test_composite_unit_2(self):
        self.assertEqual((Angstrom**2).name, "Angstrom**2")

    def test_composite_unit_3(self):
        self.assertEqual((Angstrom**-1).name, "1.0 / Angstrom")

    def test_composite_unit_4(self):
        self.assertEqual(repr(Angstrom**-2), "1.0 / Angstrom**2")

    def test_composite_unit_5(self):
        self.assertEqual((Angstrom**-1 * Kilogram**-1).name, "1.0 / (Angstrom * KiloGram)")

    def test_composite_unit_6(self):
        self.assertEqual((Bohr**2/Megabohr).reduced().name, "(1e-06 Bohr)")

    def test_composite_unit_7(self):
        self.assertEqual(((Angstrom**2*Joule)/(Kilogram*Second**-1))**-1, Kilogram/(Angstrom**2*Joule*Second))

    def test_composite_unit_8(self):
        self.assertEqual(((Angstrom**2*Joule)/(Kilogram*Second**-1))**2, Angstrom**4*Joule**2*Second**2/Kilogram**2)

    #--------------------------------------------------------------------------------#

    def test_composite_unit_to(self):
        expected = Angstroms.to(Meters)**2 * AMU.to(Kilogram)
        actual = (Angstroms**2 * AMU).to(Meter**2 * Kilogram)
        self.assertEqual(expected, actual)

    #--------------------------------------------------------------------------------#

    def test_value_with_units_1(self):
        self.assertEqual((10*Angstrom)*(10*Angstrom), 100*Angstrom**2)

    def test_value_with_units_2(self):
        self.assertEqual((10*Angstrom)*(Angstrom), 10*Angstrom**2)

    def test_value_with_units_3(self):
        self.assertEqual((10*Angstrom)*10, 100*Angstrom)

    def test_value_with_units_4(self):
        self.assertEqual(10*(10*Angstrom), 100*Angstrom)

    def test_value_with_units_5(self):
        self.assertEqual((100*Angstrom**2)/(10*Angstrom), 10*Angstrom)

    def test_value_with_units_6(self):
        self.assertEqual(10/(20.0 * Angstrom), 0.5*Angstrom**-1)

    def test_value_with_units_7(self):
        self.assertEqual(10 * (Angstrom/Joule), (Angstrom/Joule) * 10)

    def test_value_with_units_8(self):
        self.assertIsInstance(Millibohr**-1 / 35.0, ValueWithUnits)

    def test_value_with_units_9(self):
        self.assertIsInstance(10.0 / Angstrom, ValueWithUnits)

    def test_value_with_units_10(self):
        self.assertIsInstance(Angstrom * 10, ValueWithUnits)

    def test_value_with_units_11(self):
        self.assertEqual(35 / Angstrom**2, 35.0 * Angstrom**-2)

    def test_value_with_units_12(self):
        self.assertEqual(str(10*Angstrom), "10.0 Angstrom")

    def test_value_with_units_13(self):
        self.assertEqual((10*Angstrom).in_units(Kilometers), 1e-12 * Kilometers)

    def test_value_with_units_14(self):
        self.assertEqual(10*Angstrom + 15*Angstrom, 25*Angstrom)

    def test_value_with_units_15(self):
        self.assertEqual(17*Angstrom - 15*Angstrom, 2*Angstrom)

    #--------------------------------------------------------------------------------#

    def test_type_errors_1(self):
        with self.assertRaises(TypeError):
            class StupidUnit: __metaclass__ = Unit

    def test_type_errors_2(self):
        with self.assertRaises(TypeError):
            (Angstrom/KiloAtomicUnitOfEnergy).to(str)
    def test_type_errors_3(self):
        with self.assertRaises(TypeError):
            Angstrom**2 / "Hello World"

    def test_type_errors_4(self):
        with self.assertRaises(TypeError):
            Angstrom**2 * "Hello World"

    def test_type_errors_5(self):
        with self.assertRaises(TypeError):
            "Hello World" / Angstrom**2

    def test_type_errors_6(self):
        with self.assertRaisesRegexp(TypeError, r'unsupported operand'):
            10*Angstrom + Angstrom

    def test_type_errors_7(self):
        with self.assertRaisesRegexp(TypeError, r'unsupported operand'):
            Angstrom + 10*Angstrom

    def test_type_errors_8(self):
        with self.assertRaisesRegexp(TypeError, r'unsupported operand'):
            10*Angstrom - Angstrom

    def test_type_errors_9(self):
        with self.assertRaisesRegexp(TypeError, r'unsupported operand'):
            Angstrom - 10*Angstrom

