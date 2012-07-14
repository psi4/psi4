# efpmd

## Input file format

##### Specify the type of the simulation

`run_type [sp|grad|cg|nve|nvt]`

`sp` - single point calculation.

`grad` - gradient calculation.

`cg` - geometry optimization using conjugate gradient method.

`nve` - molecular dynamics in microcanonical ensemble.

`nvt` - molecular dynamics in canonical ensemble.

##### Specify the format of fragment input

`coord [xyzabc|points]`

`xyzabc` - Coordinates of the center of mass and Euler angles.

`points` - Coordinates of three atoms for each fragment.

##### Specify the units to use

`units [angs|bohr]`

`angs` - Units are Angstroms.

`bohr` - Units are Bohr.

##### Specify which energy terms to include in EFP computation

`terms [elec [pol [disp [xr]]]]`

`elec` - Include electrostatics energy.

`pol` - Include polarization energy.

`disp` - Include dispersion energy.

`xr` - Include exchange repulsion energy.

##### Specify electrostatic damping type

`elec_damp [screen|overlap|off]`

`screen` - Damping formula based on SCREEN group in the EFP potential.

`overlap` - Overlap based damping which computes charge penetration energy.

`off` - No electrostatic damping.

##### Specify dispersion damping type

`disp_damp [tt|overlap|off]`

`tt` - Damping based on the formula by Tang and Toennies.

`overlap` - Overlap-based dispersion damping.

`off` - No dispersion damping.

##### Specify the path to the directory with fragment library

`fraglib_path <path>`

The `<path>` parameter should not contain spaces or be in double quotes.

##### Specify the path to the directory with user-created fragments

`userlib_path <path>`

The `<path>` parameter should not contain spaces or be in double quotes.

##### Fragment input

One or more `fragment` groups. Each `fragment <name>` line is followed by
either a line of six numbers if `coord` is `xyzabc` or three lines with three
numbers on each line if `coord` is `points`. If `<name>` contains an `_l`
suffix the fragment parameters will be searched in the `fraglib_path`
directory. Otherwise the directory specified by the `userlib_path` option will
be used.
