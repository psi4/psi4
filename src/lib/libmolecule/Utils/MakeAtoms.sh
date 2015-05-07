#!/bin/sh
atoms=( "Hydrogen" "Helium" "Lithium" "Beryllium" "Boron" "Carbon" "Nitrogen" "Oxygen" "Fluorine" "Neon" "Sodium" "Magnesium" "Aluminium" "Silicon" "Phosphorus" "Sulfur" "Chlorine" "Argon" "Potassium" "Calcium" "Scandium" "Titanium" "Vanadium" "Chromium" "Manganese" "Iron" "Cobalt" "Nickel" "Copper" "Zinc" "Gallium" "Germanium" "Arsenic" "Selenium" "Bromine" "Krypton" "Rubidium" "Strontium" "Yttrium" "Zirconium" "Niobium" "Molybdenum" "Technetium" "Ruthenium" "Rhodium" "Palladium" "Silver" "Cadmium" "Indium" "Tin" "Antimony" "Tellurium" "Iodine" "Xenon"
 "Cesium" "Barium" "Lanthanum" "Cerium" "Praseodymium" "Neodymium" "Promethium" "Samarium" "Europium" "Gadolinium" "Terbium" "Dysprosium" "Holmium" "Erbium" "Thulium" "Ytterbium" "Lutetium" "Hafnium" "Tantalum" "Tungsten" "Rhenium" "Osmium" "Iridium" "Platinum" "Gold" "Mercury" "Thallium" "Lead" "Bismuth" "Polonium" "Astatine" "Radon" "Francium" "Radium" "Actinium" "Thorium" "Protactinium" "Uranium" "Neptunium" "Plutonium" "Americium" "Curium" "Berkelium" "Californium" "Einsteinium" "Fermium" "Mendelevium" "Nobelium" "Lawrencium" "Rutherfordium" "Dubnium" "Seaborgium" "Bohrium" "Hassium" "Meitnerium" "Darmstadtium" "Roentgenium" "Copernicium" "Ununtrium" "Flerovium" "Ununpentium" "Livermorium" "Ununseptium" "Ununoctium" )

#for i in ${atoms[@]};do
#  printf "${i^^},"
#done
#for i in ${atoms[@]};do
#  echo "temp[${i^^}]=\"Elemental ${i}\";"
#done
for i in ${atoms[@]};do
    echo "class ${i}:public FxnalGroup{"
    echo "public:"
    echo "${i}():FxnalGroup(${i^^},0){}"
    echo "${i}(const int* Members):FxnalGroup(${i^^},0,1,Members){}"
    echo "};"
done
#for i in ${atoms[@]};do
#   printf "${i},"
#done
