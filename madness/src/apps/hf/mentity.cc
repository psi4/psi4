/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id: test.cc 257 2007-06-25 19:09:38Z HartmanBaker $
*/


/// \file mentity.cc
/// \brief Simple management of molecular information and potential

#include <mra/mra.h>
#include <constants.h>
#include "mentity.h"
#include <misc/misc.h>

using namespace madness;

static const double PI = 3.1415926535897932384;

static inline double distance(double x1, double y1, double z1, double x2, double y2, double z2) {
    double xx = x1-x2;
    double yy = y1-y2;
    double zz = z1-z2;
    return sqrt(xx*xx + yy*yy + zz*zz);
}

static inline double distance_sq(double x1, double y1, double z1, double x2, double y2, double z2) {
    double xx = x1-x2;
    double yy = y1-y2;
    double zz = z1-z2;
    return xx*xx + yy*yy + zz*zz;
}

static const unsigned int NUMBER_OF_ATOMS_IN_TABLE = 110;
static const AtomicData atomic_data[NUMBER_OF_ATOMS_IN_TABLE] = {
    {"Bq",  "bq",   0  ,  0   ,  0.0               , 0.0           ,0.0             , 0.0    },
    {"H",   "h",    1  ,  1   ,  2.6569547399e-05  , 1.32234e-05   ,2.1248239171e+09, 0.30   },
    {"He",  "he",   2  ,  4   ,  3.5849373401e-05  , 2.63172e-05   ,1.1671538870e+09, 1.22   },
    {"Li",  "li",   3  ,  7   ,  4.0992133976e-05  , 2.34051e-05   ,8.9266848806e+08, 1.23   },
    {"Be",  "be",   4  ,  9   ,  4.3632829651e-05  , 3.03356e-05   ,7.8788802914e+08, 0.89   },
    {"B",   "b",    5  ,  11  ,  4.5906118608e-05  , 3.54894e-05   ,7.1178709563e+08, 0.88   },
    {"C",   "c",    6  ,  12  ,  4.6940079496e-05  , 3.76762e-05   ,6.8077502929e+08, 0.77   },
    {"N",   "n",    7  ,  14  ,  4.8847128967e-05  , 4.15204e-05   ,6.2865615725e+08, 0.70   },
    {"O",   "o",    8  ,  16  ,  5.0580178957e-05  , 4.48457e-05   ,5.8631436655e+08, 0.66   },
    {"F",   "f",    9  ,  19  ,  5.2927138943e-05  , 4.91529e-05   ,5.3546911034e+08, 0.58   },
    {"Ne",  "ne",  10  ,  20  ,  5.3654104231e-05  , 5.04494e-05   ,5.2105715255e+08, 1.60   },
    {"Na",  "na",  11  ,  23  ,  5.5699159416e-05  , 5.40173e-05   ,4.8349721509e+08, 1.66   },
    {"Mg",  "mg",  12  ,  24  ,  5.6341070732e-05  , 5.51157e-05   ,4.7254270882e+08, 1.36   },
    {"Al",  "al",  13  ,  27  ,  5.8165765928e-05  , 5.81891e-05   ,4.4335984491e+08, 1.25   },
    {"Si",  "si",  14  ,  28  ,  5.8743802504e-05  , 5.91490e-05   ,4.3467748823e+08, 1.17   },
    {"P",   "p",   15  ,  31  ,  6.0399312923e-05  , 6.18655e-05   ,4.1117553148e+08, 1.10   },
    {"S",   "s",   16  ,  32  ,  6.0927308666e-05  , 6.27224e-05   ,4.0407992047e+08, 1.04   },
    {"Cl",  "cl",  17  ,  35  ,  6.2448101115e-05  , 6.51676e-05   ,3.8463852873e+08, 0.99   },
    {"Ar",  "ar",  18  ,  40  ,  6.4800211825e-05  , 6.88887e-05   ,3.5722217300e+08, 1.91   },
    {"K",   "k",   19  ,  39  ,  6.4346167051e-05  , 6.81757e-05   ,3.6228128110e+08, 2.03   },
    {"Ca",  "ca",  20  ,  40  ,  6.4800211825e-05  , 6.88887e-05   ,3.5722217300e+08, 1.74   },
    {"Sc",  "sc",  21  ,  45  ,  6.6963627201e-05  , 7.22548e-05   ,3.3451324570e+08, 1.44   },
    {"Ti",  "ti",  22  ,  48  ,  6.8185577480e-05  , 7.41350e-05   ,3.2263108827e+08, 1.32   },
    {"V",   "v",   23  ,  51  ,  6.9357616830e-05  , 7.59254e-05   ,3.1181925878e+08, 1.22   },
    {"Cr",  "cr",  24  ,  52  ,  6.9738057221e-05  , 7.65040e-05   ,3.0842641793e+08, 1.19   },
    {"Mn",  "mn",  25  ,  55  ,  7.0850896638e-05  , 7.81897e-05   ,2.9881373610e+08, 1.17   },
    {"Fe",  "fe",  26  ,  56  ,  7.1212829817e-05  , 7.87358e-05   ,2.9578406371e+08, 1.165  },
    {"Co",  "co",  27  ,  59  ,  7.2273420879e-05  , 8.03303e-05   ,2.8716667270e+08, 1.16   },
    {"Ni",  "ni",  28  ,  58  ,  7.1923970253e-05  , 7.98058e-05   ,2.8996391416e+08, 1.15   },
    {"Cu",  "cu",  29  ,  63  ,  7.3633018675e-05  , 8.23625e-05   ,2.7665979354e+08, 1.17   },
    {"Zn",  "zn",  30  ,  64  ,  7.3963875193e-05  , 8.28551e-05   ,2.7419021043e+08, 1.25   },
    {"Ga",  "ga",  31  ,  69  ,  7.5568424848e-05  , 8.52341e-05   ,2.6267002737e+08, 1.25   },
    {"Ge",  "ge",  32  ,  74  ,  7.7097216161e-05  , 8.74862e-05   ,2.5235613399e+08, 1.22   },
    {"As",  "as",  33  ,  75  ,  7.7394645153e-05  , 8.79228e-05   ,2.5042024280e+08, 1.21   },
    {"Se",  "se",  34  ,  80  ,  7.8843427408e-05  , 9.00427e-05   ,2.4130163719e+08, 1.17   },
    {"Br",  "br",  35  ,  79  ,  7.8558604038e-05  , 8.96268e-05   ,2.4305454351e+08, 1.14   },
    {"Kr",  "kr",  36  ,  84  ,  7.9959560033e-05  , 9.16684e-05   ,2.3461213272e+08, 1.98   },
    {"Rb",  "rb",  37  ,  85  ,  8.0233033713e-05  , 9.20658e-05   ,2.3301551109e+08, 2.22   },
    {"Sr",  "sr",  38  ,  88  ,  8.1040799081e-05  , 9.32375e-05   ,2.2839354730e+08, 1.92   },
    {"Y",   "y",   39  ,  89  ,  8.1305968993e-05  , 9.36215e-05   ,2.2690621893e+08, 1.62   },
    {"Zr",  "zr",  40  ,  90  ,  8.1569159980e-05  , 9.40022e-05   ,2.2544431039e+08, 1.45   },
    {"Nb",  "nb",  41  ,  93  ,  8.2347219223e-05  , 9.51261e-05   ,2.2120420724e+08, 1.34   },
    {"Mo",  "mo",  42  ,  98  ,  8.3607614434e-05  , 9.69412e-05   ,2.1458511597e+08, 1.29   },
    {"Tc",  "tc",  43  ,  98  ,  8.3607614434e-05  , 9.69412e-05   ,2.1458511597e+08, 1.27   },
    {"Ru",  "ru",  44  , 102  ,  8.4585397905e-05  , 9.83448e-05   ,2.0965270287e+08, 1.24   },
    {"Rh",  "rh",  45  , 103  ,  8.4825835954e-05  , 9.86893e-05   ,2.0846586999e+08, 1.25   },
    {"Pd",  "pd",  46  , 106  ,  8.5537941156e-05  , 9.97084e-05   ,2.0500935221e+08, 1.28   },
    {"Ag",  "ag",  47  , 107  ,  8.5772320442e-05  , 1.00043e-04   ,2.0389047621e+08, 1.34   },
    {"Cd",  "cd",  48  , 114  ,  8.7373430179e-05  , 1.02327e-04   ,1.9648639618e+08, 1.41   },
    {"In",  "in",  49  , 115  ,  8.7596760865e-05  , 1.02644e-04   ,1.9548577691e+08, 1.50   },
    {"Sn",  "sn",  50  , 120  ,  8.8694413774e-05  , 1.04204e-04   ,1.9067718154e+08, 1.40   },
    {"Sb",  "sb",  51  , 121  ,  8.8910267995e-05  , 1.04510e-04   ,1.8975246242e+08, 1.41   },
    {"Te",  "te",  52  , 130  ,  9.0801452955e-05  , 1.07185e-04   ,1.8193056289e+08, 1.37   },
    {"I",   "i",   53  , 127  ,  9.0181040290e-05  , 1.06309e-04   ,1.8444240538e+08, 1.33   },
    {"Xe",  "xe",  54  , 132  ,  9.1209776425e-05  , 1.07762e-04   ,1.8030529331e+08, 2.09   },
    {"Cs",  "cs",  55  , 133  ,  9.1412392742e-05  , 1.08047e-04   ,1.7950688281e+08, 2.35   },
    {"Ba",  "ba",  56  , 138  ,  9.2410525664e-05  , 1.09453e-04   ,1.7565009043e+08, 1.98   },
    {"La",  "la",  57  , 139  ,  9.2607247118e-05  , 1.09730e-04   ,1.7490463170e+08, 1.69   },
    {"Ce",  "ce",  58  , 140  ,  9.2803027311e-05  , 1.10006e-04   ,1.7416744147e+08, 1.65   },
    {"Pr",  "pr",  59  , 141  ,  9.2997877424e-05  , 1.10279e-04   ,1.7343837120e+08, 1.65   },
    {"Nd",  "nd",  60  , 144  ,  9.3576955934e-05  , 1.11093e-04   ,1.7129844956e+08, 1.64   },
    {"Pm",  "pm",  61  , 145  ,  9.3768193375e-05  , 1.11361e-04   ,1.7060044589e+08, 1.65   },
    {"Sm",  "sm",  62  , 152  ,  9.5082839751e-05  , 1.13204e-04   ,1.6591550422e+08, 1.66   },
    {"Eu",  "eu",  63  , 153  ,  9.5267329183e-05  , 1.13462e-04   ,1.6527352089e+08, 1.65   },
    {"Gd",  "gd",  64  , 158  ,  9.6177915369e-05  , 1.14735e-04   ,1.6215880671e+08, 1.61   },
    {"Tb",  "tb",  65  , 159  ,  9.6357719009e-05  , 1.14986e-04   ,1.6155419421e+08, 1.59   },
    {"Dy",  "dy",  66  , 162  ,  9.6892647152e-05  , 1.15733e-04   ,1.5977529080e+08, 1.59   },
    {"Ho",  "ho",  67  , 162  ,  9.6892647152e-05  , 1.15733e-04   ,1.5977529080e+08, 1.58   },
    {"Er",  "er",  68  , 168  ,  9.7943009317e-05  , 1.17198e-04   ,1.5636673634e+08, 1.57   },
    {"Tm",  "tm",  69  , 169  ,  9.8115626740e-05  , 1.17438e-04   ,1.5581702004e+08, 1.56   },
    {"Yb",  "yb",  70  , 174  ,  9.8968651305e-05  , 1.18625e-04   ,1.5314257850e+08, 1.56   },
    {"Lu",  "lu",  71  , 175  ,  9.9137288835e-05  , 1.18859e-04   ,1.5262201512e+08, 1.56   },
    {"Hf",  "hf",  72  , 180  ,  9.9970978172e-05  , 1.20018e-04   ,1.5008710340e+08, 1.44   },
    {"Ta",  "ta",  73  , 181  ,  1.0013585755e-04  , 1.20246e-04   ,1.4959325643e+08, 1.34   },
    {"W",   "w",   74  , 184  ,  1.0062688070e-04  , 1.20928e-04   ,1.4813689532e+08, 1.30   },
    {"Re",  "re",  75  , 187  ,  1.0111259523e-04  , 1.21601e-04   ,1.4671710337e+08, 1.28   },
    {"Os",  "os",  76  , 192  ,  1.0191070333e-04  , 1.22706e-04   ,1.4442808782e+08, 1.26   },
    {"Ir",  "ir",  77  , 193  ,  1.0206865731e-04  , 1.22925e-04   ,1.4398142103e+08, 1.26   },
    {"Pt",  "pt",  78  , 195  ,  1.0238293593e-04  , 1.23360e-04   ,1.4309883584e+08, 1.29   },
    {"Au",  "au",  79  , 197  ,  1.0269507292e-04  , 1.23792e-04   ,1.4223027307e+08, 1.34   },
    {"Hg",  "hg",  80  , 202  ,  1.0346628039e-04  , 1.24857e-04   ,1.4011788914e+08, 1.44   },
    {"Tl",  "tl",  81  , 205  ,  1.0392291259e-04  , 1.25488e-04   ,1.3888925203e+08, 1.55   },
    {"Pb",  "pb",  82  , 208  ,  1.0437511130e-04  , 1.26112e-04   ,1.3768840081e+08, 1.54   },
    {"Bi",  "bi",  83  , 209  ,  1.0452487744e-04  , 1.26318e-04   ,1.3729411599e+08, 1.52   },
    {"Po",  "po",  84  , 209  ,  1.0452487744e-04  , 1.26318e-04   ,1.3729411599e+08, 1.53   },
    {"At",  "at",  85  , 210  ,  1.0467416660e-04  , 1.26524e-04   ,1.3690277000e+08, 1.50   },
    {"Rn",  "rn",  86  , 222  ,  1.0642976299e-04  , 1.28942e-04   ,1.3242350205e+08, 2.20   },
    {"Fr",  "fr",  87  , 223  ,  1.0657317899e-04  , 1.29139e-04   ,1.3206733609e+08, 3.24   },
    {"Ra",  "ra",  88  , 226  ,  1.0700087100e-04  , 1.29727e-04   ,1.3101367628e+08, 2.68   },
    {"Ac",  "ac",  89  , 227  ,  1.0714259349e-04  , 1.29922e-04   ,1.3066730974e+08, 2.25   },
    {"Th",  "th",  90  , 232  ,  1.0784503195e-04  , 1.30887e-04   ,1.2897067480e+08, 2.16   },
    {"Pa",  "pa",  91  , 231  ,  1.0770535752e-04  , 1.30695e-04   ,1.2930539512e+08, 1.93   },
    {"U",   "u",   92  , 238  ,  1.0867476102e-04  , 1.32026e-04   ,1.2700881714e+08, 3.00   },
    {"Np",  "np",  93  , 237  ,  1.0853744903e-04  , 1.31838e-04   ,1.2733038109e+08, 1.57   },
    {"Pu",  "pu",  94  , 244  ,  1.0949065967e-04  , 1.33145e-04   ,1.2512299012e+08, 1.81   },
    {"Am",  "am",  95  , 243  ,  1.0935561268e-04  , 1.32960e-04   ,1.2543221826e+08, 2.21   },
    {"Cm",  "cm",  96  , 247  ,  1.0989359973e-04  , 1.33697e-04   ,1.2420711085e+08, 1.43   },
    {"Bk",  "bk",  97  , 247  ,  1.0989359973e-04  , 1.33697e-04   ,1.2420711085e+08, 1.42   },
    {"Cf",  "cf",  98  , 251  ,  1.1042580946e-04  , 1.34426e-04   ,1.2301273547e+08, 1.40   },
    {"Es",  "es",  99  , 252  ,  1.1055797721e-04  , 1.34607e-04   ,1.2271879740e+08, 1.39   },
    {"Fm",  "fm",  100 , 257  ,  1.1121362374e-04  , 1.35504e-04   ,1.2127611477e+08, 1.38   },
    {"Md",  "md",  101 , 258  ,  1.1134373034e-04  , 1.35682e-04   ,1.2099285491e+08, 1.37   },
    {"No",  "no",  102 , 259  ,  1.1147350119e-04  , 1.35859e-04   ,1.2071131346e+08, 1.36   },
    {"Lr",  "lr",  103 , 262  ,  1.1186082063e-04  , 1.36389e-04   ,1.1987683191e+08, 1.34   },
    {"Db",  "db",  104 , 261  ,  1.1173204420e-04  , 1.36213e-04   ,1.2015331850e+08, 1.40   },
    {"Jl",  "jl",  105 , 262  ,  1.1186082063e-04  , 1.36389e-04   ,1.1987683191e+08, 1.40   },
    {"Rf",  "rf",  106 , 263  ,  1.1198926979e-04  , 1.36565e-04   ,1.1960199758e+08, 1.40   },
    {"Bh",  "bh",  107 , 262  ,  1.1186082063e-04  , 1.36389e-04   ,1.1987683191e+08, 1.40   },
    {"Hn",  "hn",  108 , 265  ,  1.1224519460e-04  , 1.36914e-04   ,1.1905722195e+08, 1.40   },
    {"Mt",  "mt",  109 , 266  ,  1.1237267433e-04  , 1.37088e-04   ,1.1878724932e+08, 1.40   } };

const AtomicData& get_atomic_data(unsigned int atomic_number) {
    if (atomic_number >= NUMBER_OF_ATOMS_IN_TABLE) throw "I am not an alchemist";
    return atomic_data[atomic_number];
}


unsigned int symbol_to_atomic_number(const std::string& symbol) {
    std::string tlow = madness::lowercase(symbol);
    for (unsigned int i=0; i<NUMBER_OF_ATOMS_IN_TABLE; i++) {
        if (tlow.compare(atomic_data[i].symbol_lowercase) == 0) return i;
    }
    std::string msg = "unknown atom -- " + symbol;
    throw msg;
}


/// Returns radius for smoothing nuclear potential with energy precision eprec
static double smoothing_parameter(double Z, double eprec) {
    // The min is since asymptotic form not so good at low acc.
    // The 2 is from two electrons in 1s closed shell.
    if (Z == 0.0) return 1.0;
    double Z5 = Z*Z*Z*Z*Z;
    double c = pow(std::min(1e-3,eprec)/2.0/0.00435/Z5,1.0/3.0);
    return c;
}


/// Regularized 1/r potential.

/// Invoke as \c u(r/c)/c where \c c is the radius of the
/// smoothed volume.
static double smoothed_potential(double r) {
    const double THREE_SQRTPI = 5.31736155271654808184;
    double r2 = r*r, pot;
    if (r > 6.5){
        pot = 1.0/r;
    } else if (r > 1e-8){
        pot = erf(r)/r + (exp(-r2) + 16.0*exp(-4.0*r2))/(THREE_SQRTPI);
    } else{
        pot = (2.0 + 17.0/3.0)/sqrt(PI);
    }

    return pot;
}

/// Derivative of the regularized 1/r potential

/// dV/dx = (x/r) * du(r/c)/(c*c)
static double dsmoothed_potential(double r)
{
    const double THREE_SQRTPI = 5.31736155271654808184;
    const double SQRTPI = 1.77245385090551602728;
    double value;
    double r2 = r*r;

    if (r > 6.5)
        value = -1.0/r2;
    else {
        if (r > 1e-8){
            value =  (-2.0*r*exp(-r2)-128.0*r*exp(-4.0*r2))/(THREE_SQRTPI);
            value += 2.0*exp(-r2)/(r*SQRTPI) - erf(r)/r2;
        }else if (r > 0.1) {
            value =  (-2.0*exp(-r2)-128.0*exp(-4.0*r2))/(THREE_SQRTPI);
            value += (-4./3.+(4./5.+(-2./7.+(2./27.-1./66.*r2)*r2)*r2)*r2)/SQRTPI;
            value *= r;
        }else if (r != 0.0){
            value =  (-2.0*exp(-r2)-128.0*exp(-4.0*r2))/(THREE_SQRTPI);
            value += (-4./3.+(4./5.+(-2./7.+(2./27.-1./66.*r2)*r2)*r2)*r2)/SQRTPI;
            value *= r;
        }else
            value = 0.0;
    }
    return value;
}

/// Charge density corresponding to smoothed 1/r potential

/// Invoke as \c rho(r/c)/c^3 where \c c is the radius of the
/// smoothed volume.
static double smoothed_density(double r)
{
    const double RPITO1P5 = 0.1795871221251665617; // 1.0/Pi^1.5
    double rsquared = r*r;
    double tmp1 = exp(-rsquared);
    double tmp2 = tmp1*tmp1;
    double tmp4 = tmp2*tmp2;
    return ((-3.0/2.0+(1.0/3.0)*rsquared)*tmp1+(-32.0+(256.0/3.0)*rsquared)*tmp4)*RPITO1P5;
}

std::ostream& operator<<(std::ostream& s, const Atom& atom) {
    s << "Atom([" << atom.x << ", " << atom.y << ", " << atom.z << "], " << atom.q << "," << atom.atomic_number << ")";
    return s;
}

/// Read coordinates from a file

/// Scans the file for the first geometry block in the format
/// \code
///    geometry
///       tag x y z
///       ...
///    end
/// \endcode
/// The charge \c q is inferred from the tag which is
/// assumed to be the standard symbol for an element.
/// Same as the simplest NWChem format.  For ghost
/// atoms (\c bq ) the  charge is read as a fifth field
/// on the line.
///
/// This code is just for the examples ... don't trust it!
MolecularEntity::MolecularEntity(const std::string& filename, bool fractional = false) {
    read_file(filename, fractional);
}

void MolecularEntity::read_file(const std::string& filename, bool fractional = false) {
    atoms.clear(); rcut.clear(); rsqasymptotic.clear();
    std::ifstream f(filename.c_str());
    madness::position_stream(f, "geometry");
    double scale = 1.0;

    std::string s;
    while (std::getline(f,s))
    {
        std::istringstream ss(s);
        std::string tag;
        ss >> tag;
        if (tag == "end")
        {
            goto finished;
        }
        else if (tag == "units")
        {
            if (natom()) throw "MolecularEntity: read_file: presently units must be the first line of the geometry block";
            ss >> tag;
            if (tag == "a.u." || tag == "au" || tag == "atomic")
            {
                std::cout << "\nAtomic units being used " << scale << "\n\n";
                scale = 1.0;
            }
            else if (tag == "angstrom" || tag == "angs")
            {
                scale = 1e-10 / madness::constants::atomic_unit_of_length;
                std::cout << "\nAngstrom being used " << scale << "\n\n";
            }
            else
            {
                throw "MolecularEntity: read_file: unknown units requested";
            }
        }
        else
        {
            Tensor<double> factor = FunctionDefaults<3>::get_cell_width();
            double xx, yy, zz;
            ss >> xx >> yy >> zz;
            if (fractional)
            {
              // If using fractional coordinates, the restrict x, y, and z to be between 0.0 and 1.0
              MADNESS_ASSERT(xx <= 1.0);
              MADNESS_ASSERT(yy <= 1.0);
              MADNESS_ASSERT(zz <= 1.0);
              MADNESS_ASSERT(xx >= 0.0);
              MADNESS_ASSERT(yy >= 0.0);
              MADNESS_ASSERT(zz >= 0.0);
              xx *= factor[0]; yy *= factor[1]; zz *= factor[2];
            }
            else
            {
              xx *= scale; yy *= scale; zz *= scale;
            }
            int atn = symbol_to_atomic_number(tag); // Charge of ghost atom
            double qq = atn;
            if (atn == 0) ss >> qq;
            add_atom(xx,yy,zz,qq,atn);
        }
    }
    throw "No end to the geometry in the input file";
 finished:
    ;
}

void MolecularEntity::add_atom(double x, double y, double z, int atomic_number, double q) {
    atoms.push_back(Atom(x,y,z,atomic_number,q));
    double c = smoothing_parameter(q, 1e-5); // This is error per atom
    rsqasymptotic.push_back(36.0*c*c);
    rcut.push_back(1.0/c);
}

void MolecularEntity::set_atom_coords(unsigned int i, double x, double y, double z) {
    if (i>=atoms.size()) throw "trying to set coords of invalid atom";
    atoms[i].x = x;
    atoms[i].y = y;
    atoms[i].z = z;
}

const Atom& MolecularEntity::get_atom(unsigned int i) const {
    if (i>=atoms.size()) throw "trying to get coords of invalid atom";
    return atoms[i];
}

void MolecularEntity::print() const {
    std::cout.flush();
    printf(" geometry\n");
    for (int i=0; i<natom(); i++) {
        printf("   %-2s  %20.8f %20.8f %20.8f", atomic_data[atoms[i].atomic_number].symbol,
               atoms[i].x, atoms[i].y, atoms[i].z);
        if (atoms[i].atomic_number == 0) printf("     %20.8f", atoms[i].q);
        printf("\n");
    }
    printf(" end\n");
}

double MolecularEntity::inter_atomic_distance(unsigned int i,unsigned int j) const {
    if (i>=atoms.size()) throw "trying to compute distance with invalid atom";
    if (j>=atoms.size()) throw "trying to compute distance with invalid atom";
    return distance(atoms[i].x, atoms[i].y, atoms[i].z,
                    atoms[j].x, atoms[j].y, atoms[j].z);
}

double MolecularEntity::nuclear_repulsion_energy() const {
    double sum = 0.0;
    for (unsigned int i=0; i<atoms.size(); i++) {
        for (unsigned int j=i+1; j<atoms.size(); j++) {
            sum += atoms[i].atomic_number * atoms[j].atomic_number / inter_atomic_distance(i,j);
        }
    }
    return sum;
}

double MolecularEntity::smallest_length_scale() const {
    double rcmax = 0.0;
    for (unsigned int i=0; i<atoms.size(); i++) {
        rcmax = std::max(rcmax,rcut[i]);
    }
    return 1.0/rcmax;
}


/// Moves the center of nuclear charge to the origin
void MolecularEntity::center() {
    double xx=0.0, yy=0.0, zz=0.0, qq=0.0;
    for (unsigned int i=0; i<atoms.size(); i++) {
        xx += atoms[i].x*atoms[i].q;
        yy += atoms[i].y*atoms[i].q;
        zz += atoms[i].z*atoms[i].q;
        qq += atoms[i].q;
    }
    xx /= qq;
    yy /= qq;
    zz /= qq;
    for (unsigned int i=0; i<atoms.size(); i++) {
        atoms[i].x -= xx;
        atoms[i].y -= yy;
        atoms[i].z -= zz;
    }
}

/// Returns the half width of the bounding cube

/// The MolecularEntity will be contained in the cube [-L,+L].
double MolecularEntity::bounding_cube() const {
    double L = 0.0;
    for (unsigned int i=0; i<atoms.size(); i++) {
        L = std::max(L, fabs(atoms[i].x));
        L = std::max(L, fabs(atoms[i].y));
        L = std::max(L, fabs(atoms[i].z));
    }
    return L;
}

double MolecularEntity::total_nuclear_charge() const {
    double sum = 0.0;
    for (unsigned int i=0; i<atoms.size(); i++) {
        sum += atoms[i].q;
    }
    return sum;
}

double MolecularEntity::nuclear_attraction_potential(double x, double y, double z) const {
    // This is very inefficient since it scales as O(ngrid*natom)
    // ... we can easily make an O(natom) version using
    // the integral operator and sparse projection of an effective
    // density ... its potential can be evaluated at the same
    // time as the electronic Coulomb potential so it will be
    // essentially free.

    double sum = 0.0;
    for (unsigned int i=0; i<atoms.size(); i++) {
        double r = distance(atoms[i].x, atoms[i].y, atoms[i].z, x, y, z);
        //sum -= atoms[i].q/(r+1e-8);
        sum -= atoms[i].q * smoothed_potential(r*rcut[i])*rcut[i];
    }
    return sum;
}

double MolecularEntity::nuclear_charge_density(double x, double y, double z) const {
  // Only one atom will contribute due to the short range of the nuclear
  // charge density
  for (unsigned int i=0; i<atoms.size(); i++) {
      double big = rsqasymptotic[i];
      double xx = atoms[i].x - x;
      double rsq = xx*xx;
      if (rsq <  big) {
          double yy = atoms[i].y - y;
          rsq += yy*yy;
          if (rsq < big) {
              double zz = atoms[i].z - z;
              rsq += zz*zz;
              if (rsq < big) {
                  double r = sqrt(rsq);
                  return atoms[i].atomic_number * smoothed_density(r*rcut[i])*rcut[i]*rcut[i]*rcut[i];
              }
          }
      }
  }
  return 0.0;

//  double sum = 0.0;
//  for (unsigned int i = 0; i < atoms.size(); i++)
//  {
//    double r = distance(atoms[i].x, atoms[i].y, atoms[i].z, x, y, z);
//    double e1 = 50.0;
//    double coeff = pow(e1/PI, 1.5);
//    sum -= atoms[i].atomic_number * coeff * exp(-e1 * r * r);
//  }
//  return sum;
}


