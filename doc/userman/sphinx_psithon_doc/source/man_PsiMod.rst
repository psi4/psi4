Help on module PsiMod ::

    NAME
        PsiMod
    
    FILE
        (built-in)
    
    CLASSES
        Boost.Python.enum(__builtin__.int)
            DiagonalizeOrder
            PsiReturnType
        Boost.Python.instance(__builtin__.object)
            Arguments
            BasisSet
            BasisSetParser
                Gaussian94BasisSetParser
            CdSalcList
            Checkpoint
            DFChargeFitter
            Environment
            ExternalPotential
            FittingMetric
            Functional
            GridProp
            IO
            IOManager
            IntVector
            Matrix
            MatrixFactory
            MintsHelper
            MoldenWriter
            Molecule
            MultipoleSymmetry
            NBOWriter
            OEProp
            PetiteList
            PointGroup
            Process
            PseudoTrial
            SOBasisSet
            SuperFunctional
            SymmetryOperation
            Vector
            Vector3
            Wavefunction
                HF
                    RHF(HF, Wavefunction)
            matrix_vector
        
        class Arguments(Boost.Python.instance)
         |  Method resolution order:
         |      Arguments
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __getitem__(...)
         |      __getitem__( (Arguments)arg1, (int)arg2) -> str :
         |          docstring
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 40
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class BasisSet(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      BasisSet
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  max_am(...)
         |      max_am( (BasisSet)arg1) -> int :
         |          docstring
         |  
         |  nao(...)
         |      nao( (BasisSet)arg1) -> int :
         |          docstring
         |  
         |  nbf(...)
         |      nbf( (BasisSet)arg1) -> int :
         |          docstring
         |  
         |  nprimitive(...)
         |      nprimitive( (BasisSet)arg1) -> int :
         |          docstring
         |  
         |  nshell(...)
         |      nshell( (BasisSet)arg1) -> int :
         |          docstring
         |  
         |  print_detail_out(...)
         |      print_detail_out( (BasisSet)arg1) -> None :
         |          docstring
         |  
         |  print_out(...)
         |      print_out( (BasisSet)arg1) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Static methods defined here:
         |  
         |  construct(...)
         |      construct( (BasisSetParser)arg1, (Molecule)arg2, (str)arg3) -> BasisSet :
         |          docstring
         |  
         |  make_filename(...)
         |      make_filename( (str)arg1) -> str :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __init__ = <built-in function __init__>
         |      Raises an exception
         |      This class cannot be instantiated from Python
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class BasisSetParser(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      BasisSetParser
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __init__ = <built-in function __init__>
         |      Raises an exception
         |      This class cannot be instantiated from Python
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class CdSalcList(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      CdSalcList
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  matrix(...)
         |      matrix( (CdSalcList)arg1) -> Matrix :
         |          docstring
         |  
         |  print_out(...)
         |      print_out( (CdSalcList)arg1) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __init__ = <built-in function __init__>
         |      Raises an exception
         |      This class cannot be instantiated from Python
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class Checkpoint(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      Checkpoint
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1, (IO)arg2, (int)arg3) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  ----------------------------------------------------------------------
         |  Static methods defined here:
         |  
         |  shared_object(...)
         |      shared_object() -> Checkpoint :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors defined here:
         |  
         |  disp
         |      docstring
         |  
         |  e_t
         |      docstring
         |  
         |  eccsd
         |      docstring
         |  
         |  ecorr
         |      docstring
         |  
         |  efzc
         |      docstring
         |  
         |  emp2
         |      docstring
         |  
         |  enuc
         |      docstring
         |  
         |  eref
         |      docstring
         |  
         |  escf
         |      docstring
         |  
         |  etot
         |      docstring
         |  
         |  label
         |      docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class DFChargeFitter(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      DFChargeFitter
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  d(...)
         |      d( (DFChargeFitter)arg1) -> Vector :
         |          docstring
         |  
         |  fit(...)
         |      fit( (DFChargeFitter)arg1) -> Vector :
         |          docstring
         |  
         |  setAuxiliary(...)
         |      setAuxiliary( (DFChargeFitter)arg1, (BasisSet)arg2) -> None :
         |          docstring
         |  
         |  setD(...)
         |      setD( (DFChargeFitter)arg1, (Matrix)arg2) -> None :
         |          docstring
         |  
         |  setPrimary(...)
         |      setPrimary( (DFChargeFitter)arg1, (BasisSet)arg2) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class DiagonalizeOrder(Boost.Python.enum)
         |  docstring
         |  
         |  Method resolution order:
         |      DiagonalizeOrder
         |      Boost.Python.enum
         |      __builtin__.int
         |      __builtin__.object
         |  
         |  Data and other attributes defined here:
         |  
         |  Ascending = PsiMod.DiagonalizeOrder.Ascending
         |  
         |  Descending = PsiMod.DiagonalizeOrder.Descending
         |  
         |  names = {'Ascending': PsiMod.DiagonalizeOrder.Ascending, 'Descending':...
         |  
         |  values = {1: PsiMod.DiagonalizeOrder.Ascending, 3: PsiMod.DiagonalizeO...
         |  
         |  ----------------------------------------------------------------------
         |  Methods inherited from Boost.Python.enum:
         |  
         |  __repr__(...)
         |      x.__repr__() <==> repr(x)
         |  
         |  __str__(...)
         |      x.__str__() <==> str(x)
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.enum:
         |  
         |  name
         |  
         |  ----------------------------------------------------------------------
         |  Methods inherited from __builtin__.int:
         |  
         |  __abs__(...)
         |      x.__abs__() <==> abs(x)
         |  
         |  __add__(...)
         |      x.__add__(y) <==> x+y
         |  
         |  __and__(...)
         |      x.__and__(y) <==> x&y
         |  
         |  __cmp__(...)
         |      x.__cmp__(y) <==> cmp(x,y)
         |  
         |  __coerce__(...)
         |      x.__coerce__(y) <==> coerce(x, y)
         |  
         |  __div__(...)
         |      x.__div__(y) <==> x/y
         |  
         |  __divmod__(...)
         |      x.__divmod__(y) <==> divmod(x, y)
         |  
         |  __float__(...)
         |      x.__float__() <==> float(x)
         |  
         |  __floordiv__(...)
         |      x.__floordiv__(y) <==> x//y
         |  
         |  __format__(...)
         |  
         |  __getattribute__(...)
         |      x.__getattribute__('name') <==> x.name
         |  
         |  __getnewargs__(...)
         |  
         |  __hash__(...)
         |      x.__hash__() <==> hash(x)
         |  
         |  __hex__(...)
         |      x.__hex__() <==> hex(x)
         |  
         |  __index__(...)
         |      x[y:z] <==> x[y.__index__():z.__index__()]
         |  
         |  __int__(...)
         |      x.__int__() <==> int(x)
         |  
         |  __invert__(...)
         |      x.__invert__() <==> ~x
         |  
         |  __long__(...)
         |      x.__long__() <==> long(x)
         |  
         |  __lshift__(...)
         |      x.__lshift__(y) <==> x<<y
         |  
         |  __mod__(...)
         |      x.__mod__(y) <==> x%y
         |  
         |  __mul__(...)
         |      x.__mul__(y) <==> x*y
         |  
         |  __neg__(...)
         |      x.__neg__() <==> -x
         |  
         |  __nonzero__(...)
         |      x.__nonzero__() <==> x != 0
         |  
         |  __oct__(...)
         |      x.__oct__() <==> oct(x)
         |  
         |  __or__(...)
         |      x.__or__(y) <==> x|y
         |  
         |  __pos__(...)
         |      x.__pos__() <==> +x
         |  
         |  __pow__(...)
         |      x.__pow__(y[, z]) <==> pow(x, y[, z])
         |  
         |  __radd__(...)
         |      x.__radd__(y) <==> y+x
         |  
         |  __rand__(...)
         |      x.__rand__(y) <==> y&x
         |  
         |  __rdiv__(...)
         |      x.__rdiv__(y) <==> y/x
         |  
         |  __rdivmod__(...)
         |      x.__rdivmod__(y) <==> divmod(y, x)
         |  
         |  __rfloordiv__(...)
         |      x.__rfloordiv__(y) <==> y//x
         |  
         |  __rlshift__(...)
         |      x.__rlshift__(y) <==> y<<x
         |  
         |  __rmod__(...)
         |      x.__rmod__(y) <==> y%x
         |  
         |  __rmul__(...)
         |      x.__rmul__(y) <==> y*x
         |  
         |  __ror__(...)
         |      x.__ror__(y) <==> y|x
         |  
         |  __rpow__(...)
         |      y.__rpow__(x[, z]) <==> pow(x, y[, z])
         |  
         |  __rrshift__(...)
         |      x.__rrshift__(y) <==> y>>x
         |  
         |  __rshift__(...)
         |      x.__rshift__(y) <==> x>>y
         |  
         |  __rsub__(...)
         |      x.__rsub__(y) <==> y-x
         |  
         |  __rtruediv__(...)
         |      x.__rtruediv__(y) <==> y/x
         |  
         |  __rxor__(...)
         |      x.__rxor__(y) <==> y^x
         |  
         |  __sub__(...)
         |      x.__sub__(y) <==> x-y
         |  
         |  __truediv__(...)
         |      x.__truediv__(y) <==> x/y
         |  
         |  __trunc__(...)
         |      Truncating an Integral returns itself.
         |  
         |  __xor__(...)
         |      x.__xor__(y) <==> x^y
         |  
         |  bit_length(...)
         |      int.bit_length() -> int
         |      
         |      Number of bits necessary to represent self in binary.
         |      >>> bin(37)
         |      '0b100101'
         |      >>> (37).bit_length()
         |      6
         |  
         |  conjugate(...)
         |      Returns self, the complex conjugate of any int.
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from __builtin__.int:
         |  
         |  denominator
         |      the denominator of a rational number in lowest terms
         |  
         |  imag
         |      the imaginary part of a complex number
         |  
         |  numerator
         |      the numerator of a rational number in lowest terms
         |  
         |  real
         |      the real part of a complex number
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from __builtin__.int:
         |  
         |  __new__ = <built-in method __new__ of type object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class Environment(Boost.Python.instance)
         |  Method resolution order:
         |      Environment
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __getitem__(...)
         |      __getitem__( (Environment)arg1, (str)arg2) -> str :
         |          docstring
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 352
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class ExternalPotential(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      ExternalPotential
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  addBasis(...)
         |      addBasis( (ExternalPotential)arg1, (BasisSet)arg2, (Vector)arg3) -> None :
         |          docstring
         |  
         |  addCharge(...)
         |      addCharge( (ExternalPotential)arg1, (float)arg2, (float)arg3, (float)arg4, (float)arg5) -> None :
         |          docstring
         |  
         |  clear(...)
         |      clear( (ExternalPotential)arg1) -> None :
         |          docstring
         |  
         |  computePotentialMatrix(...)
         |      computePotentialMatrix( (ExternalPotential)arg1, (BasisSet)arg2) -> Matrix :
         |          docstring
         |  
         |  print_out(...)
         |      print_out( (ExternalPotential)arg1) -> None :
         |          docstring
         |  
         |  setName(...)
         |      setName( (ExternalPotential)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class FittingMetric(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      FittingMetric
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  form_QR_inverse(...)
         |      form_QR_inverse( (FittingMetric)arg1, (float)arg2) -> None :
         |          docstring
         |  
         |  form_cholesky_inverse(...)
         |      form_cholesky_inverse( (FittingMetric)arg1) -> None :
         |          docstring
         |  
         |  form_eig_inverse(...)
         |      form_eig_inverse( (FittingMetric)arg1, (float)arg2) -> None :
         |          docstring
         |  
         |  form_fitting_metric(...)
         |      form_fitting_metric( (FittingMetric)arg1) -> None :
         |          docstring
         |  
         |  form_full_inverse(...)
         |      form_full_inverse( (FittingMetric)arg1) -> None :
         |          docstring
         |  
         |  get_algorithm(...)
         |      get_algorithm( (FittingMetric)arg1) -> str :
         |          docstring
         |  
         |  get_metric(...)
         |      get_metric( (FittingMetric)arg1) -> Matrix :
         |          docstring
         |  
         |  get_pivots(...)
         |      get_pivots( (FittingMetric)arg1) -> IntVector :
         |          docstring
         |  
         |  get_reverse_pivots(...)
         |      get_reverse_pivots( (FittingMetric)arg1) -> IntVector :
         |          docstring
         |  
         |  is_inverted(...)
         |      is_inverted( (FittingMetric)arg1) -> bool :
         |          docstring
         |  
         |  is_poisson(...)
         |      is_poisson( (FittingMetric)arg1) -> bool :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class Functional(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      Functional
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  computeRKSFunctional(...)
         |      computeRKSFunctional( (Functional)arg1, (object)arg2) -> None :
         |          docstring
         |  
         |  computeUKSFunctional(...)
         |      computeUKSFunctional( (Functional)arg1, (object)arg2) -> None :
         |          docstring
         |  
         |  get_citation(...)
         |      get_citation( (Functional)arg1) -> str :
         |          docstring
         |  
         |  get_density_cutoff(...)
         |      get_density_cutoff( (Functional)arg1) -> float :
         |          docstring
         |  
         |  get_deriv(...)
         |      get_deriv( (Functional)arg1) -> int :
         |          docstring
         |  
         |  get_description(...)
         |      get_description( (Functional)arg1) -> str :
         |          docstring
         |  
         |  get_name(...)
         |      get_name( (Functional)arg1) -> str :
         |          docstring
         |  
         |  get_npoints(...)
         |      get_npoints( (Functional)arg1) -> int :
         |          docstring
         |  
         |  get_parameters(...)
         |      get_parameters( (Functional)arg1) -> object :
         |          docstring
         |  
         |  get_parameters_string(...)
         |      get_parameters_string( (Functional)arg1) -> str :
         |          docstring
         |  
         |  is_gga(...)
         |      is_gga( (Functional)arg1) -> bool :
         |          docstring
         |  
         |  is_meta(...)
         |      is_meta( (Functional)arg1) -> bool :
         |          docstring
         |  
         |  set_citation(...)
         |      set_citation( (Functional)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  set_density_cutoff(...)
         |      set_density_cutoff( (Functional)arg1, (float)arg2) -> None :
         |          docstring
         |  
         |  set_deriv(...)
         |      set_deriv( (Functional)arg1, (int)arg2) -> None :
         |          docstring
         |  
         |  set_description(...)
         |      set_description( (Functional)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  set_name(...)
         |      set_name( (Functional)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  set_npoints(...)
         |      set_npoints( (Functional)arg1, (int)arg2) -> None :
         |          docstring
         |  
         |  set_parameter(...)
         |      set_parameter( (Functional)arg1, (str)arg2, (float)arg3) -> None :
         |          docstring
         |  
         |  set_parameters(...)
         |      set_parameters( (Functional)arg1, (object)arg2) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Static methods defined here:
         |  
         |  available_functionals(...)
         |      available_functionals() -> str :
         |          docstring
         |  
         |  available_names(...)
         |      available_names() -> object :
         |          docstring
         |  
         |  create_functional(...)
         |      create_functional( (str)arg1, (int)arg2, (int)arg3) -> Functional :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __init__ = <built-in function __init__>
         |      Raises an exception
         |      This class cannot be instantiated from Python
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class Gaussian94BasisSetParser(BasisSetParser)
         |  docstring
         |  
         |  Method resolution order:
         |      Gaussian94BasisSetParser
         |      BasisSetParser
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class GridProp(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      GridProp
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  add(...)
         |      add( (GridProp)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  add_alpha_mo(...)
         |      add_alpha_mo( (GridProp)arg1, (int)arg2, (int)arg3) -> None :
         |          docstring
         |  
         |  add_basis_fun(...)
         |      add_basis_fun( (GridProp)arg1, (int)arg2, (int)arg3) -> None :
         |          docstring
         |  
         |  add_beta_mo(...)
         |      add_beta_mo( (GridProp)arg1, (int)arg2, (int)arg3) -> None :
         |          docstring
         |  
         |  build_grid_overages(...)
         |      build_grid_overages( (GridProp)arg1, (float)arg2) -> None :
         |          docstring
         |  
         |  compute(...)
         |      compute( (GridProp)arg1) -> None :
         |          docstring
         |  
         |  get_l(...)
         |      get_l( (GridProp)arg1, (int)arg2) -> float :
         |          docstring
         |  
         |  get_n(...)
         |      get_n( (GridProp)arg1, (int)arg2) -> int :
         |          docstring
         |  
         |  get_o(...)
         |      get_o( (GridProp)arg1, (int)arg2) -> float :
         |          docstring
         |  
         |  set_caxis(...)
         |      set_caxis( (GridProp)arg1, (float)arg2, (float)arg3) -> None :
         |          docstring
         |  
         |  set_filename(...)
         |      set_filename( (GridProp)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  set_format(...)
         |      set_format( (GridProp)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  set_l(...)
         |      set_l( (GridProp)arg1, (float)arg2, (float)arg3, (float)arg4) -> None :
         |          docstring
         |  
         |  set_n(...)
         |      set_n( (GridProp)arg1, (int)arg2, (int)arg3, (int)arg4) -> None :
         |          docstring
         |  
         |  set_o(...)
         |      set_o( (GridProp)arg1, (float)arg2, (float)arg3, (float)arg4) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class HF(Wavefunction)
         |  docstring
         |  
         |  Method resolution order:
         |      HF
         |      Wavefunction
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __init__ = <built-in function __init__>
         |      Raises an exception
         |      This class cannot be instantiated from Python
         |  
         |  ----------------------------------------------------------------------
         |  Methods inherited from Wavefunction:
         |  
         |  Ca(...)
         |      Ca( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Cb(...)
         |      Cb( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Da(...)
         |      Da( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Db(...)
         |      Db( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Fa(...)
         |      Fa( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Fb(...)
         |      Fb( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  add_postiteration_callback(...)
         |      add_postiteration_callback( (Wavefunction)arg1, (object)arg2) -> None :
         |          docstring
         |  
         |  add_preiteration_callback(...)
         |      add_preiteration_callback( (Wavefunction)arg1, (object)arg2) -> None :
         |          docstring
         |  
         |  basisset(...)
         |      basisset( (Wavefunction)arg1) -> BasisSet :
         |          docstring
         |  
         |  energy(...)
         |      energy( (Wavefunction)arg1) -> float :
         |          docstring
         |  
         |  epsilon_a(...)
         |      epsilon_a( (Wavefunction)arg1) -> Vector :
         |          docstring
         |  
         |  epsilon_b(...)
         |      epsilon_b( (Wavefunction)arg1) -> Vector :
         |          docstring
         |  
         |  frequencies(...)
         |      frequencies( (Wavefunction)arg1) -> Vector :
         |          docstring
         |  
         |  gradient(...)
         |      gradient( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  nirrep(...)
         |      nirrep( (Wavefunction)arg1) -> int :
         |          docstring
         |  
         |  nmo(...)
         |      nmo( (Wavefunction)arg1) -> int :
         |          docstring
         |  
         |  nso(...)
         |      nso( (Wavefunction)arg1) -> int :
         |          docstring
         |  
         |  sobasisset(...)
         |      sobasisset( (Wavefunction)arg1) -> SOBasisSet :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class IO(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      IO
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  close(...)
         |      close( (IO)arg1, (int)arg2, (int)arg3) -> None :
         |          docstring
         |  
         |  open(...)
         |      open( (IO)arg1, (int)arg2, (int)arg3) -> None :
         |          docstring
         |  
         |  open_check(...)
         |      open_check( (IO)arg1, (int)arg2) -> int :
         |          docstring
         |  
         |  rehash(...)
         |      rehash( (IO)arg1, (int)arg2) -> None :
         |          docstring
         |  
         |  state(...)
         |      state( (IO)arg1) -> int :
         |          docstring
         |  
         |  tocclean(...)
         |      tocclean( (IO)arg1, (int)arg2, (str)arg3) -> None :
         |          docstring
         |  
         |  tocprint(...)
         |      tocprint( (IO)arg1, (int)arg2) -> None :
         |          docstring
         |  
         |  tocwrite(...)
         |      tocwrite( (IO)arg1, (int)arg2) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Static methods defined here:
         |  
         |  change_file_namespace(...)
         |      change_file_namespace( (int)arg1, (str)arg2, (str)arg3) -> None :
         |          docstring
         |  
         |  get_default_namespace(...)
         |      get_default_namespace() -> str :
         |          docstring
         |  
         |  set_default_namespace(...)
         |      set_default_namespace( (str)arg1) -> None :
         |          docstring
         |  
         |  shared_object(...)
         |      shared_object() -> IO
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class IOManager(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      IOManager
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  crashclean(...)
         |      crashclean( (IOManager)arg1) -> None :
         |          docstring
         |  
         |  get_default_path(...)
         |      get_default_path( (IOManager)arg1) -> str :
         |          docstring
         |  
         |  get_file_path(...)
         |      get_file_path( (IOManager)arg1, (int)arg2) -> str :
         |          docstring
         |  
         |  mark_file_for_retention(...)
         |      mark_file_for_retention( (IOManager)arg1, (str)arg2, (bool)arg3) -> None :
         |          docstring
         |  
         |  print_out(...)
         |      print_out( (IOManager)arg1) -> None :
         |          docstring
         |  
         |  psiclean(...)
         |      psiclean( (IOManager)arg1) -> None :
         |          docstring
         |  
         |  set_default_path(...)
         |      set_default_path( (IOManager)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  set_specific_path(...)
         |      set_specific_path( (IOManager)arg1, (int)arg2, (str)arg3) -> None :
         |          docstring
         |  
         |  set_specific_retention(...)
         |      set_specific_retention( (IOManager)arg1, (int)arg2, (bool)arg3) -> None :
         |          docstring
         |  
         |  write_scratch_file(...)
         |      write_scratch_file( (IOManager)arg1, (str)arg2, (str)arg3) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Static methods defined here:
         |  
         |  shared_object(...)
         |      shared_object() -> IOManager :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class IntVector(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      IntVector
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |      
         |      __init__( (object)arg1, (int)arg2) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  dim(...)
         |      dim( (IntVector)arg1, (int)arg2) -> int :
         |          docstring
         |  
         |  get(...)
         |      get( (IntVector)arg1, (int)arg2, (int)arg3) -> int :
         |          docstring
         |  
         |  nirrep(...)
         |      nirrep( (IntVector)arg1) -> int :
         |          docstring
         |  
         |  print_out(...)
         |      print_out( (IntVector)arg1) -> None :
         |          docstring
         |  
         |  set(...)
         |      set( (IntVector)arg1, (int)arg2, (int)arg3, (int)arg4) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class Matrix(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      Matrix
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __getitem__(...)
         |      __getitem__( (Matrix)arg1, (tuple)arg2) -> float :
         |          docstring
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |      
         |      __init__( (object)arg1, (int)arg2, (int)arg3) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  __setitem__(...)
         |      __setitem__( (Matrix)arg1, (tuple)arg2, (float)arg3) -> None :
         |          docstring
         |  
         |  accumulate_product(...)
         |      accumulate_product( (Matrix)arg1, (Matrix)arg2, (Matrix)arg3) -> None :
         |          docstring
         |  
         |  add(...)
         |      add( (Matrix)arg1, (Matrix)arg2) -> None :
         |          docstring
         |  
         |  back_transform(...)
         |      back_transform( (Matrix)arg1, (Matrix)arg2, (Matrix)arg3) -> None :
         |          docstring
         |  
         |  cholesky_factorize(...)
         |      cholesky_factorize( (Matrix)arg1) -> None :
         |          docstring
         |  
         |  cols(...)
         |      cols( (Matrix)arg1, (int)arg2) -> int :
         |          docstring
         |  
         |  copy_lower_to_upper(...)
         |      copy_lower_to_upper( (Matrix)arg1) -> None :
         |          docstring
         |  
         |  copy_upper_to_lower(...)
         |      copy_upper_to_lower( (Matrix)arg1) -> None :
         |          docstring
         |  
         |  diagonalize(...)
         |      diagonalize( (Matrix)arg1, (Matrix)arg2, (Vector)arg3, (DiagonalizeOrder)arg4) -> None :
         |          docstring
         |  
         |  exp(...)
         |      exp( (Matrix)arg1) -> None :
         |          docstring
         |  
         |  gemm(...)
         |      gemm( (Matrix)arg1, (bool)arg2, (bool)arg3, (float)arg4, (Matrix)arg5, (Matrix)arg6, (float)arg7) -> None :
         |          docstring
         |  
         |  get(...)
         |      get( (Matrix)arg1, (int)arg2, (int)arg3 [, (int)arg4]) -> float :
         |          docstring
         |  
         |  identity(...)
         |      identity( (Matrix)arg1) -> None :
         |          docstring
         |  
         |  invert(...)
         |      invert( (Matrix)arg1) -> None :
         |          docstring
         |  
         |  load(...)
         |      load( (Matrix)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  name(...)
         |      name( (Matrix)arg1) -> str :
         |          docstring
         |  
         |  nirrep(...)
         |      nirrep( (Matrix)arg1) -> int :
         |          docstring
         |  
         |  partial_cholesky_factorize(...)
         |      partial_cholesky_factorize( (Matrix)arg1, (float)arg2, (bool)arg3) -> Matrix :
         |          docstring
         |  
         |  power(...)
         |      power( (Matrix)arg1, (float)arg2, (float)arg3) -> None :
         |          docstring
         |  
         |  print_out(...)
         |      print_out( (Matrix)arg1) -> None :
         |          docstring
         |  
         |  project_out(...)
         |      project_out( (Matrix)arg1, (Matrix)arg2) -> None :
         |          docstring
         |  
         |  remove_symmetry(...)
         |      remove_symmetry( (Matrix)arg1, (Matrix)arg2, (Matrix)arg3) -> None :
         |          docstring
         |  
         |  rms(...)
         |      rms( (Matrix)arg1) -> float :
         |          docstring
         |  
         |  rows(...)
         |      rows( (Matrix)arg1, (int)arg2) -> int :
         |          docstring
         |  
         |  save(...)
         |      save( (Matrix)arg1, (str)arg2, (bool)arg3, (bool)arg4, (bool)arg5) -> None :
         |          docstring
         |  
         |  scale(...)
         |      scale( (Matrix)arg1, (float)arg2) -> None :
         |          docstring
         |  
         |  scale_column(...)
         |      scale_column( (Matrix)arg1, (int)arg2, (int)arg3, (float)arg4) -> None :
         |          docstring
         |  
         |  scale_row(...)
         |      scale_row( (Matrix)arg1, (int)arg2, (int)arg3, (float)arg4) -> None :
         |          docstring
         |  
         |  set(...)
         |      set( (Matrix)arg1, (int)arg2, (int)arg3, (float)arg4) -> None :
         |          docstring
         |      
         |      set( (Matrix)arg1, (int)arg2, (int)arg3, (int)arg4, (float)arg5) -> None :
         |          docstring
         |      
         |      set( (Matrix)arg1, (list)arg2) -> None :
         |          docstring
         |  
         |  set_name(...)
         |      set_name( (Matrix)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  subtract(...)
         |      subtract( (Matrix)arg1, (Matrix)arg2) -> None :
         |          docstring
         |  
         |  sum_of_squares(...)
         |      sum_of_squares( (Matrix)arg1) -> float :
         |          docstring
         |  
         |  symmetry(...)
         |      symmetry( (Matrix)arg1) -> int :
         |          docstring
         |  
         |  trace(...)
         |      trace( (Matrix)arg1) -> float :
         |          docstring
         |  
         |  transform(...)
         |      transform( (Matrix)arg1, (Matrix)arg2) -> None :
         |          docstring
         |      
         |      transform( (Matrix)arg1, (Matrix)arg2 [, (Matrix)arg3]) -> None :
         |          docstring
         |  
         |  vector_dot(...)
         |      vector_dot( (Matrix)arg1, (Matrix)arg2) -> float :
         |          docstring
         |  
         |  zero(...)
         |      zero( (Matrix)arg1) -> None :
         |          docstring
         |  
         |  zero_diagonal(...)
         |      zero_diagonal( (Matrix)arg1) -> None :
         |          docstring
         |  
         |  zero_lower(...)
         |      zero_lower( (Matrix)arg1) -> None :
         |          docstring
         |  
         |  zero_upper(...)
         |      zero_upper( (Matrix)arg1) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class MatrixFactory(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      MatrixFactory
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  create_matrix(...)
         |      create_matrix( (MatrixFactory)arg1) -> Matrix :
         |          docstring
         |      
         |      create_matrix( (MatrixFactory)arg1, (str)arg2) -> Matrix :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Static methods defined here:
         |  
         |  shared_object(...)
         |      shared_object() -> MatrixFactory :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class MintsHelper(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      MintsHelper
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  ao_angular_momentum(...)
         |      ao_angular_momentum( (MintsHelper)arg1) -> matrix_vector :
         |          docstring
         |  
         |  ao_erf_eri(...)
         |      ao_erf_eri( (MintsHelper)arg1, (float)arg2) -> Matrix :
         |          docstring
         |  
         |  ao_eri(...)
         |      ao_eri( (MintsHelper)arg1) -> Matrix :
         |          docstring
         |  
         |  ao_kinetic(...)
         |      ao_kinetic( (MintsHelper)arg1) -> Matrix :
         |          docstring
         |  
         |  ao_nabla(...)
         |      ao_nabla( (MintsHelper)arg1) -> matrix_vector :
         |          docstring
         |  
         |  ao_overlap(...)
         |      ao_overlap( (MintsHelper)arg1) -> Matrix :
         |          docstring
         |  
         |  ao_potential(...)
         |      ao_potential( (MintsHelper)arg1) -> Matrix :
         |          docstring
         |  
         |  basisset(...)
         |      basisset( (MintsHelper)arg1) -> BasisSet :
         |          docstring
         |  
         |  cdsalcs(...)
         |      cdsalcs( (MintsHelper)arg1, (int)arg2, (bool)arg3, (bool)arg4) -> CdSalcList :
         |          docstring
         |  
         |  factory(...)
         |      factory( (MintsHelper)arg1) -> MatrixFactory :
         |          docstring
         |  
         |  integrals(...)
         |      integrals( (MintsHelper)arg1) -> None :
         |          docstring
         |  
         |  one_electron_integrals(...)
         |      one_electron_integrals( (MintsHelper)arg1) -> None :
         |          docstring
         |      
         |      one_electron_integrals( (MintsHelper)arg1) -> None :
         |          docstring
         |  
         |  petite_list(...)
         |      petite_list( (MintsHelper)arg1) -> PetiteList :
         |          docstring
         |  
         |  play(...)
         |      play( (MintsHelper)arg1) -> None :
         |          docstring
         |  
         |  so_angular_momentum(...)
         |      so_angular_momentum( (MintsHelper)arg1) -> matrix_vector :
         |          docstring
         |  
         |  so_dipole(...)
         |      so_dipole( (MintsHelper)arg1) -> matrix_vector :
         |          docstring
         |  
         |  so_kinetic(...)
         |      so_kinetic( (MintsHelper)arg1) -> Matrix :
         |          docstring
         |  
         |  so_nabla(...)
         |      so_nabla( (MintsHelper)arg1) -> matrix_vector :
         |          docstring
         |  
         |  so_overlap(...)
         |      so_overlap( (MintsHelper)arg1) -> Matrix :
         |          docstring
         |  
         |  so_potential(...)
         |      so_potential( (MintsHelper)arg1) -> Matrix :
         |          docstring
         |  
         |  so_quadrupole(...)
         |      so_quadrupole( (MintsHelper)arg1) -> matrix_vector :
         |          docstring
         |  
         |  so_traceless_quadrupole(...)
         |      so_traceless_quadrupole( (MintsHelper)arg1) -> matrix_vector :
         |          docstring
         |  
         |  sobasisset(...)
         |      sobasisset( (MintsHelper)arg1) -> SOBasisSet :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class MoldenWriter(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      MoldenWriter
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1, (Wavefunction)arg2) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  write(...)
         |      write( (MoldenWriter)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class Molecule(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      Molecule
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  Z(...)
         |      Z( (Molecule)arg1, (int)arg2) -> float :
         |          docstring
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  activate_all_fragments(...)
         |      activate_all_fragments( (Molecule)arg1) -> None :
         |          docstring
         |  
         |  add_atom(...)
         |      add_atom( (Molecule)arg1, (int)arg2, (float)arg3, (float)arg4, (float)arg5, (str)arg6, (float)arg7, (float)arg8, (int)arg9) -> None :
         |          docstring
         |  
         |  atom_at_position(...)
         |      atom_at_position( (Molecule)arg1, (float)arg2, (float)arg3) -> int :
         |          docstring
         |  
         |  center_of_mass(...)
         |      center_of_mass( (Molecule)arg1) -> Vector3 :
         |          docstring
         |  
         |  charge(...)
         |      charge( (Molecule)arg1, (int)arg2) -> float :
         |          docstring
         |  
         |  deactivate_all_fragments(...)
         |      deactivate_all_fragments( (Molecule)arg1) -> None :
         |          docstring
         |  
         |  extract_subsets(...)
         |      extract_subsets( (Molecule)arg1, (list)arg2, (list)arg3) -> Molecule :
         |          docstring
         |      
         |      extract_subsets( (Molecule)arg1, (list)arg2, (int)arg3) -> Molecule :
         |          docstring
         |      
         |      extract_subsets( (Molecule)arg1, (int)arg2, (list)arg3) -> Molecule :
         |          docstring
         |      
         |      extract_subsets( (Molecule)arg1, (int)arg2, (int)arg3) -> Molecule :
         |          docstring
         |      
         |      extract_subsets( (Molecule)arg1, (list)arg2) -> Molecule :
         |          docstring
         |      
         |      extract_subsets( (Molecule)arg1, (int)arg2) -> Molecule :
         |          docstring
         |  
         |  find_point_group(...)
         |      find_point_group( (Molecule)arg1, (float)arg2) -> PointGroup :
         |          docstring
         |  
         |  fix_orientation(...)
         |      fix_orientation( (Molecule)arg1, (bool)arg2) -> None :
         |          docstring
         |  
         |  form_symmetry_information(...)
         |      form_symmetry_information( (Molecule)arg1, (float)arg2) -> None :
         |          docstring
         |  
         |  get_variable(...)
         |      get_variable( (Molecule)arg1, (str)arg2) -> float :
         |          docstring
         |  
         |  init_with_checkpoint(...)
         |      init_with_checkpoint( (Molecule)arg1, (Checkpoint)arg2) -> None :
         |          docstring
         |  
         |  init_with_io(...)
         |      init_with_io( (Molecule)arg1, (IO)arg2) -> None :
         |          docstring
         |  
         |  is_variable(...)
         |      is_variable( (Molecule)arg1, (str)arg2) -> bool :
         |          docstring
         |  
         |  label(...)
         |      label( (Molecule)arg1, (int)arg2) -> str :
         |          docstring
         |  
         |  mass(...)
         |      mass( (Molecule)arg1, (int)arg2) -> float :
         |          docstring
         |  
         |  molecular_charge(...)
         |      molecular_charge( (Molecule)arg1) -> int :
         |          docstring
         |  
         |  move_to_com(...)
         |      move_to_com( (Molecule)arg1) -> None :
         |          docstring
         |  
         |  multiplicity(...)
         |      multiplicity( (Molecule)arg1) -> int :
         |          docstring
         |  
         |  name(...)
         |      name( (Molecule)arg1) -> str :
         |          docstring
         |  
         |  natom(...)
         |      natom( (Molecule)arg1) -> int :
         |          docstring
         |  
         |  nfragments(...)
         |      nfragments( (Molecule)arg1) -> int :
         |          docstring
         |  
         |  nuclear_repulsion_energy(...)
         |      nuclear_repulsion_energy( (Molecule)arg1) -> float :
         |          docstring
         |  
         |  point_group(...)
         |      point_group( (Molecule)arg1) -> PointGroup :
         |          docstring
         |  
         |  print_in_input_format(...)
         |      print_in_input_format( (Molecule)arg1) -> None :
         |          docstring
         |  
         |  print_out(...)
         |      print_out( (Molecule)arg1) -> None :
         |          docstring
         |      
         |      print_out( (Molecule)arg1) -> None :
         |          docstring
         |  
         |  print_out_in_bohr(...)
         |      print_out_in_bohr( (Molecule)arg1) -> None :
         |          docstring
         |  
         |  reinterpret_coordentry(...)
         |      reinterpret_coordentry( (Molecule)arg1, (bool)arg2) -> None :
         |          docstring
         |  
         |  reset_point_group(...)
         |      reset_point_group( (Molecule)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  save_string_xyz(...)
         |      save_string_xyz( (Molecule)arg1) -> str :
         |          docstring
         |  
         |  save_to_checkpoint(...)
         |      save_to_checkpoint( (Molecule)arg1, (Checkpoint)arg2, (str)arg3) -> None :
         |          docstring
         |  
         |  save_xyz(...)
         |      save_xyz( (Molecule)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  schoenflies_symbol(...)
         |      schoenflies_symbol( (Molecule)arg1) -> str :
         |          docstring
         |  
         |  set_active_fragment(...)
         |      set_active_fragment( (Molecule)arg1, (int)arg2) -> None :
         |          docstring
         |  
         |  set_active_fragments(...)
         |      set_active_fragments( (Molecule)arg1, (list)arg2) -> None :
         |          docstring
         |  
         |  set_basis_all_atoms(...)
         |      set_basis_all_atoms( (Molecule)arg1, (str)arg2, (str)arg3) -> None :
         |          docstring
         |  
         |  set_basis_by_label(...)
         |      set_basis_by_label( (Molecule)arg1, (str)arg2, (str)arg3, (str)arg4) -> None :
         |          docstring
         |  
         |  set_basis_by_number(...)
         |      set_basis_by_number( (Molecule)arg1, (int)arg2, (str)arg3, (str)arg4) -> None :
         |          docstring
         |  
         |  set_basis_by_symbol(...)
         |      set_basis_by_symbol( (Molecule)arg1, (str)arg2, (str)arg3, (str)arg4) -> None :
         |          docstring
         |  
         |  set_geometry(...)
         |      set_geometry( (Molecule)arg1, (Matrix)arg2) -> None :
         |          docstring
         |  
         |  set_ghost_fragment(...)
         |      set_ghost_fragment( (Molecule)arg1, (int)arg2) -> None :
         |          docstring
         |  
         |  set_ghost_fragments(...)
         |      set_ghost_fragments( (Molecule)arg1, (list)arg2) -> None :
         |          docstring
         |  
         |  set_molecular_charge(...)
         |      set_molecular_charge( (Molecule)arg1, (int)arg2) -> None :
         |          docstring
         |  
         |  set_multiplicity(...)
         |      set_multiplicity( (Molecule)arg1, (int)arg2) -> None :
         |          docstring
         |  
         |  set_name(...)
         |      set_name( (Molecule)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  set_point_group(...)
         |      set_point_group( (Molecule)arg1, (PointGroup)arg2) -> None :
         |          docstring
         |  
         |  set_variable(...)
         |      set_variable( (Molecule)arg1, (str)arg2, (float)arg3) -> None :
         |          docstring
         |  
         |  symbol(...)
         |      symbol( (Molecule)arg1, (int)arg2) -> str :
         |          docstring
         |  
         |  translate(...)
         |      translate( (Molecule)arg1, (Vector3)arg2) -> None :
         |          docstring
         |  
         |  update_geometry(...)
         |      update_geometry( (Molecule)arg1) -> None :
         |          docstring
         |      
         |      update_geometry( (Molecule)arg1) -> None :
         |          docstring
         |  
         |  x(...)
         |      x( (Molecule)arg1, (int)arg2) -> float :
         |          docstring
         |  
         |  y(...)
         |      y( (Molecule)arg1, (int)arg2) -> float :
         |          docstring
         |  
         |  z(...)
         |      z( (Molecule)arg1, (int)arg2) -> float :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Static methods defined here:
         |  
         |  create_molecule_from_string(...)
         |      create_molecule_from_string( (str)arg1) -> Molecule :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class MultipoleSymmetry(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      MultipoleSymmetry
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1, (int)arg2, (Molecule)arg3, (object)arg4, (MatrixFactory)arg5) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  create_matrices(...)
         |      create_matrices( (MultipoleSymmetry)arg1, (str)arg2) -> matrix_vector :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class NBOWriter(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      NBOWriter
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1, (Wavefunction)arg2) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  write(...)
         |      write( (NBOWriter)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class OEProp(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      OEProp
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  add(...)
         |      add( (OEProp)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  compute(...)
         |      compute( (OEProp)arg1) -> None :
         |          docstring
         |  
         |  set_title(...)
         |      set_title( (OEProp)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class PetiteList(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      PetiteList
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  aotoso(...)
         |      aotoso( (PetiteList)arg1) -> Matrix :
         |          docstring
         |  
         |  print(...)
         |      print( (PetiteList)arg1, (object)arg2) -> None :
         |          docstring
         |  
         |  sotoao(...)
         |      sotoao( (PetiteList)arg1) -> Matrix :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __init__ = <built-in function __init__>
         |      Raises an exception
         |      This class cannot be instantiated from Python
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class PointGroup(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      PointGroup
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |      
         |      __init__( (object)arg1, (str)arg2) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  symbol(...)
         |      symbol( (PointGroup)arg1) -> str :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class Process(Boost.Python.instance)
         |  Method resolution order:
         |      Process
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors defined here:
         |  
         |  environment
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 24
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class PseudoTrial(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      PseudoTrial
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  getA(...)
         |      getA( (PseudoTrial)arg1) -> Matrix :
         |          docstring
         |  
         |  getI(...)
         |      getI( (PseudoTrial)arg1) -> Matrix :
         |          docstring
         |  
         |  getIPS(...)
         |      getIPS( (PseudoTrial)arg1) -> Matrix :
         |          docstring
         |  
         |  getQ(...)
         |      getQ( (PseudoTrial)arg1) -> Matrix :
         |          docstring
         |  
         |  getR(...)
         |      getR( (PseudoTrial)arg1) -> Matrix :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class PsiReturnType(Boost.Python.enum)
         |  docstring
         |  
         |  Method resolution order:
         |      PsiReturnType
         |      Boost.Python.enum
         |      __builtin__.int
         |      __builtin__.object
         |  
         |  Data and other attributes defined here:
         |  
         |  Balk = PsiMod.PsiReturnType.Balk
         |  
         |  EndLoop = PsiMod.PsiReturnType.EndLoop
         |  
         |  Failure = PsiMod.PsiReturnType.Failure
         |  
         |  Success = PsiMod.PsiReturnType.Success
         |  
         |  names = {'Balk': PsiMod.PsiReturnType.Balk, 'EndLoop': PsiMod.PsiRetur...
         |  
         |  values = {0: PsiMod.PsiReturnType.Success, 1: PsiMod.PsiReturnType.Fai...
         |  
         |  ----------------------------------------------------------------------
         |  Methods inherited from Boost.Python.enum:
         |  
         |  __repr__(...)
         |      x.__repr__() <==> repr(x)
         |  
         |  __str__(...)
         |      x.__str__() <==> str(x)
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.enum:
         |  
         |  name
         |  
         |  ----------------------------------------------------------------------
         |  Methods inherited from __builtin__.int:
         |  
         |  __abs__(...)
         |      x.__abs__() <==> abs(x)
         |  
         |  __add__(...)
         |      x.__add__(y) <==> x+y
         |  
         |  __and__(...)
         |      x.__and__(y) <==> x&y
         |  
         |  __cmp__(...)
         |      x.__cmp__(y) <==> cmp(x,y)
         |  
         |  __coerce__(...)
         |      x.__coerce__(y) <==> coerce(x, y)
         |  
         |  __div__(...)
         |      x.__div__(y) <==> x/y
         |  
         |  __divmod__(...)
         |      x.__divmod__(y) <==> divmod(x, y)
         |  
         |  __float__(...)
         |      x.__float__() <==> float(x)
         |  
         |  __floordiv__(...)
         |      x.__floordiv__(y) <==> x//y
         |  
         |  __format__(...)
         |  
         |  __getattribute__(...)
         |      x.__getattribute__('name') <==> x.name
         |  
         |  __getnewargs__(...)
         |  
         |  __hash__(...)
         |      x.__hash__() <==> hash(x)
         |  
         |  __hex__(...)
         |      x.__hex__() <==> hex(x)
         |  
         |  __index__(...)
         |      x[y:z] <==> x[y.__index__():z.__index__()]
         |  
         |  __int__(...)
         |      x.__int__() <==> int(x)
         |  
         |  __invert__(...)
         |      x.__invert__() <==> ~x
         |  
         |  __long__(...)
         |      x.__long__() <==> long(x)
         |  
         |  __lshift__(...)
         |      x.__lshift__(y) <==> x<<y
         |  
         |  __mod__(...)
         |      x.__mod__(y) <==> x%y
         |  
         |  __mul__(...)
         |      x.__mul__(y) <==> x*y
         |  
         |  __neg__(...)
         |      x.__neg__() <==> -x
         |  
         |  __nonzero__(...)
         |      x.__nonzero__() <==> x != 0
         |  
         |  __oct__(...)
         |      x.__oct__() <==> oct(x)
         |  
         |  __or__(...)
         |      x.__or__(y) <==> x|y
         |  
         |  __pos__(...)
         |      x.__pos__() <==> +x
         |  
         |  __pow__(...)
         |      x.__pow__(y[, z]) <==> pow(x, y[, z])
         |  
         |  __radd__(...)
         |      x.__radd__(y) <==> y+x
         |  
         |  __rand__(...)
         |      x.__rand__(y) <==> y&x
         |  
         |  __rdiv__(...)
         |      x.__rdiv__(y) <==> y/x
         |  
         |  __rdivmod__(...)
         |      x.__rdivmod__(y) <==> divmod(y, x)
         |  
         |  __rfloordiv__(...)
         |      x.__rfloordiv__(y) <==> y//x
         |  
         |  __rlshift__(...)
         |      x.__rlshift__(y) <==> y<<x
         |  
         |  __rmod__(...)
         |      x.__rmod__(y) <==> y%x
         |  
         |  __rmul__(...)
         |      x.__rmul__(y) <==> y*x
         |  
         |  __ror__(...)
         |      x.__ror__(y) <==> y|x
         |  
         |  __rpow__(...)
         |      y.__rpow__(x[, z]) <==> pow(x, y[, z])
         |  
         |  __rrshift__(...)
         |      x.__rrshift__(y) <==> y>>x
         |  
         |  __rshift__(...)
         |      x.__rshift__(y) <==> x>>y
         |  
         |  __rsub__(...)
         |      x.__rsub__(y) <==> y-x
         |  
         |  __rtruediv__(...)
         |      x.__rtruediv__(y) <==> y/x
         |  
         |  __rxor__(...)
         |      x.__rxor__(y) <==> y^x
         |  
         |  __sub__(...)
         |      x.__sub__(y) <==> x-y
         |  
         |  __truediv__(...)
         |      x.__truediv__(y) <==> x/y
         |  
         |  __trunc__(...)
         |      Truncating an Integral returns itself.
         |  
         |  __xor__(...)
         |      x.__xor__(y) <==> x^y
         |  
         |  bit_length(...)
         |      int.bit_length() -> int
         |      
         |      Number of bits necessary to represent self in binary.
         |      >>> bin(37)
         |      '0b100101'
         |      >>> (37).bit_length()
         |      6
         |  
         |  conjugate(...)
         |      Returns self, the complex conjugate of any int.
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from __builtin__.int:
         |  
         |  denominator
         |      the denominator of a rational number in lowest terms
         |  
         |  imag
         |      the imaginary part of a complex number
         |  
         |  numerator
         |      the numerator of a rational number in lowest terms
         |  
         |  real
         |      the real part of a complex number
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from __builtin__.int:
         |  
         |  __new__ = <built-in method __new__ of type object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class RHF(HF, Wavefunction)
         |  docstring
         |  
         |  Method resolution order:
         |      RHF
         |      HF
         |      Wavefunction
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __init__ = <built-in function __init__>
         |      Raises an exception
         |      This class cannot be instantiated from Python
         |  
         |  ----------------------------------------------------------------------
         |  Methods inherited from Wavefunction:
         |  
         |  Ca(...)
         |      Ca( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Cb(...)
         |      Cb( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Da(...)
         |      Da( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Db(...)
         |      Db( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Fa(...)
         |      Fa( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Fb(...)
         |      Fb( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  add_postiteration_callback(...)
         |      add_postiteration_callback( (Wavefunction)arg1, (object)arg2) -> None :
         |          docstring
         |  
         |  add_preiteration_callback(...)
         |      add_preiteration_callback( (Wavefunction)arg1, (object)arg2) -> None :
         |          docstring
         |  
         |  basisset(...)
         |      basisset( (Wavefunction)arg1) -> BasisSet :
         |          docstring
         |  
         |  energy(...)
         |      energy( (Wavefunction)arg1) -> float :
         |          docstring
         |  
         |  epsilon_a(...)
         |      epsilon_a( (Wavefunction)arg1) -> Vector :
         |          docstring
         |  
         |  epsilon_b(...)
         |      epsilon_b( (Wavefunction)arg1) -> Vector :
         |          docstring
         |  
         |  frequencies(...)
         |      frequencies( (Wavefunction)arg1) -> Vector :
         |          docstring
         |  
         |  gradient(...)
         |      gradient( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  nirrep(...)
         |      nirrep( (Wavefunction)arg1) -> int :
         |          docstring
         |  
         |  nmo(...)
         |      nmo( (Wavefunction)arg1) -> int :
         |          docstring
         |  
         |  nso(...)
         |      nso( (Wavefunction)arg1) -> int :
         |          docstring
         |  
         |  sobasisset(...)
         |      sobasisset( (Wavefunction)arg1) -> SOBasisSet :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class SOBasisSet(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      SOBasisSet
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  petite_list(...)
         |      petite_list( (SOBasisSet)arg1) -> PetiteList :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __init__ = <built-in function __init__>
         |      Raises an exception
         |      This class cannot be instantiated from Python
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class SuperFunctional(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      SuperFunctional
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  computeRKSFunctional(...)
         |      computeRKSFunctional( (SuperFunctional)arg1, (object)arg2) -> None :
         |          docstring
         |  
         |  computeUKSFunctional(...)
         |      computeUKSFunctional( (SuperFunctional)arg1, (object)arg2) -> None :
         |          docstring
         |  
         |  get_citation(...)
         |      get_citation( (SuperFunctional)arg1) -> str :
         |          docstring
         |  
         |  get_composition(...)
         |      get_composition( (SuperFunctional)arg1) -> str :
         |          docstring
         |  
         |  get_dash_d(...)
         |      get_dash_d( (SuperFunctional)arg1) -> object :
         |          docstring
         |  
         |  get_deriv(...)
         |      get_deriv( (SuperFunctional)arg1) -> int :
         |          docstring
         |  
         |  get_description(...)
         |      get_description( (SuperFunctional)arg1) -> str :
         |          docstring
         |  
         |  get_exact_exchange(...)
         |      get_exact_exchange( (SuperFunctional)arg1) -> float :
         |          docstring
         |  
         |  get_functional(...)
         |      get_functional( (SuperFunctional)arg1, (int)arg2) -> Functional :
         |          docstring
         |  
         |  get_name(...)
         |      get_name( (SuperFunctional)arg1) -> str :
         |          docstring
         |  
         |  get_npoints(...)
         |      get_npoints( (SuperFunctional)arg1) -> int :
         |          docstring
         |  
         |  get_omega(...)
         |      get_omega( (SuperFunctional)arg1) -> float :
         |          docstring
         |  
         |  get_pt2(...)
         |      get_pt2( (SuperFunctional)arg1) -> float :
         |          docstring
         |  
         |  get_size(...)
         |      get_size( (SuperFunctional)arg1) -> int :
         |          docstring
         |  
         |  get_weight(...)
         |      get_weight( (SuperFunctional)arg1, (int)arg2) -> float :
         |          docstring
         |  
         |  is_dash_d(...)
         |      is_dash_d( (SuperFunctional)arg1) -> bool :
         |          docstring
         |  
         |  is_double_hybrid(...)
         |      is_double_hybrid( (SuperFunctional)arg1) -> bool :
         |          docstring
         |  
         |  is_gga(...)
         |      is_gga( (SuperFunctional)arg1) -> bool :
         |          docstring
         |  
         |  is_hybrid(...)
         |      is_hybrid( (SuperFunctional)arg1) -> bool :
         |          docstring
         |  
         |  is_meta(...)
         |      is_meta( (SuperFunctional)arg1) -> bool :
         |          docstring
         |  
         |  is_range_corrected(...)
         |      is_range_corrected( (SuperFunctional)arg1) -> bool :
         |          docstring
         |  
         |  set_citation(...)
         |      set_citation( (SuperFunctional)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  set_dash_d(...)
         |      set_dash_d( (SuperFunctional)arg1, (object)arg2, (float)arg3) -> None :
         |          docstring
         |  
         |  set_deriv(...)
         |      set_deriv( (SuperFunctional)arg1, (int)arg2) -> None :
         |          docstring
         |  
         |  set_description(...)
         |      set_description( (SuperFunctional)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  set_exact_exchange(...)
         |      set_exact_exchange( (SuperFunctional)arg1, (float)arg2) -> None :
         |          docstring
         |  
         |  set_name(...)
         |      set_name( (SuperFunctional)arg1, (str)arg2) -> None :
         |          docstring
         |  
         |  set_npoints(...)
         |      set_npoints( (SuperFunctional)arg1, (int)arg2) -> None :
         |          docstring
         |  
         |  set_omega(...)
         |      set_omega( (SuperFunctional)arg1, (float)arg2) -> None :
         |          docstring
         |  
         |  set_parameter(...)
         |      set_parameter( (SuperFunctional)arg1, (str)arg2, (str)arg3, (float)arg4) -> None :
         |          docstring
         |  
         |  set_pt2(...)
         |      set_pt2( (SuperFunctional)arg1, (float)arg2) -> None :
         |          docstring
         |  
         |  set_size(...)
         |      set_size( (SuperFunctional)arg1) -> int :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Static methods defined here:
         |  
         |  available_names(...)
         |      available_names() -> object :
         |          docstring
         |  
         |  available_superfunctionals(...)
         |      available_superfunctionals() -> str :
         |          docstring
         |  
         |  build_superfunctional(...)
         |      build_superfunctional( (str)arg1, (int)arg2, (int)arg3) -> SuperFunctional :
         |          docstring
         |  
         |  create_superfunctional(...)
         |      create_superfunctional( (str)arg1, (int)arg2, (int)arg3) -> SuperFunctional :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class SymmetryOperation(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      SymmetryOperation
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  E(...)
         |      E( (SymmetryOperation)arg1) -> None :
         |          docstring
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |      
         |      __init__( (object)arg1, (SymmetryOperation)arg2) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  c2_x(...)
         |      c2_x( (SymmetryOperation)arg1) -> None :
         |          docstring
         |  
         |  c2_y(...)
         |      c2_y( (SymmetryOperation)arg1) -> None :
         |          docstring
         |  
         |  i(...)
         |      i( (SymmetryOperation)arg1) -> None :
         |          docstring
         |  
         |  operate(...)
         |      operate( (SymmetryOperation)arg1, (SymmetryOperation)arg2) -> SymmetryOperation :
         |          docstring
         |  
         |  rotate_n(...)
         |      rotate_n( (SymmetryOperation)arg1, (int)arg2) -> None :
         |          docstring
         |  
         |  rotate_theta(...)
         |      rotate_theta( (SymmetryOperation)arg1, (float)arg2) -> None :
         |          docstring
         |  
         |  sigma_xy(...)
         |      sigma_xy( (SymmetryOperation)arg1) -> None :
         |          docstring
         |  
         |  sigma_xz(...)
         |      sigma_xz( (SymmetryOperation)arg1) -> None :
         |          docstring
         |  
         |  sigma_yz(...)
         |      sigma_yz( (SymmetryOperation)arg1) -> None :
         |          docstring
         |  
         |  trace(...)
         |      trace( (SymmetryOperation)arg1) -> float :
         |          docstring
         |  
         |  transform(...)
         |      transform( (SymmetryOperation)arg1, (SymmetryOperation)arg2) -> SymmetryOperation :
         |          docstring
         |  
         |  transpose(...)
         |      transpose( (SymmetryOperation)arg1) -> None :
         |          docstring
         |  
         |  unit(...)
         |      unit( (SymmetryOperation)arg1) -> None :
         |          docstring
         |  
         |  zero(...)
         |      zero( (SymmetryOperation)arg1) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 96
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class Vector(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      Vector
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __getitem__(...)
         |      __getitem__( (Vector)arg1, (int)arg2) -> float :
         |          docstring
         |      
         |      __getitem__( (Vector)arg1, (tuple)arg2) -> float :
         |          docstring
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |      
         |      __init__( (object)arg1, (int)arg2) -> None
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  __setitem__(...)
         |      __setitem__( (Vector)arg1, (int)arg2, (float)arg3) -> None :
         |          docstring
         |      
         |      __setitem__( (Vector)arg1, (tuple)arg2, (float)arg3) -> None :
         |          docstring
         |  
         |  dim(...)
         |      dim( (Vector)arg1, (int)arg2) -> int :
         |          docstring
         |  
         |  get(...)
         |      get( (Vector)arg1, (int)arg2) -> float :
         |          docstring
         |      
         |      get( (Vector)arg1, (int)arg2, (int)arg3) -> float :
         |          docstring
         |  
         |  nirrep(...)
         |      nirrep( (Vector)arg1) -> int :
         |          docstring
         |  
         |  print_out(...)
         |      print_out( (Vector)arg1) -> None :
         |          docstring
         |  
         |  scale(...)
         |      scale( (Vector)arg1, (float)arg2) -> None :
         |          docstring
         |  
         |  set(...)
         |      set( (Vector)arg1, (int)arg2, (float)arg3) -> None :
         |          docstring
         |      
         |      set( (Vector)arg1, (int)arg2, (int)arg3, (float)arg4) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 32
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class Vector3(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      Vector3
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __add__(...)
         |      __add__( (Vector3)arg1, (Vector3)arg2) -> object
         |  
         |  __getitem__(...)
         |      __getitem__( (Vector3)arg1, (int)arg2) -> float :
         |          docstring
         |  
         |  __iadd__(...)
         |      __iadd__( (object)arg1, (Vector3)arg2) -> object
         |  
         |  __imul__(...)
         |      __imul__( (object)arg1, (float)arg2) -> object
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |      
         |      __init__( (object)arg1, (float)arg2) -> None
         |      
         |      __init__( (object)arg1, (float)arg2, (float)arg3, (float)arg4) -> None
         |      
         |      __init__( (object)arg1, (Vector3)arg2) -> None
         |  
         |  __isub__(...)
         |      __isub__( (object)arg1, (Vector3)arg2) -> object
         |  
         |  __neg__(...)
         |      __neg__( (Vector3)arg1) -> object
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  __str__(...)
         |      __str__( (Vector3)arg1) -> str :
         |          docstring
         |  
         |  __sub__(...)
         |      __sub__( (Vector3)arg1, (Vector3)arg2) -> object
         |  
         |  cross(...)
         |      cross( (Vector3)arg1, (Vector3)arg2) -> Vector3 :
         |          docstring
         |  
         |  distance(...)
         |      distance( (Vector3)arg1, (Vector3)arg2) -> float :
         |          docstring
         |  
         |  dot(...)
         |      dot( (Vector3)arg1, (Vector3)arg2) -> float :
         |          docstring
         |  
         |  norm(...)
         |      norm( (Vector3)arg1) -> float :
         |          docstring
         |  
         |  normalize(...)
         |      normalize( (Vector3)arg1) -> None :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 40
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class Wavefunction(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      Wavefunction
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  Ca(...)
         |      Ca( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Cb(...)
         |      Cb( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Da(...)
         |      Da( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Db(...)
         |      Db( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Fa(...)
         |      Fa( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  Fb(...)
         |      Fb( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  add_postiteration_callback(...)
         |      add_postiteration_callback( (Wavefunction)arg1, (object)arg2) -> None :
         |          docstring
         |  
         |  add_preiteration_callback(...)
         |      add_preiteration_callback( (Wavefunction)arg1, (object)arg2) -> None :
         |          docstring
         |  
         |  basisset(...)
         |      basisset( (Wavefunction)arg1) -> BasisSet :
         |          docstring
         |  
         |  energy(...)
         |      energy( (Wavefunction)arg1) -> float :
         |          docstring
         |  
         |  epsilon_a(...)
         |      epsilon_a( (Wavefunction)arg1) -> Vector :
         |          docstring
         |  
         |  epsilon_b(...)
         |      epsilon_b( (Wavefunction)arg1) -> Vector :
         |          docstring
         |  
         |  frequencies(...)
         |      frequencies( (Wavefunction)arg1) -> Vector :
         |          docstring
         |  
         |  gradient(...)
         |      gradient( (Wavefunction)arg1) -> Matrix :
         |          docstring
         |  
         |  nirrep(...)
         |      nirrep( (Wavefunction)arg1) -> int :
         |          docstring
         |  
         |  nmo(...)
         |      nmo( (Wavefunction)arg1) -> int :
         |          docstring
         |  
         |  nso(...)
         |      nso( (Wavefunction)arg1) -> int :
         |          docstring
         |  
         |  sobasisset(...)
         |      sobasisset( (Wavefunction)arg1) -> SOBasisSet :
         |          docstring
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __init__ = <built-in function __init__>
         |      Raises an exception
         |      This class cannot be instantiated from Python
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
        
        class matrix_vector(Boost.Python.instance)
         |  docstring
         |  
         |  Method resolution order:
         |      matrix_vector
         |      Boost.Python.instance
         |      __builtin__.object
         |  
         |  Methods defined here:
         |  
         |  __contains__(...)
         |      __contains__( (matrix_vector)arg1, (object)arg2) -> bool
         |  
         |  __delitem__(...)
         |      __delitem__( (matrix_vector)arg1, (object)arg2) -> None
         |  
         |  __getitem__(...)
         |      __getitem__( (object)arg1, (object)arg2) -> object
         |  
         |  __init__(...)
         |      __init__( (object)arg1) -> None
         |  
         |  __iter__(...)
         |      __iter__( (object)arg1) -> object
         |  
         |  __len__(...)
         |      __len__( (matrix_vector)arg1) -> int
         |  
         |  __reduce__ = <unnamed Boost.Python function>(...)
         |  
         |  __setitem__(...)
         |      __setitem__( (matrix_vector)arg1, (object)arg2, (object)arg3) -> None
         |  
         |  append(...)
         |      append( (matrix_vector)arg1, (object)arg2) -> None
         |  
         |  extend(...)
         |      extend( (matrix_vector)arg1, (object)arg2) -> None
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes defined here:
         |  
         |  __instance_size__ = 40
         |  
         |  ----------------------------------------------------------------------
         |  Data descriptors inherited from Boost.Python.instance:
         |  
         |  __dict__
         |  
         |  __weakref__
         |  
         |  ----------------------------------------------------------------------
         |  Data and other attributes inherited from Boost.Python.instance:
         |  
         |  __new__ = <built-in method __new__ of Boost.Python.class object>
         |      T.__new__(S, ...) -> a new object with type S, a subtype of T
    
    FUNCTIONS
        DASUM(...)
            DASUM( (int)arg1, (int)arg2, (Vector)arg3, (int)arg4) -> float :
                docstring
        
        DAXPY(...)
            DAXPY( (int)arg1, (int)arg2, (float)arg3, (Vector)arg4, (int)arg5, (Vector)arg6, (int)arg7) -> None :
                docstring
        
        DCOPY(...)
            DCOPY( (int)arg1, (int)arg2, (Vector)arg3, (int)arg4, (Vector)arg5, (int)arg6) -> None :
                docstring
        
        DDOT(...)
            DDOT( (int)arg1, (int)arg2, (Vector)arg3, (int)arg4, (Vector)arg5, (int)arg6) -> float :
                docstring
        
        DGBMV(...)
            DGBMV( (int)arg1, (str)arg2, (int)arg3, (int)arg4, (int)arg5, (int)arg6, (float)arg7, (Matrix)arg8, (int)arg9, (Vector)arg10, (int)arg11, (float)arg12, (Vector)arg13, (int)arg14) -> None :
                docstring
        
        DGEEV(...)
            DGEEV( (int)arg1, (str)arg2, (str)arg3, (int)arg4, (Matrix)arg5, (int)arg6, (Vector)arg7, (Vector)arg8, (Matrix)arg9, (int)arg10, (Matrix)arg11, (int)arg12, (Vector)arg13, (int)arg14) -> int :
                docstring
        
        DGEMM(...)
            DGEMM( (int)arg1, (str)arg2, (str)arg3, (int)arg4, (int)arg5, (int)arg6, (float)arg7, (Matrix)arg8, (int)arg9, (Matrix)arg10, (int)arg11, (float)arg12, (Matrix)arg13, (int)arg14) -> None :
                docstring
        
        DGEMV(...)
            DGEMV( (int)arg1, (str)arg2, (int)arg3, (int)arg4, (float)arg5, (Matrix)arg6, (int)arg7, (Vector)arg8, (int)arg9, (float)arg10, (Vector)arg11, (int)arg12) -> None :
                docstring
        
        DGER(...)
            DGER( (int)arg1, (int)arg2, (int)arg3, (float)arg4, (Vector)arg5, (int)arg6, (Vector)arg7, (int)arg8, (Matrix)arg9, (int)arg10) -> None :
                docstring
        
        DGETRF(...)
            DGETRF( (int)arg1, (int)arg2, (int)arg3, (Matrix)arg4, (int)arg5, (IntVector)arg6) -> int :
                docstring
        
        DGETRI(...)
            DGETRI( (int)arg1, (int)arg2, (Matrix)arg3, (int)arg4, (IntVector)arg5, (Vector)arg6, (int)arg7) -> int :
                docstring
        
        DGETRS(...)
            DGETRS( (int)arg1, (str)arg2, (int)arg3, (int)arg4, (Matrix)arg5, (int)arg6, (IntVector)arg7, (Matrix)arg8, (int)arg9) -> int :
                docstring
        
        DNRM2(...)
            DNRM2( (int)arg1, (int)arg2, (Vector)arg3, (int)arg4) -> float :
                docstring
        
        DPOTRF(...)
            DPOTRF( (int)arg1, (str)arg2, (int)arg3, (Matrix)arg4, (int)arg5) -> int :
                docstring
        
        DPOTRI(...)
            DPOTRI( (int)arg1, (str)arg2, (int)arg3, (Matrix)arg4, (int)arg5) -> int :
                docstring
        
        DPOTRS(...)
            DPOTRS( (int)arg1, (str)arg2, (int)arg3, (int)arg4, (Matrix)arg5, (int)arg6, (Matrix)arg7, (int)arg8) -> int :
                docstring
        
        DROT(...)
            DROT( (int)arg1, (int)arg2, (Vector)arg3, (int)arg4, (Vector)arg5, (int)arg6, (float)arg7, (float)arg8) -> None :
                docstring
        
        DSBMV(...)
            DSBMV( (int)arg1, (str)arg2, (int)arg3, (int)arg4, (float)arg5, (Matrix)arg6, (int)arg7, (Vector)arg8, (int)arg9, (float)arg10, (Vector)arg11, (int)arg12) -> None :
                docstring
        
        DSCAL(...)
            DSCAL( (int)arg1, (int)arg2, (float)arg3, (Vector)arg4, (int)arg5) -> None :
                docstring
        
        DSWAP(...)
            DSWAP( (int)arg1, (int)arg2, (Vector)arg3, (int)arg4, (Vector)arg5, (int)arg6) -> None :
                docstring
        
        DSYEV(...)
            DSYEV( (int)arg1, (str)arg2, (str)arg3, (int)arg4, (Matrix)arg5, (int)arg6, (Vector)arg7, (Vector)arg8, (int)arg9) -> int :
                docstring
        
        DSYMM(...)
            DSYMM( (int)arg1, (str)arg2, (str)arg3, (int)arg4, (int)arg5, (float)arg6, (Matrix)arg7, (int)arg8, (Matrix)arg9, (int)arg10, (float)arg11, (Matrix)arg12, (int)arg13) -> None :
                docstring
        
        DSYMV(...)
            DSYMV( (int)arg1, (str)arg2, (int)arg3, (float)arg4, (Matrix)arg5, (int)arg6, (Vector)arg7, (int)arg8, (float)arg9, (Vector)arg10, (int)arg11) -> None :
                docstring
        
        DSYR(...)
            DSYR( (int)arg1, (str)arg2, (int)arg3, (float)arg4, (Vector)arg5, (int)arg6, (Matrix)arg7, (int)arg8) -> None :
                docstring
        
        DSYR2(...)
            DSYR2( (int)arg1, (str)arg2, (int)arg3, (float)arg4, (Vector)arg5, (int)arg6, (Vector)arg7, (int)arg8, (Matrix)arg9, (int)arg10) -> None :
                docstring
        
        DSYR2K(...)
            DSYR2K( (int)arg1, (str)arg2, (str)arg3, (int)arg4, (int)arg5, (float)arg6, (Matrix)arg7, (int)arg8, (Matrix)arg9, (int)arg10, (float)arg11, (Matrix)arg12, (int)arg13) -> None :
                docstring
        
        DSYRK(...)
            DSYRK( (int)arg1, (str)arg2, (str)arg3, (int)arg4, (int)arg5, (float)arg6, (Matrix)arg7, (int)arg8, (float)arg9, (Matrix)arg10, (int)arg11) -> None :
                docstring
        
        DSYSV(...)
            DSYSV( (int)arg1, (str)arg2, (int)arg3, (int)arg4, (Matrix)arg5, (int)arg6, (IntVector)arg7, (Matrix)arg8, (int)arg9, (Vector)arg10, (int)arg11) -> int :
                docstring
        
        DTBMV(...)
            DTBMV( (int)arg1, (str)arg2, (str)arg3, (str)arg4, (int)arg5, (int)arg6, (Matrix)arg7, (int)arg8, (Vector)arg9, (int)arg10) -> None :
                docstring
        
        DTBSV(...)
            DTBSV( (int)arg1, (str)arg2, (str)arg3, (str)arg4, (int)arg5, (int)arg6, (Matrix)arg7, (int)arg8, (Vector)arg9, (int)arg10) -> None :
                docstring
        
        DTRMM(...)
            DTRMM( (int)arg1, (str)arg2, (str)arg3, (str)arg4, (str)arg5, (int)arg6, (int)arg7, (float)arg8, (Matrix)arg9, (int)arg10, (Matrix)arg11, (int)arg12) -> None :
                docstring
        
        DTRMV(...)
            DTRMV( (int)arg1, (str)arg2, (str)arg3, (str)arg4, (int)arg5, (Matrix)arg6, (int)arg7, (Vector)arg8, (int)arg9) -> None :
                docstring
        
        DTRSM(...)
            DTRSM( (int)arg1, (str)arg2, (str)arg3, (str)arg4, (str)arg5, (int)arg6, (int)arg7, (float)arg8, (Matrix)arg9, (int)arg10, (Matrix)arg11, (int)arg12) -> None :
                docstring
        
        DTRSV(...)
            DTRSV( (int)arg1, (str)arg2, (str)arg3, (str)arg4, (int)arg5, (Matrix)arg6, (int)arg7, (Vector)arg8, (int)arg9) -> None :
                docstring
        
        IDAMAX(...)
            IDAMAX( (int)arg1, (int)arg2, (Vector)arg3, (int)arg4) -> int :
                docstring
        
        adc(...)
            adc() -> float :
                docstring
        
        add_user_basis_file(...)
            add_user_basis_file( (str)arg1) -> None :
                docstring
        
        benchmark_blas1(...)
            benchmark_blas1( (int)arg1, (float)arg2) -> None :
                docstring
        
        benchmark_blas2(...)
            benchmark_blas2( (int)arg1, (float)arg2) -> None :
                docstring
        
        benchmark_blas3(...)
            benchmark_blas3( (int)arg1, (float)arg2, (int)arg3) -> None :
                docstring
        
        benchmark_disk(...)
            benchmark_disk( (int)arg1, (float)arg2) -> None :
                docstring
        
        benchmark_integrals(...)
            benchmark_integrals( (int)arg1, (float)arg2) -> None :
                docstring
        
        benchmark_math(...)
            benchmark_math( (float)arg1) -> None :
                docstring
        
        ccdensity(...)
            ccdensity() -> float :
                docstring
        
        ccenergy(...)
            ccenergy() -> float :
                docstring
        
        cceom(...)
            cceom() -> float :
                docstring
        
        cchbar(...)
            cchbar() -> float :
                docstring
        
        cclambda(...)
            cclambda() -> float :
                docstring
        
        ccresponse(...)
            ccresponse() -> float :
                docstring
        
        ccsort(...)
            ccsort() -> float :
                docstring
        
        cctriples(...)
            cctriples() -> float :
                docstring
        
        clean(...)
            clean() -> None :
                Function to remove scratch files. Call between independent jobs.
        
        close_outfile(...)
            close_outfile() -> None :
                docstring
        
        dcft(...)
            dcft() -> float :
                docstring
        
        deriv(...)
            deriv() -> int :
                docstring
        
        detci(...)
            detci() -> float :
                docstring
        
        dfcc(...)
            dfcc() -> float :
                docstring
        
        dfmp2(...)
            dfmp2() -> float :
                docstring
        
        fd_1_0(...)
            fd_1_0( (list)arg1) -> PsiReturnType :
                docstring
        
        fd_freq_0(...)
            fd_freq_0( (list)arg1, (int)arg2) -> PsiReturnType :
                docstring
        
        fd_freq_1(...)
            fd_freq_1( (list)arg1, (int)arg2) -> PsiReturnType :
                docstring
        
        fd_geoms_1_0(...)
            fd_geoms_1_0() -> matrix_vector :
                docstring
        
        fd_geoms_freq_0(...)
            fd_geoms_freq_0( (int)arg1) -> matrix_vector :
                docstring
        
        fd_geoms_freq_1(...)
            fd_geoms_freq_1( (int)arg1) -> matrix_vector :
                docstring
        
        fd_geoms_hessian_0(...)
            fd_geoms_hessian_0() -> matrix_vector :
                docstring
        
        fd_hessian_0(...)
            fd_hessian_0( (list)arg1) -> PsiReturnType :
                docstring
        
        flush_outfile(...)
            flush_outfile() -> None :
                docstring
        
        get_active_molecule(...)
            get_active_molecule() -> Molecule :
                docstring
        
        get_global_option(...)
            get_global_option( (str)arg1) -> object :
                docstring
        
        get_global_option_list(...)
            get_global_option_list() -> list :
                docstring
        
        get_gradient(...)
            get_gradient() -> Matrix :
                docstring
        
        get_input_directory(...)
            get_input_directory() -> str :
                docstring
        
        get_local_option(...)
            get_local_option( (str)arg1, (str)arg2) -> object :
                docstring
        
        get_memory(...)
            get_memory() -> int :
                docstring
        
        get_option(...)
            get_option( (str)arg1) -> object :
                docstring
        
        get_variable(...)
            get_variable( (str)arg1) -> float :
                docstring
        
        has_global_option_changed(...)
            has_global_option_changed( (str)arg1) -> bool :
                docstring
        
        has_local_option_changed(...)
            has_local_option_changed( (str)arg1, (str)arg2) -> bool :
                docstring
        
        has_option_changed(...)
            has_option_changed( (str)arg1) -> bool :
                docstring
        
        libfock(...)
            libfock() -> int :
                docstring
        
        lmp2(...)
            lmp2() -> float :
                docstring
        
        mcscf(...)
            mcscf() -> float :
                docstring
        
        me(...)
            me() -> int :
                docstring
        
        mints(...)
            mints() -> int :
                docstring
        
        mp2(...)
            mp2() -> float :
                docstring
        
        mrcc_generate_input(...)
            mrcc_generate_input( (dict)arg1) -> PsiReturnType :
                docstring
        
        mrcc_load_densities(...)
            mrcc_load_densities( (dict)arg1) -> PsiReturnType :
                docstring
        
        nproc(...)
            nproc() -> int :
                docstring
        
        nthread(...)
            nthread() -> int :
                docstring
        
        nuclear_dipole(...)
            nuclear_dipole( (Molecule)arg1) -> Vector :
                docstring
        
        omp2(...)
            omp2() -> int :
                docstring
        
        opt_clean(...)
            opt_clean() -> None :
                docstring
        
        optking(...)
            optking() -> int :
                docstring
        
        outfile_name(...)
            outfile_name() -> str :
                docstring
        
        plugin(...)
            plugin( (str)arg1) -> int :
                docstring
        
        plugin_close(...)
            plugin_close( (str)arg1) -> None :
                docstring
        
        plugin_close_all(...)
            plugin_close_all() -> None :
                docstring
        
        plugin_load(...)
            plugin_load( (str)arg1) -> int :
                docstring
        
        prepare_options_for_module(...)
            prepare_options_for_module( (str)arg1) -> None :
                docstring
        
        print_global_options(...)
            print_global_options() -> None :
                docstring
        
        print_options(...)
            print_options() -> None :
                docstring
        
        print_out(...)
            print_out( (str)arg1) -> None :
                docstring
        
        print_variables(...)
            print_variables() -> None :
                docstring
        
        psi_top_srcdir(...)
            psi_top_srcdir() -> str :
                docstring
        
        psimrcc(...)
            psimrcc() -> float :
                docstring
        
        reference_wavefunction(...)
            reference_wavefunction() -> Wavefunction :
                docstring
        
        reopen_outfile(...)
            reopen_outfile() -> None :
                docstring
        
        revoke_global_option_changed(...)
            revoke_global_option_changed( (str)arg1) -> None :
                docstring
        
        revoke_local_option_changed(...)
            revoke_local_option_changed( (str)arg1, (str)arg2) -> None :
                docstring
        
        revoke_option_changed(...)
            revoke_option_changed( (str)arg1) -> None :
                docstring
        
        sapt(...)
            sapt() -> float :
                docstring
        
        scf(...)
            scf( (object)arg1, (object)arg2) -> float :
                docstring
            
            scf() -> float :
                docstring
        
        set_active_molecule(...)
            set_active_molecule( (Molecule)arg1) -> None :
                docstring
        
        set_global_option(...)
            set_global_option( (str)arg1, (str)arg2) -> bool :
                docstring
            
            set_global_option( (str)arg1, (float)arg2) -> bool :
                docstring
            
            set_global_option( (str)arg1, (int)arg2) -> bool :
                docstring
            
            set_global_option( (str)arg1, (list)arg2 [, (object)arg3]) -> bool
        
        set_global_option_python(...)
            set_global_option_python( (str)arg1, (object)arg2) -> bool :
                docstring
        
        set_gradient(...)
            set_gradient( (Matrix)arg1) -> None :
                docstring
        
        set_local_option(...)
            set_local_option( (str)arg1, (str)arg2, (str)arg3) -> bool :
                docstring
            
            set_local_option( (str)arg1, (str)arg2, (float)arg3) -> bool :
                docstring
            
            set_local_option( (str)arg1, (str)arg2, (int)arg3) -> bool :
                docstring
            
            set_local_option( (str)arg1, (str)arg2, (list)arg3 [, (object)arg4]) -> bool
        
        set_local_option_python(...)
            set_local_option_python( (str)arg1, (object)arg2) -> None :
                docstring
        
        set_memory(...)
            set_memory( (int)arg1) -> None :
                docstring
        
        set_nthread(...)
            set_nthread( (int)arg1) -> None :
                docstring
        
        set_variable(...)
            set_variable( (str)arg1, (float)arg2) -> None :
                docstring
        
        transqt(...)
            transqt() -> float :
                docstring
        
        transqt2(...)
            transqt2() -> float :
                docstring
        
        version(...)
            version() -> str :
                docstring
    
    DATA
        Ascending = PsiMod.DiagonalizeOrder.Ascending
        Balk = PsiMod.PsiReturnType.Balk
        Descending = PsiMod.DiagonalizeOrder.Descending
        EndLoop = PsiMod.PsiReturnType.EndLoop
        Failure = PsiMod.PsiReturnType.Failure
        Success = PsiMod.PsiReturnType.Success

