/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */
#ifndef ATOMSET_H_
#define ATOMSET_H_
#include "Set.h"
#include "LibFragTypes.h"
#include "libmints/vector3.h"
#include <boost/shared_ptr.hpp>
#include "masses.h"
namespace psi {
namespace LibFrag {

///A class for things with Cartesian coordinates
class CartObject {
   private:
      ///The actual copy function
      void Copy(const CartObject& other) {
         this->Carts_=other.Carts_;
      }

   protected:
      Vector3 Carts_;

   public:
      ///Returns the coordinates
      Vector3 Carts() const {
         return Carts_;
      }

      ///Sets carts via pointer
      CartObject(const double *Carts) :
            Carts_(Carts) {
      }

      ///Sets carts via components
      CartObject(const double x, const double y, const double z) :
            Carts_(x, y, z) {
      }

      ///Deep copies Carts_
      CartObject(const CartObject& other) {
         this->Copy(other);
      }

      ///No memory to free up
      virtual ~CartObject() {
      }

      ///Sets this Carts object equal to other
      CartObject& operator=(const CartObject& other) {
         if (this!=&other) this->Copy(other);
         return *this;
      }

      ///Returns the distance between this Carts object and other
      double Distance(const CartObject& other) const {
         return this->Carts_.distance(other.Carts_);
      }
};

///A handy typedef for a boost pointer to a CartObject
typedef boost::shared_ptr<CartObject> SharedCarts;

/** \brief A class that holds data associated with Atoms
 *
 *  An atom is defined by it's chemical identity and mass.  It also
 *  has Cartesian coordinates, the properties of which it gets from
 *  CartObject.
 */
class Atom:public virtual CartObject {
   private:
      ///The atomic number of the atom
      int Z_;
      ///The mass of the atom
      double mass_;

      ///Actual function that copies an atom
      void Copy(const Atom& other);

   public:

      ///Constructor for Cart array
      Atom(const int Z, const double m, const double *Carts);

      ///Constructor if coordinates are not in an array
      Atom(const int Z, const double m, const double x, const double y,
            const double z);

      ///Copy constructor
      Atom(const Atom& other) :
            CartObject(other) {
         this->Copy(other);
      }

      ///Destructor, does nothing
      virtual ~Atom() {
      }

      ///Assignment operator, checks for self assignment
      const Atom& operator=(const Atom& other);

      ///Returns the atomic number of the atom
      int Z() const {
         return Z_;
      }

      ///Returns the atomic symbol of the atom
      std::string AtSym() const {
         return atomic_labels[this->Z()];
      }

      ///Returns the mass of the atom
      double Mass() const {
         return mass_;
      }
};

///A handy typedef for a boost pointer to an Atom object
typedef boost::shared_ptr<Atom> SharedAtom;

/** \brief StandIns are objects that "stand in" for actual atoms
 *
 *  Things like ghosts, caps, and charges replace actual atoms.  They may
 *  possess some or all of the properties of an Atom.  Minimally though
 *  they have Cartesian coordinates, which we get from CartObject
 */
class StandIn:public virtual CartObject {
   private:

      ///Copies AtomIReplace_
      void Copy(const StandIn& other) {
         this->AtomIReplace_=other.AtomIReplace_;
      }

   protected:

      ///The identity of the Atom that this StandIn replaced
      int AtomIReplace_;

   public:

      ///Constructor for explicit components
      StandIn(const int ReplacedAtom, const double x, const double y,
            const double z);

      ///Constructor for array of components
      StandIn(const int ReplacedAtom, const double* Carts);

      ///Copy constructor
      StandIn(const StandIn& other) :
            CartObject(other) {
         this->Copy(other);
      }

      ///No memory to free up
      virtual ~StandIn() {
      }

      ///Assignment operator, checks for self-assignment
      const StandIn& operator=(const StandIn& other);

      ///Returns the replaced atom
      int ReplacedAtom() const {
         return AtomIReplace_;
      }
};

///A handy typdef for a boost pointer to a StandIn
typedef boost::shared_ptr<StandIn> SharedStandIn;

/** \brief Caps take care of severed bonds.
 *
 *  Caps are both StandIns and Atoms. They also need to know the atom they
 *  are bonded to.
 */
class Cap:public StandIn,public Atom {
   private:

      ///Copies the bonded atom
      void Copy(const Cap& other) {
         this->AtomIBonded2_=other.AtomIBonded2_;
      }

   protected:

      ///The atom this cap is bonded to
      int AtomIBonded2_;

   public:

      ///Constructor if we have Carts array
      Cap(const int BondAtom, const int ReplaceAtom, const int Z, const int m,
            const double* Carts);

      ///Constructor for components
      Cap(const int BondAt, const int RepAt, const int Z, const int m,
            const double x, const double y, const double z);

      ///Constructor for copying a cap
      Cap(const Cap& other);

      ///Assignment operator
      const Cap& operator=(const Cap& other);

      ///No memory to free up
      ~Cap() {
      }

      ///Returns the atom this cap is bonded to
      int BondedAtom() const {
         return AtomIBonded2_;
      }

      ///Returns a string that can be fed to python
      std::string print_out() const;

};

///Convenient typedef of a shared pointer to a cap
typedef boost::shared_ptr<Cap> SharedCap;

/** \brief Ghosts are added to fragments for the purposes of BSSE
 *   corrections
 *
 *   We define the mass of a ghost atom to be 0.  Thus the only thing
 *   needed for a ghost is the atom it replaces, it's atomic number, and
 *   it's spatial location.
 */
class Ghost:public StandIn,public Atom {
   public:
      ///Makes a ghost via Cart array
      Ghost(const int ReplacedAtom, const int Z, const double *Carts);

      ///Makes a ghost via Cart comps
      Ghost(const int ReplacedAtom, const int Z, const double x, const double y,
            const double z);

      ///Makes a ghost via a copy
      Ghost(const Ghost& other);

      ///Makes this ghost via assignment
      const Ghost& operator=(const Ghost& other);

      ///Destructor, does nothing
      ~Ghost() {
      }
};

///Another typedef, because I hate typing, this time for Ghosts
typedef boost::shared_ptr<Ghost> SharedGhost;

/** \brief Charges are atomic-centered point charges
 *
 *  Charges possess no mass, or other atomic properties, but they
 *  do replace atoms and (get this) have a charge.
 */
class Charge:public StandIn {
   private:

      ///The actual charge of this Charge
      double q_;

      ///Copies the charge
      void Copy(const Charge& other) {
         this->q_=other.q_;
      }
   public:

      ///Construction with Carts array
      Charge(const double q, const int ReplacedAtom, const double* Carts) :
            StandIn(ReplacedAtom, Carts), CartObject(Carts), q_(q) {
      }

      ///Construction with components
      Charge(const double q, const int ReplacedAtom, const double x,
            const double y, const double z) :
            StandIn(ReplacedAtom, x, y, z), CartObject(x, y, z), q_(q) {
      }

      ///Copy constructor
      Charge(const Charge& other) :
            StandIn(other), CartObject(other) {
         this->Copy(other);
      }

      ///Assignment operator
      const Charge& operator=(const Charge& other) {
         StandIn::operator=(other);
         if (this!=&other) this->Copy(other);
         return *this;
      }

      ///Destructor, does nothing
      ~Charge() {
      }

      ///Returns the charge
      double Q() const {
         return q_;
      }
};

typedef boost::shared_ptr<Charge> SharedCharge;

/** \brief A set for things with mass and Cartesian coordinates
 *
 *   In addition to basic set stuff this class also can compute
 *   some additional properties like center of mass, and the distance
 *   between it and other object's centers of mass.
 *
 *   \param[in] T Type of object stored in the set, must have a Mass()
 *                and Carts() fxn.  T is assumed to be a boost::pointer.
 *                If you want T to be something else pass it in as pointer
 *                or else the dereference calls won't work.
 *
 */
template <typename T>
class CartSet:public Set<T>,public CartObject {
   private:

      ///The mass of this Set
      double Mass_;

      ///Call that actually does the copy
      void Copy(const CartSet<T>& other) {
         this->Mass_=other.Mass_;
      }

      ///Calculates the CoM
      void CalcCoM() {
         Mass_=0.0;
         Carts_*=0.0;
         for (int i=0; i<this->size(); i++) {
            double massi=this->Object((*this)[i])->Mass();
            Mass_+=massi;
            Carts_+=massi*this->Object((*this)[i])->Carts();
         }
         if (Mass_>=0.1) Carts_/=Mass_;
      }

      ///Calculates the moments of inertia, as if the fragment was at the center of mass
      void CalcMoI(SharedMatrix& MoI_) const{
         SharedMatrix temp(new Matrix(3, 3));
         MoI_=temp;
         Vector3 com=CoM();
         for (int i=0; i<this->size(); i++) {
            T atom=this->Object((*this)[i]);
            double mass=atom->Mass(),z=atom->Carts()[2]-com[2];
            double x=atom->Carts()[0]-com[0],y=atom->Carts()[1]-com[1];
            (*MoI_)(0, 0)+=mass*(y*y+z*z);
            (*MoI_)(1, 1)+=mass*(x*x+z*z);
            (*MoI_)(2, 2)+=mass*(x*x+y*y);
            (*MoI_)(0, 1)-=mass*x*y;
            (*MoI_)(0, 2)-=mass*x*z;
            (*MoI_)(1, 2)-=mass*y*z;
         }
         (*MoI_)(1, 0)=(*MoI_)(0, 1);
         (*MoI_)(2, 0)=(*MoI_)(0, 2);
         (*MoI_)(2, 1)=(*MoI_)(1, 2);
      }
   public:
      ///Returns the moments of inertia (stolen from libmints/molecule)
      SharedMatrix MoI() const{
         SharedMatrix MoI_;
         CalcMoI(MoI_);
         return MoI_;
      }
            ///Returns the mass
      double Mass() {
         return Mass_;
      }

      ///Wrapper to avoid confusion
      Vector3 CoM() const{return Carts();}

      ///Calculates the distance between this set's CoM and other's
      double Distance(const CartSet<T>& other) const {
         return this->Carts_.distance(other.Carts_);
      }

      CartSet<T>() :CartObject(0.0,0.0,0.0),Mass_(0.0) {
      }

      ///Construction via copy
      CartSet<T>(const CartSet<T>& other) :
            Set<T>(other),CartObject(other) {
         this->Copy(other);
      }

      ///Assignment operator
      const CartSet<T>& operator=(const CartSet<T>& other) {
         Set<T>::operator=(other);
         CartObject::operator=(other);
         if (this!=&other) this->Copy(other);
         return *this;
      }

      ///Nothing to free up
      virtual ~CartSet<T>() {
      }

      ///Operators for adding elements, need to update CoM
      //@{
      void operator+=(const int newelem) {
         Set<T>::operator+=(newelem);
         CalcCoM();
      }

      Set<T> operator<<(const int newelem) {
         Set<T>::operator<<(newelem);
         CalcCoM();
         return (*(dynamic_cast<Set<T>*>(this)));
      }
      //@}
};

/*A set of atoms, knows basic atom-y stuff like carts and masses
 class AtomSet:public Set{
 private:

 /** \brief Makes this a copy of other
 *
 *   Although it may be counter intuitive at first, this constructor
 *   does not copy all AtomSet members.  It only copies the
 /
 void Copy(const AtomSet& other);

 /** \brief This maps the atom number in the input file to data.
 *
 * Atoms in Elem2Atoms are specified by atom number from the input file
 * so as to prevent problems from arising from reordering the elements
 * of the set, which is necessary because they must always be in
 * numeric order for the intersection and union operations to work.
 /
 std::map<int,Atom> Elem2Atoms;

 ///This set's center of mass (a.u.) with coords given to set
 std::vector<double> CoM;

 ///Calculates the CoM
 void CalcCoM();

 public:

 ///Returns the distance (a.u.) between the CoM of this set and "other"
 double Distance(AtomSet& other);

 ///Basic constructor
 AtomSet():Set(){}

 ///Default destructor
 ~AtomSet(){}

 ///Copy constructor
 AtomSet(const AtomSet& other):Set(other){Copy(other);}

 ///Assignment operator
 const AtomSet& operator=(const AtomSet&other);

 /**For BSSE purposes ghosts attached to this set, made it public for
 calls like Set.Ghosts[i] to distinguish from Set[i], which gives
 the real atom i/
 std::vector<int> Ghosts;

 std::vector<double> Charges;
 std::vector<int> Atom2Charge;


 std::vector<Cap> Caps;
 ///Wrapper function to add ghost atoms
 void AddGhost(const int i){Ghosts.push_back(i);}

 ///Sets the carts of atom i, {x,y,z}
 void AddCarts(const int i, const double x,const double y,
 const double z);

 ///Sets the mass of atom i to m
 void AddMass(const int i,const double m);

 ///Returns the mass of atom i
 double Mass(const int i){return Elem2Atoms[Atoms[i]].mass;}

 ///Returns a vector of the carts
 std::vector<double> Carts(const int i){
 return Elem2Atoms[Atoms[i]].carts;
 }
 };
 */
}
}      //End namespaces

#endif /* ATOMSET_H_ */
