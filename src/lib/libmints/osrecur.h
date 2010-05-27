#ifndef _psi_src_lib_libmints_osrecur_h
#define _psi_src_lib_libmints_osrecur_h

namespace psi {

/*! \ingroup MINTS
 *  \class ObaraSaikaTwoCenterRecursion
 *  \brief Generic Obara and Saika recursion object.
 */
class ObaraSaikaTwoCenterRecursion
{
    int max_am1_;
    int max_am2_;
    
    double **x_;
    double **y_;
    double **z_;
    
    // No default constructor
    ObaraSaikaTwoCenterRecursion();
    // No assignment operator
    ObaraSaikaTwoCenterRecursion& operator=(const ObaraSaikaTwoCenterRecursion&);
    
public:
    /// Constructor, max_am1 and max_am2 are the max angular momentum on center 1 and 2.
    /// Needed to allocate enough memory.
    ObaraSaikaTwoCenterRecursion(int max_am1, int max_am2);
    ~ObaraSaikaTwoCenterRecursion();
    
    /// Returns the x recursion matrix.
    double **x() const { return x_; }
    /// Returns the y recursion matrix.
    double **y() const { return y_; }
    /// Returns the z recursion matrix.
    double **z() const { return z_; }
    
    /// Computes the recursion matrices for the data provided.
    void compute(double PA[3], double PB[3], double gamma, int am1, int am2);
};

/*! \ingroup MINTS
 *  \class ObaraSaikaTwoCenterMIRecursion
 *  \brief Obara and Saika recursion object for moment integrals. Currently not used by DipoleInt, hopefully soon.
 *  THIS CLASS HAS NOT BEEN TESTED!!!
 */
class ObaraSaikaTwoCenterMIRecursion
{
    int max_am1_;
    int max_am2_;
    int max_m_;
    
    double ***x_;
    double ***y_;
    double ***z_;
    
    // No default constructor
    ObaraSaikaTwoCenterMIRecursion();
    // No assignment operator
    ObaraSaikaTwoCenterMIRecursion& operator=(const ObaraSaikaTwoCenterMIRecursion&);
    
public:
    ObaraSaikaTwoCenterMIRecursion(int max_am1, int max_am2, int max_m);
    ~ObaraSaikaTwoCenterMIRecursion();
    
    double ***x() const { return x_; }
    double ***y() const { return y_; }
    double ***z() const { return z_; }
    void compute(double PA[3], double PB[3], double gamma, int am1, int am2);
};

/*! \ingroup MINTS
 *  \class ObaraSaikaTwoCenterVIRecursion
 *  \brief Obara and Saika recursion object for potential integrals.
 */
class ObaraSaikaTwoCenterVIRecursion
{
protected:
    int max_am1_;
    int max_am2_;
    int size_;
    
    double ***vi_;

    // Forms Fm(U) from A20 (OS 1986)
    void calculate_f(double *F, int n, double t);
    
private:    
    // No default constructor
    ObaraSaikaTwoCenterVIRecursion();
    // No assignment operator
    ObaraSaikaTwoCenterVIRecursion& operator=(const ObaraSaikaTwoCenterVIRecursion&);

public:
    /// Constructor, max_am1 and max_am2 are the max angular momentum on center 1 and 2.
    /// Needed to allocate enough memory.
    ObaraSaikaTwoCenterVIRecursion(int max_am1, int max_am2);
    virtual ~ObaraSaikaTwoCenterVIRecursion();

    /// Returns the potential integral 3D matrix
    double ***vi() const { return vi_; }
    
    /// Computes the potential integral 3D matrix using the data provided.
    virtual void compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2);
};

/*! \ingroup MINTS
 *  \class ObaraSaikaTwoCenterVIDerivRecursion
 *  \brief Obara and Saika recursion object for computing potential derivatives.
 */
class ObaraSaikaTwoCenterVIDerivRecursion : public ObaraSaikaTwoCenterVIRecursion
{
protected:
    double ***vx_;
    double ***vy_;
    double ***vz_;
    
private:
    // No default constructor();
    ObaraSaikaTwoCenterVIDerivRecursion();
    // No assignment operator
    ObaraSaikaTwoCenterVIDerivRecursion& operator=(const ObaraSaikaTwoCenterVIDerivRecursion&);
    
public:
    ObaraSaikaTwoCenterVIDerivRecursion(int max_am1, int max_am2);
    virtual ~ObaraSaikaTwoCenterVIDerivRecursion();
    
    double ***vx() const { return vx_; }
    double ***vy() const { return vy_; }
    double ***vz() const { return vz_; }
    
    virtual void compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2);
};

/*! \ingroup MINTS
 *  \class ObaraSaikaTwoCenterElectricField
 *  \brief Obara and Saika recursion object for computing electric field integrals.
 */
class ObaraSaikaTwoCenterElectricField : public ObaraSaikaTwoCenterVIRecursion
{
protected:
    double ***ex_;
    double ***ey_;
    double ***ez_;
    
private:
    // No default constructor
    ObaraSaikaTwoCenterElectricField();
    // No assignment operator
    ObaraSaikaTwoCenterElectricField& operator=(const ObaraSaikaTwoCenterElectricField&);
    
public:
    ObaraSaikaTwoCenterElectricField(int max_am1, int max_am2);
    virtual ~ObaraSaikaTwoCenterElectricField();
    
    double ***ex() const { return ex_; }
    double ***ey() const { return ey_; }
    double ***ez() const { return ez_; }
    
    virtual void compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2);
};

/*! \ingroup MINTS
 *  \class ObaraSaikaTwoCenterElectricFieldGradient
 *  \brief Obara and Saika recursion object for computing electric field gradient integrals.
 */
class ObaraSaikaTwoCenterElectricFieldGradient : public ObaraSaikaTwoCenterElectricField
{
protected:
    double ***exx_;
    double ***eyy_;
    double ***ezz_;
    double ***exy_;
    double ***exz_;
    double ***eyz_;
    
private:
    // No default constructor
    ObaraSaikaTwoCenterElectricFieldGradient();
    // No assignment operator
    ObaraSaikaTwoCenterElectricFieldGradient& operator=(const ObaraSaikaTwoCenterElectricFieldGradient&);
    
public:
    ObaraSaikaTwoCenterElectricFieldGradient(int max_am1, int max_am2);
    virtual ~ObaraSaikaTwoCenterElectricFieldGradient();
    
    double ***ex() const { return ex_; }
    double ***ey() const { return ey_; }
    double ***ez() const { return ez_; }
    
    virtual void compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2);
};

}

#endif
