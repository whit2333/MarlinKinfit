/*! \file
 *  \brief Declares class LeptonFitObject
 *
 *
 */

#ifndef __LEPTONFITOBJECT_H
#define __LEPTONFITOBJECT_H

#include "ParticleFitObject.h"

// Class LeptonFitObject
/// Class for leptons with (1/pt, eta, phi) in kinematic fits

class LeptonFitObject : public ParticleFitObject {
  public:

    LeptonFitObject(double ptinv, double theta, double phi,
                  double Dptinv, double Dtheta, double Dphi,
                  double m = 0);
                 
    /// Copy constructor
    LeptonFitObject (const LeptonFitObject& rhs              ///< right hand side
                   );
    /// Assignment               
    LeptonFitObject& operator= (const LeptonFitObject& rhs   ///< right hand side
                             );

    virtual ~LeptonFitObject();

    /// Return a new copy of itself
    virtual LeptonFitObject *copy() const;

    /// Assign from anther object, if of same type
    virtual LeptonFitObject& assign (const BaseFitObject& source   ///< The source object
                                   );

    /// Get name of parameter ilocal
    virtual const char *getParamName (int ilocal     ///< Local parameter number
                                      ) const;

    /// Read values from global vector, readjust vector; return: significant change
    virtual bool   updateParams (double p[],   ///< The parameter vector
                                 int idim      ///< Length of the vector
                                 );

    // these depend on actual parametrisation!
    virtual double getDPx(int ilocal) const;
    virtual double getDPy(int ilocal) const;
    virtual double getDPz(int ilocal) const;
    virtual double getDE(int ilocal) const;

    /// Get chi squared from measured and fitted parameters
    //  virtual double getChi2() const;

    double getFirstDerivative( int iMeta, int ilocal , int metaSet ) const; // derivative of intermediate variable iMeta wrt local parameter ilocal
    double getSecondDerivative( int iMeta, int ilocal , int jlocal, int metaSet ) const; // derivative of intermediate variable iMeta wrt local parameter ilocal

    virtual int getNPar() const {return NPAR;}

  protected:

    //  virtual void initCov();

    void updateCache() const;

    //    mutable bool cachevalid;

    mutable double ctheta, stheta, stheta2, cphi, sphi, cottheta,
      p2, p, e, e2, pt, pt2, pt3, px, py, pz, dpdptinv, dpdtheta, dptdptinv,
      dpxdptinv, dpydptinv, dpzdptinv, dpxdtheta, dpydtheta, dpzdtheta, dpxdphi, dpydphi,
      chi2, dEdptinv, dEdtheta, dEdp;

    static bool adjustPtinvThetaPhi (double& m, double &ptinv, double& theta, double& phi);

    enum {NPAR=3};

};



#endif // __LEPTONFITOBJECT_H

