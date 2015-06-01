/*! \file
 *  \brief Declares class LeptonFitObject
 *
 *
 */

#ifndef __TRACKPARTICLEFITOBJECT_H
#define __TRACKPARTICLEFITOBJECT_H

#include "ParticleFitObject.h"

// Class TrackParticleFitObject
/// Class for lcio tracks

#include "EVENT/TrackState.h"
#include "EVENT/Track.h"

class TrackParticleFitObject : public ParticleFitObject {
 public:

  TrackParticleFitObject( const EVENT::Track*      trk, double m);
  TrackParticleFitObject( const EVENT::TrackState* trk, double m);
  TrackParticleFitObject( const double* _ppars, const double* _cov, double m, const double* refPt_=0);

  virtual ~TrackParticleFitObject();

  TrackParticleFitObject (const TrackParticleFitObject& rhs
			  );

  TrackParticleFitObject& operator= (const TrackParticleFitObject& rhs
				     );

  /// Return a new copy of itself
  virtual TrackParticleFitObject *copy() const;

  /// Assign from anther object, if of same type
  virtual TrackParticleFitObject& assign (const BaseFitObject& source   ///< The source object
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

  virtual double getFirstDerivative ( int iMeta, int ilocal , int metaSet ) const; // derivative of intermediate variable iMeta wrt local parameter ilocal
  virtual double getSecondDerivative( int iMeta, int ilocal , int jlocal, int metaSet ) const; // derivative of intermediate variable iMeta wrt local parameter ilocal

  virtual double getChi2 () const;

  virtual int getNPar() const {return NPAR;}

  virtual ThreeVector getTrackPlaneNormal() const {return trackPlaneNormal;}

  /// Set the B field for all tracks
  static double setBfield (double bfield_             ///< New Value of B field (in Tesla)
			   );
  
  /// Get the B field for all tracks (in Tesla)
  inline static double getBfield ();
  
  /// Global B field in Tesla(!)
  static double bfield;

  static const double omega_pt_conv;
  static const double maxpt;

  enum {iD0=0, iPhi0, iOmega, iZ0, iTanL, NPAR};
  
 protected:

  static const double parfact[NPAR];

  virtual void initialise( const double* _pars, const double* _cov, double m);

  void updateCache() const;
  void updateMomentumDerivatives() const;
  void updateNormalDerivatives() const;

  mutable ThreeVector trackReferencePoint;
  mutable ThreeVector trackPlaneNormal;

  mutable double momentumFirstDerivatives[4][NPAR];
  mutable double momentumSecondDerivatives[4][NPAR][NPAR];

  mutable double normalFirstDerivatives[3][NPAR];
  mutable double normalSecondDerivatives[3][NPAR][NPAR];

  mutable double phi0 ;
  mutable double omega;
  mutable double tanl ;
  mutable double d0   ;
  mutable double z0   ;

  mutable double chi2;

  void   resetMomentumFirstDerivatives() const;
  void   resetMomentumSecondDerivatives() const;

  void   resetNormalFirstDerivatives() const;
  void   resetNormalSecondDerivatives() const;

  void   setMomentumFirstDerivatives(int i, int j, double x) const;
  void   setMomentumSecondDerivatives(int i, int j, int k, double x) const;

  void   setNormalFirstDerivatives(int i, int j, double x) const;
  void   setNormalSecondDerivatives(int i, int j, int k, double x) const;

  virtual double getMomentumFirstDerivatives(int i, int j) const;
  virtual double getMomentumSecondDerivatives(int i, int j, int k) const;

  virtual double getNormalFirstDerivatives(int i, int j) const;
  virtual double getNormalSecondDerivatives(int i, int j, int k) const;

};



#endif // __TRACKPARTICLEFITOBJECT_H

