/*! \file
 *  \brief Implements class TrackParticleFitObject
 *  TrackParticleFitObject works similiar to JetFitObject, but it uses a 1/pt, eta, phi parametrization for the
 *  leptons, which is e.g. more appropriate for electrons.
 *  Especially the covarianz matrix differs from a common E, theta, phi parametrization.
 */

#include "TrackParticleFitObject.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>

using std::sqrt;
using std::sin;
using std::cos;
using std::cout;
using std::endl;

const double TrackParticleFitObject::omega_pt_conv = 2.99792458e-4; // for mm, Tesla, GeV
const double TrackParticleFitObject::maxpt = 500; // GeV

// parameter scalings, to get in approx range 1-100 { iD0,iPhi0,iOmega,iZ0,iTanL }
// DJeans is not sure if we really need these: set them to 1 for now (in which case they have no effect)
const double TrackParticleFitObject::parfact[NPAR] = {1., 1., 1., 1., 1.};

TrackParticleFitObject::TrackParticleFitObject( const EVENT::Track* trk, double m) {
  double ppar[NPAR];
  ppar[ iD0    ] =  trk->getD0()       ;
  ppar[ iPhi0  ] =  trk->getPhi()      ;
  ppar[ iOmega ] =  trk->getOmega()    ;
  ppar[ iZ0    ] =  trk->getZ0()       ;
  ppar[ iTanL  ] =  trk->getTanLambda();

  double cov[15];
  for (int i=0; i<15; i++) cov[i]=trk->getCovMatrix()[i];

  trackReferencePoint.setValues( trk->getReferencePoint()[0],
                                 trk->getReferencePoint()[1],
                                 trk->getReferencePoint()[2] );

  initialise( ppar , cov, m );
}

TrackParticleFitObject::TrackParticleFitObject( const EVENT::TrackState* trk, double m) {
  double ppar[NPAR];
  ppar[ iD0    ] =  trk->getD0()       ;
  ppar[ iPhi0  ] =  trk->getPhi()      ;
  ppar[ iOmega ] =  trk->getOmega()    ;
  ppar[ iZ0    ] =  trk->getZ0()       ;
  ppar[ iTanL  ] =  trk->getTanLambda();

  double cov[15];
  for (int i=0; i<15; i++) cov[i]=trk->getCovMatrix()[i];

  trackReferencePoint.setValues( trk->getReferencePoint()[0],
                                 trk->getReferencePoint()[1],
                                 trk->getReferencePoint()[2] );

  initialise( ppar , cov, m );
}

TrackParticleFitObject::TrackParticleFitObject( const double* _ppars, const double* _cov, double m, const double* refPt_) {
  assert( int(NPAR) <= int(BaseDefs::MAXPAR) );
  if ( refPt_ ) trackReferencePoint.setValues(refPt_[0],refPt_[1],refPt_[2]);
  else          trackReferencePoint.setValues(0,0,0);
  initialise(_ppars, _cov, m);
}


void TrackParticleFitObject::initialise( const double* _ppars, const double* _cov, double m) {
  assert( int(NPAR) <= int(BaseDefs::MAXPAR) );

  initCov();
  setMass (m);

  //  cout << "parameters: ";
  for (int i=0; i<NPAR; i++) {
    //setParam( i, _ppars[i], true, false ); // measured, un-fixed
    //setMParam( i, _ppars[i] );
    setParam ( i, _ppars[i]/parfact[i], true, false );
    setMParam( i, _ppars[i]/parfact[i] );
    //    cout << std::setw(10) << _ppars[i] << " ";
  }
  //  cout << endl;

  setCov( iD0   , iD0    , _cov[ 0] / (parfact[iD0   ]*parfact[iD0   ]) ); // d0  d0
  setCov( iPhi0 , iD0    , _cov[ 1] / (parfact[iPhi0 ]*parfact[iD0   ]) ); // phi d0
  setCov( iPhi0 , iPhi0  , _cov[ 2] / (parfact[iPhi0 ]*parfact[iPhi0 ]) ); // phi phi
  setCov( iOmega, iD0    , _cov[ 3] / (parfact[iOmega]*parfact[iD0   ]) ); // ome d0
  setCov( iOmega, iPhi0  , _cov[ 4] / (parfact[iOmega]*parfact[iPhi0 ]) ); // ome phi
  setCov( iOmega, iOmega , _cov[ 5] / (parfact[iOmega]*parfact[iOmega]) ); // ome ome
  setCov( iZ0   , iD0    , _cov[ 6] / (parfact[iZ0   ]*parfact[iD0   ]) ); // z0  d0
  setCov( iZ0   , iPhi0  , _cov[ 7] / (parfact[iZ0   ]*parfact[iPhi0 ]) ); // z0  phi
  setCov( iZ0   , iOmega , _cov[ 8] / (parfact[iZ0   ]*parfact[iOmega]) ); // z0  ome
  setCov( iZ0   , iZ0    , _cov[ 9] / (parfact[iZ0   ]*parfact[iZ0   ]) ); // z0  z0
  setCov( iTanL , iD0    , _cov[10] / (parfact[iTanL ]*parfact[iD0   ]) ); // tan d0
  setCov( iTanL , iPhi0  , _cov[11] / (parfact[iTanL ]*parfact[iPhi0 ]) ); // tan phi
  setCov( iTanL , iOmega , _cov[12] / (parfact[iTanL ]*parfact[iOmega]) ); // tan ome
  setCov( iTanL , iZ0    , _cov[13] / (parfact[iTanL ]*parfact[iZ0   ]) ); // tan z0
  setCov( iTanL , iTanL  , _cov[14] / (parfact[iTanL ]*parfact[iTanL ]) ); // tan tan

  // parameter iPhi0 repeats every 2*pi
  paramCycl[iPhi0]=2.*M_PI/parfact[iPhi0 ];

  invalidateCache();

  //  cout << "fourmom = " << getFourMomentum() << endl;

  // cout << "errors:     ";
  // for (int i=0; i<NPAR; i++) {
  //   cout << std::setw(10) << sqrt(getCov(i,i)) << " ";
  // }
  // cout << endl;
  // 
  // cout << "correl matrix:" << endl;
  // for (int i=0; i<NPAR; i++) {
  //   for (int j=0; j<NPAR; j++) {
  //     cout << std::setw(10) << getCov(i,j)/sqrt(getCov(i,i)*getCov(j,j)) << " ";
  //   }
  //   cout << endl;
  // }

  return;
}


// destructor
TrackParticleFitObject::~TrackParticleFitObject() {}



TrackParticleFitObject::TrackParticleFitObject (const TrackParticleFitObject& rhs)
{
  //std::cout << "copying TrackParticleFitObject with name " << rhs.name << std::endl;
  TrackParticleFitObject::assign (rhs);
}

TrackParticleFitObject& TrackParticleFitObject::operator= (const TrackParticleFitObject& rhs) {
  if (this != &rhs) {
    assign (rhs); // calls virtual function assign of derived class
  }
  return *this;
}



TrackParticleFitObject *TrackParticleFitObject::copy() const {
  return new TrackParticleFitObject (*this);
}

TrackParticleFitObject& TrackParticleFitObject::assign (const BaseFitObject& source) {
  if (const TrackParticleFitObject *psource = dynamic_cast<const TrackParticleFitObject *>(&source)) {
    if (psource != this) {
      ParticleFitObject::assign (source);
      // only mutable data members, need not to be copied, if cache is invalid
    }
  }
  else {
    assert (0);
  }
  return *this;
}

const char *TrackParticleFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
  case iD0   : return "d0" ;
  case iPhi0 : return "phi0" ;
  case iOmega: return "omega" ;
  case iZ0   : return "z0" ;
  case iTanL : return "tanLambda" ;
  }
  return "undefined";
}

bool TrackParticleFitObject::updateParams (double p[], int idim) {

  //  cout << "TrackParticleFitObject::updateParams " << getName() << endl;

  invalidateCache();

  double tempPar[NPAR]={0};

  // check that omega is not too small (pt too large)
  double omegaMin = fabs(omega_pt_conv*getBfield()/maxpt);

  bool result=false;

  for (int i=0; i<getNPar(); i++) {
    int iglobal = getGlobalParNum(i);
    if (iglobal>=0) {
      tempPar[i] = p[ iglobal ];

      // check that pt is not unphysically large
      if ( i==iOmega ) {
        double omega = tempPar[iOmega]*parfact[iOmega];
        if ( fabs( omega ) < omegaMin ) {
          int signO = par[i]>0 ? 1 : -1;
          tempPar[i] = signO*omegaMin/parfact[iOmega];
          cout << "TrackParticleFitObject::updateParams INFO: regularising Omega to: " << tempPar[i] << endl;
        }
      }

      // check is there has been a significant parameter update
      if ( pow( tempPar[i] - par[i], 2) >  eps2 * cov[i][i] ) result=true; // check if any have been updated
      p[iglobal]=tempPar[i];  // update the global vars
      par[i]    =tempPar[i];  // update local variables
    }
  }

  return result;
}

double TrackParticleFitObject::getDPx(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  return getMomentumFirstDerivatives(1, ilocal);
}

double TrackParticleFitObject::getDPy(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  return getMomentumFirstDerivatives(2, ilocal);
}

double TrackParticleFitObject::getDPz(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  return getMomentumFirstDerivatives(3, ilocal);
}

double TrackParticleFitObject::getDE(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  return getMomentumFirstDerivatives(0, ilocal);
}

double TrackParticleFitObject::getFirstDerivative( int iMeta, int ilocal , int metaSet ) const {
  // iMeta = intermediate variable (i.e. E,px,py,pz)
  // ilocal = local variable (ptinv, theta, phi)
  // metaSet = which set of intermediate varlables
  if (!cachevalid) updateCache();
  switch ( metaSet ) {
  case BaseDefs::VARBASIS_EPXYZ:
    return getMomentumFirstDerivatives(iMeta, ilocal);
    break;
  case BaseDefs::VARBASIS_TRKNORMAL:
    return getNormalFirstDerivatives(iMeta, ilocal);
    break;
  default:
    assert(0);
  }
  return 0; // should never get here
}

double TrackParticleFitObject::getSecondDerivative( int iMeta, int ilocal , int jlocal , int metaSet ) const {
  if (!cachevalid) updateCache();
  switch ( metaSet ) {
  case BaseDefs::VARBASIS_EPXYZ:
    return getMomentumSecondDerivatives(iMeta, ilocal, jlocal);
    break;
  case BaseDefs::VARBASIS_TRKNORMAL:
    return getNormalSecondDerivatives(iMeta, ilocal, jlocal);
    break;
  default:
    assert(0);
  }
  return 0; // should never get here
}


double TrackParticleFitObject::getChi2 () const {
  if (!cachevalid) updateCache();
  return chi2;
}

void TrackParticleFitObject::updateCache() const {
  //  std::cout << "TrackParticleFitObject::updateCache" << std::endl;

  chi2 = ParticleFitObject::getChi2 ();

  phi0  = getParam(iPhi0 )*parfact[iPhi0 ] ; // rescale to physical units
  omega = getParam(iOmega)*parfact[iOmega] ;
  tanl  = getParam(iTanL )*parfact[iTanL ] ;
  d0    = getParam(iD0   )*parfact[iD0   ] ;
  z0    = getParam(iZ0   )*parfact[iZ0   ] ;

  double aB = omega_pt_conv*getBfield();
  double pt = aB/fabs( omega );
  double p  = pt * sqrt ( 1 + pow( tanl, 2 ) );

  fourMomentum.setValues( sqrt ( pow( p, 2 ) + pow ( mass, 2 ) ) ,
                          pt*cos( phi0 ),
                          pt*sin( phi0 ),
                          pt*tanl );

  updateNormalDerivatives();
  updateMomentumDerivatives();

  //  cout << "FourMomentum = " << fourMomentum << endl;
  //  cout << "Normal vector = " << trackPlaneNormal << " " << trackPlaneNormal.getMag() << endl;

  cachevalid = true;
}

void TrackParticleFitObject::updateMomentumDerivatives() const {
  // -------------------------------
  // the momentum derivatives
  // -------------------------------
  resetMomentumFirstDerivatives();
  resetMomentumSecondDerivatives();

  /*
    pt = aB/|omega|

    p  = pt ( 1 + tanl^2 )

    e = ( p^2 + m^2 )^(1/2)
    px = p cos (phi) ( 1 + tanl^2 )^-1
    py = p sin (phi) ( 1 + tanl^2 )^-1
    pz = p tanl ( 1 + tanl^2 )^-1

    (tanl = t)

    intermediate vars: p, phi, tanl

    omega, phi, tanl, z0, d0

    p  = aB ( 1 + t^2 ) / |omega|

    dp / dt = aB * 2 * t / |omega|
    dp / domega = - sign(omega) aB ( 1 + t^2 ) / |omega|^2

    d2p / dt2 = 2 aB / |omega|
    d2p / domega dt = - sign (omega) * 2 ab t / |omega|^2
    d2p / domega2 = 2 ab ( 1 + t^2 ) / |omega|^3

  */

  // intermediate vars: p=0, phi=1, tanl=2
  const int iP = 0;
  const int iPh= 1;
  const int iT = 2;
  const int nInt = 3;

  // track vars: iD0=0, iPhi0, iOmega, iZ0, iTanL, NPAR

  double aB = omega_pt_conv*getBfield();
  double pt = aB/fabs( omega );
  double p  = pt * sqrt ( 1 + pow( tanl, 2 ) );
  double e = sqrt( p*p + mass*mass );
  double one_tan2 = 1 + pow(tanl,2);
  
  double interFirstDerivs[nInt][NPAR];
  double interSecondDerivs[nInt][NPAR][NPAR];
  for (int i=0; i<nInt; i++) {
    for (int j=0; j<NPAR; j++) {
      interFirstDerivs[i][j]=0;
      for (int k=0; k<NPAR; k++) {
	interSecondDerivs[i][j][k]=0;
      }
    }
  }

  int osign = omega>0 ? +1 : -1 ;

  interFirstDerivs[iP][iOmega] = - osign * aB * one_tan2 / pow(omega,2) ;
  interFirstDerivs[iP][iTanL]  = 2*aB*tanl/fabs(omega);

  interFirstDerivs[iPh][iPhi0] = 1;
  interFirstDerivs[iT][iTanL] = 1;

  //  interSecondDerivs[iP][iTanL ][iTanL ] = 2*aB*(1+pow(tanl,2))/pow(fabs(omega),2);
  interSecondDerivs[iP][iTanL ][iTanL ] = 2*aB/fabs(omega); // DJeans fixed 28May2015
  interSecondDerivs[iP][iOmega][iTanL ] = interSecondDerivs[iP][iTanL][iOmega] = -osign*2*aB*tanl/pow(omega,2);
  interSecondDerivs[iP][iOmega][iOmega] = 2*aB*one_tan2/pow(fabs(omega),3);


  // intermediate vars: p=0, phi=1, tanl=2
  double momentumInterFirstDerivs[4][nInt];
  double momentumInterSecondDerivs[4][nInt][nInt];

  for (int i=0; i<4; i++) {
    for (int j=0; j<nInt; j++) {
      momentumInterFirstDerivs[i][j]=0;
      for (int k=0; k<nInt; k++) {
	momentumInterSecondDerivs[i][j][k]=0;
      }
    }
  }

  /*
    e = ( p^2 + m^2 )^(1/2)
    de / dp = (1/2) 2p ( p^2 + m^2 )^(-1/2) = p ( p^2 + m^2 )^(-1/2)
    de / dphi = 0
    de / dt = 0

    d2e / dp2 = ( p^2 + m^2 )^(-1/2) - p*(1/2)*2p*( p^2 + m^2 )^(-3/2) = ( p^2 + m^2 )^(-1/2) - p^2 ( p^2 + m^2 )^(-3/2)
  */

  momentumInterFirstDerivs[0][iP] = p/e;
  momentumInterSecondDerivs[0][iP][iP] = 1./e - pow(p,2)/pow(e,3);

  /*
    px = p cos (phi) ( 1 + t^2 )^-1
    dpx / dp = cos (phi) ( 1 + t^2 )^-1
    dpx / dphi = -p sin(phi) ( 1 + t^2 )^-1
    dpx / dt = p cos (phi) * -2 t (1+t^2)^-2

    d2px/ dp2 = 0
    d2px/ dphi dp  = -sin (phi) ( 1 + t^2 )^-1
    d2px/ dt   dp  = -2t*cos (phi) ( 1 + t^2 )^-2
    d2px/ dp dphi  = - sin(phi) ( 1 + t^2 )^-1
    d2px/ dphi2 = -p cos(phi) ( 1 + t^2 )^-1
    d2px/ dt dphi  = p sin(phi) *2t * ( 1 + t^2 )^-2
    d2px / dp dt   = - cos (phi) * 2 t (1+t^2)^-2
    d2px / dphi dt = p sin(phi) * 2 t (1+t^2)^-2
    d2px / dt2     = -2 p cos (phi) (   (1+t^2)^-2  - 4 t^2 (1+t^2)^-3 )
  */

  momentumInterFirstDerivs[1][iP] = cos(phi0) / one_tan2;
  momentumInterFirstDerivs[1][iPh] = -p*sin(phi0)/one_tan2;
  momentumInterFirstDerivs[1][iT] = -2*p*tanl*cos(phi0)/pow(one_tan2,2);

  momentumInterSecondDerivs[1][iP ][iPh] = momentumInterSecondDerivs[1][iPh][iP ] = -sin(phi0)/one_tan2;
  momentumInterSecondDerivs[1][iP ][iT ] = momentumInterSecondDerivs[1][iT ][iP ] = -2*tanl*cos(phi0)/pow(one_tan2,2);
  momentumInterSecondDerivs[1][iPh][iPh] = -p*cos(phi0)/one_tan2;
  momentumInterSecondDerivs[1][iP ][iT ] = momentumInterSecondDerivs[1][iT ][iP ] = -2*tanl*cos(phi0)/pow(one_tan2,2);
  momentumInterSecondDerivs[1][iPh][iT ] = momentumInterSecondDerivs[1][iT ][iPh] = 2*tanl*p*sin(phi0)/pow(one_tan2,2);
  momentumInterSecondDerivs[1][iT ][iT ] = -2*p*cos(phi0)*( 1./pow(one_tan2,2) - 4*pow(tanl,2)/pow(one_tan2,3) );

  /*
    py = p sin (phi) ( 1 + t^2 )^-1
    dpy / dp = sin (phi) ( 1 + t^2 )^-1
    dpy / dphi = p cos (phi) ( 1 + t^2 )^-1
    dpy / dt = p sin (phi) * -2 t (1+t^2)^-2

    d2py / dp2 = 0
    d2py / dphi dp = cos(phi) ( 1 + t^2 )^-1
    d2py / dt dp   = -2t sin (phi) ( 1 + t^2 )^-2

    d2py / dp dphi = cos (phi) ( 1 + t^2 )^-1
    d2py / dphi2   = -p sin (phi) ( 1 + t^2 )^-1
    d2py / dt dphi = -2t p cos (phi) ( 1 + t^2 )^-2

    d2py / dp dt = sin (phi) * -2 t (1+t^2)^-2
    d2py / dphi dt = p cos (phi) * -2 t (1+t^2)^-2
    d2py / dt2 = -2 p sin (phi) ( (1+t^2)^-2 - 4 t^2 (1+t^2)^-3 )
  */

  momentumInterFirstDerivs[2][iP ] = sin(phi0) / one_tan2;
  momentumInterFirstDerivs[2][iPh] = p*cos(phi0)/one_tan2;
  momentumInterFirstDerivs[2][iT ] = -2*p*tanl*sin(phi0)/pow(one_tan2,2);

  momentumInterSecondDerivs[2][iP ][iPh] = momentumInterSecondDerivs[1][iPh][iP ] = cos(phi0)/one_tan2;
  momentumInterSecondDerivs[2][iP ][iT ] = momentumInterSecondDerivs[1][iT ][iP ] = -2*tanl*sin(phi0)/pow(one_tan2,2);
  momentumInterSecondDerivs[2][iPh][iPh] = -p*sin(phi0)/one_tan2;
  momentumInterSecondDerivs[2][iP ][iT ] = momentumInterSecondDerivs[1][iT ][iP ] = -2*tanl*sin(phi0)/pow(one_tan2,2);
  momentumInterSecondDerivs[2][iPh][iT ] = momentumInterSecondDerivs[1][iT ][iPh] = -2*tanl*p*cos(phi0)/pow(one_tan2,2);
  momentumInterSecondDerivs[2][iT ][iT ] = -2*p*sin(phi0)*( 1./pow(one_tan2,2) - 4*pow(tanl,2)/pow(one_tan2,3) );

  /*

    pz = p t ( 1 + t^2 )^-1
    dpz / dp = t ( 1 + t^2 )^-1
    dpz / dphi = 0
    dpz / dt = p ( 1 + t^2 )^-1 + p*t*-2t*(1+t^2)-2 = p ( ( 1 + t^2 )^-1 - 2t^2 (1+t^2)^-2 );

    d2pz / dp2 = 0
    d2pz / dphi dp = 0
    d2pz / dt dp = ( 1 + t^2 )^-1 - 2 t^2 ( 1 + t^2 )^-2

    d2pz / dp dt = ( ( 1 + t^2 )^-1 - 2t^2 (1+t^2)^-2 )
    d2pz / dt2 = p ( -2t( 1 + t^2 )^-2 - 4t(1+t^2)^-2 + 8 t^3 (1+t^2)^-3 ) = p ( -6t (1+t^2)^-2 + 8 t^3 (1+t^2)^-3 )


  */

  momentumInterFirstDerivs[3][iP ] = tanl/one_tan2;
  momentumInterFirstDerivs[3][iT ] = p*( 1./one_tan2 - 2*pow(tanl/one_tan2,2) );

  momentumInterSecondDerivs[3][iP ][iT ] = momentumInterSecondDerivs[3][iT ][iP ] = 
    1./one_tan2 - 2.*pow(tanl/one_tan2, 2);
  momentumInterSecondDerivs[3][iT ][iT ] = p*( -6.*tanl/pow(one_tan2,2) + 8.*pow(tanl/one_tan2, 3) );

  // now calculate the total derivatives of Epxpypz wrt trk params using chain rule
  for (int ipe=0; ipe<4; ipe++) {
    for (int ipar=0; ipar<NPAR; ipar++) {
      // the first derivs (chain rule)
      double dd(0);
      for (int j=0; j<nInt; j++) {
	dd+=momentumInterFirstDerivs[ipe][j]*interFirstDerivs[j][ipar];
      }
      setMomentumFirstDerivatives(ipe, ipar, dd);
      // the second derivs
      for (int jpar=0; jpar<NPAR; jpar++) {
	double dd2(0);
	for (int j=0; j<nInt; j++) {
	  dd2+=momentumInterFirstDerivs[ipe][j]*interSecondDerivs[j][ipar][jpar];
	  for (int k=0; k<nInt; k++) {
	    dd2+=momentumInterSecondDerivs[ipe][j][k]*interFirstDerivs[j][ipar]*interFirstDerivs[k][jpar];
	  }
	}
	setMomentumSecondDerivatives(ipe, ipar, jpar, dd2);	
      }
    }
  }

  return;
}


void TrackParticleFitObject::updateNormalDerivatives() const {
  // ------------------------------
  // the derivatives of normal to track-IP plane wrt track parameters
  // ------------------------------
  // this is rather messy
  // we use the chain rule to simplify somewhat
  // with intermediate parameters (a,b,c), proportional to the
  // components of the cross product of the line from IP->PCA and the momentum vector @ PCA

  // as of June 2015, this part is not thoroughly tested. DJeans.

  resetNormalFirstDerivatives();
  resetNormalSecondDerivatives();

  // (x,y,z) is the reference point of the track parameters
  double x = trackReferencePoint.getX();
  double y = trackReferencePoint.getY();
  double z = trackReferencePoint.getZ();

  // vector from ref point -> PCA is ( -d0 sin(phi), d0 cos(phi), z0 )
  // PCA vector: PCA = (x,y,z) + ( -d0 sin(phi), d0 cos(phi), z0 )
  // momentum 3-vector at PCA: MOM = pt*( cos(phi0), sin(phi0), tanl )

  //PCA cross MOM
  //  = pt ( (x-d0 sin(phi), y+d0 cos(phi), z+z0) cross ( cos(phi), sin(phi), tanl ) )
  //  = pt ( a , b , c ) <--- a,b,c are "intermediate" parameters, used in applying chain rule
  //
  //a = (y+d0 cos(phi))*tanl - (z+z0)*sin(phi)
  //b = (z+z0)*cos(phi) - (x-d0 sin(phi))*tanl
  //c = (x-d0 sin(phi))*sin(phi) - (y+d0 cos(phi))*cos(phi)
  //  = x sin(phi) - y cos(phi) - d0 (sin2(phi) + cos2(phi) )
  //  = x sin(phi) - y cos(phi) - d0


  double ABC[3] = { (y+d0*cos(phi0))*tanl - (z+z0)*sin(phi0),
                    (z+z0)*cos(phi0) - (x-d0*sin(phi0))*tanl,
                    x*sin(phi0) - y*cos(phi0) - d0};

  // this is the ThreeVector perpendicular to the plane defined by IP, PCA, and track momentum @ PCA
  trackPlaneNormal.setValues( ABC[0], ABC[1], ABC[2] );
  trackPlaneNormal*=1./trackPlaneNormal.getMag();

  //
  // first and second derivatives of intermediate parameters abc wrt track parameters
  //
  //da/d(d0)   = cos(phi)*tanl
  //da/d(z0)   = -sin(phi)
  //da/d(phi)  = -d0 sin(phi)*tanl - (z+z0)*cos(phi)
  //da/d(tanl) = (y+d0 cos(phi))
  //
  //db/d(d0)   = sin(phi)*tanl
  //bd/d(z0)   = cos(phi)
  //db/d(phi)  = -(z+z0)*sin(phi) + d0 cos(phi)*tanl
  //db/d(tanl) = - (x-d0 sin(phi))
  //
  //dc/d(d0)   = -1
  //bc/d(z0)   = 0
  //dc/d(phi)  = x cos(phi) + y sin(phi)
  //dc/d(tanl) = 0

  double ABCderivs[3][NPAR];
  double ABCsecondderivs[3][NPAR][NPAR];

  ABCderivs[0][iPhi0 ]= -d0*sin(phi0)*tanl - (z+z0)*cos(phi0);
  ABCderivs[0][iOmega]= 0;
  ABCderivs[0][iTanL ]= y + d0*cos(phi0);
  ABCderivs[0][iD0   ]= cos(phi0)*tanl;
  ABCderivs[0][iZ0   ]= -sin(phi0);

  for (int i=0; i<3; i++)
    for (int j=0; j<NPAR; j++)
      for (int k=0; k<NPAR; k++)
        ABCsecondderivs[i][j][k]=0;

  ABCsecondderivs[0][iPhi0 ][iPhi0 ] = -d0*cos(phi0)*tanl + (z+z0)*sin(phi0);
  ABCsecondderivs[0][iPhi0 ][iTanL ] = -d0*sin(phi0);
  ABCsecondderivs[0][iPhi0 ][iD0   ] = -sin(phi0)*tanl;
  ABCsecondderivs[0][iPhi0 ][iZ0   ] = -cos(phi0);

  ABCsecondderivs[0][iTanL ][iPhi0 ] = -d0*sin(phi0);
  ABCsecondderivs[0][iTanL ][iD0   ] = cos(phi0);

  ABCsecondderivs[0][iD0   ][iPhi0 ] = -sin(phi0)*tanl;
  ABCsecondderivs[0][iD0   ][iTanL ] = cos(phi0);

  ABCsecondderivs[0][iZ0   ][iPhi0 ] = -cos(phi0);

  ABCderivs[1][iPhi0 ]= -(z+z0)*sin(phi0) + d0*cos(phi0)*tanl;
  ABCderivs[1][iOmega]= 0;
  ABCderivs[1][iTanL ]= -(x-d0*sin(phi0));
  ABCderivs[1][iD0   ]= sin(phi0)*tanl;
  ABCderivs[1][iZ0   ]= cos(phi0);

  ABCsecondderivs[1][iPhi0 ][iPhi0 ] = -(z+z0)*cos(phi0) - d0*sin(phi0)*tanl;
  ABCsecondderivs[1][iPhi0 ][iTanL ] = d0*cos(phi0);
  ABCsecondderivs[1][iPhi0 ][iD0   ] = cos(phi0)*tanl;
  ABCsecondderivs[1][iPhi0 ][iZ0   ] = -sin(phi0);

  ABCsecondderivs[1][iTanL ][iPhi0 ] = d0*cos(phi0);
  ABCsecondderivs[1][iTanL ][iD0   ] = sin(phi0);

  ABCsecondderivs[1][iD0   ][iPhi0 ] = cos(phi0)*tanl;
  ABCsecondderivs[1][iD0   ][iTanL ] = sin(phi0);

  ABCsecondderivs[1][iZ0   ][iPhi0 ] = -sin(phi0);


  ABCderivs[2][iPhi0 ]= x*cos(phi0) + y*sin(phi0);
  ABCderivs[2][iOmega]= 0;
  ABCderivs[2][iTanL ]= 0;
  ABCderivs[2][iD0   ]= -1;
  ABCderivs[2][iZ0   ]= 0;

  ABCsecondderivs[2][iPhi0 ][iPhi0 ] = -x*sin(phi0) + y*cos(phi0);


  //
  // derivatives of the normal vector N wrt the intermediate parameters abc
  //

  double NderivsABC[3][3];

  double sqabc = sqrt( pow(ABC[0],2) + pow(ABC[1],2) + pow(ABC[2],2) );

  for (int j=0; j<3; j++) { // <---- a,b,c
    for (int i=0; i<3; i++) { // <--- 3-vector
      NderivsABC[j][i] = 0;
      if (i==j) NderivsABC[j][i]+=1./sqabc;
      NderivsABC[j][i]-=ABC[j]*ABC[i]/pow(sqabc,3);
    }
  }

  // the normal's second derivatives wrt ABC
  // a little messy. i think this is ok...
  double NsecondderivsABC[3][3][3];
  for (int i=0; i<3; i++) { // <---- a,b,c
    for (int j=0; j<3; j++) { // <---- a,b,c
      for (int k=0; k<3; k++) { // <--- 3-vector
        NsecondderivsABC[i][j][k]=0;
        NsecondderivsABC[i][j][k]+=3*ABC[i]*ABC[j]*ABC[k]/pow(sqabc,5);
        if ( k==i ) {
          NsecondderivsABC[i][j][k]-=ABC[j]/pow(sqabc,3);
        }
        if ( k==j ) {
          NsecondderivsABC[i][j][k]-=ABC[i]/pow(sqabc,3);
        }
        if ( i==j ) {
          NsecondderivsABC[i][j][k]-=ABC[k]/pow(sqabc,3);
        }
      }
    }
  }


  // now sum over abc to get the derivatives of the normal vector N wrt the track parameters.

  // derivatives of vector wrt to track parameters (chain rule)
  //dN/d(d0) = dN/da da/dd0 + dN/db db/dd0 + dN/dc dc/dd0
  for (int i=0; i<NPAR; i++) { // <-- the object's parameters
    for (int j=0; j<3; j++) {  // <-- the normal's parameters
      double totsum(0);
      for (int k=0; k<3; k++) { // <-- sum over intermediate ABC params
        totsum+=NderivsABC[k][j]*ABCderivs[k][i];
      }
      setNormalFirstDerivatives(j,i,totsum);
    }
  }


  for (int i=0; i<NPAR; i++) { // <-- the object's parameters1
    for (int j=0; j<NPAR; j++) { // <-- the object's parameters2
      for (int m=0; m<3; m++) { // <-- three vector
        double sumtot(0);
        for (int k=0; k<3; k++) { // <-- sum over intermediate ABC params
          sumtot+=NderivsABC[k][m]*ABCsecondderivs[k][i][j];
          for (int l=0; l<3; l++) { // <-- sum over intermediate ABC params
            sumtot += NsecondderivsABC[k][l][m]*ABCderivs[k][i]*ABCderivs[l][j];
          }
        }
        setNormalSecondDerivatives(m, i, j, sumtot);
      }
    }
  }

  return;
}


void TrackParticleFitObject::setNormalFirstDerivatives(int i, int j, double x) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR );
  normalFirstDerivatives[i][j]=x * parfact[j];
}

void TrackParticleFitObject::setNormalSecondDerivatives(int i, int j, int k, double x) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR && k>=0 && k<NPAR);
  normalSecondDerivatives[i][j][k]=x * parfact[j] * parfact[k];
  normalSecondDerivatives[i][k][j]=normalSecondDerivatives[i][j][k];
}

void TrackParticleFitObject::setMomentumFirstDerivatives(int i, int j, double x) const {
  assert ( i>=0 && i<4 && j>=0 && j<NPAR );
  momentumFirstDerivatives[i][j]=x * parfact[j];
}

void TrackParticleFitObject::setMomentumSecondDerivatives(int i, int j, int k, double x) const {
  assert ( i>=0 && i<4 && j>=0 && j<NPAR && k>=0 && k<NPAR);
  momentumSecondDerivatives[i][j][k]=x * parfact[j] * parfact[k];
  momentumSecondDerivatives[i][k][j]=momentumSecondDerivatives[i][j][k];
}

double TrackParticleFitObject::getNormalFirstDerivatives(int i, int j) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR );
  return normalFirstDerivatives[i][j];
}

double TrackParticleFitObject::getNormalSecondDerivatives(int i, int j, int k) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR && k>=0 && k<NPAR);
  return normalSecondDerivatives[i][j][k];
}

double TrackParticleFitObject::getMomentumFirstDerivatives(int i, int j) const {
  assert ( i>=0 && i<4 && j>=0 && j<NPAR );
  return momentumFirstDerivatives[i][j];
}

double TrackParticleFitObject::getMomentumSecondDerivatives(int i, int j, int k) const {
  assert ( i>=0 && i<4 && j>=0 && j<NPAR && k>=0 && k<NPAR);
  return momentumSecondDerivatives[i][j][k];
}

void TrackParticleFitObject::resetMomentumFirstDerivatives() const {
  for (int i=0; i<4; i++)
    for (int j=0; j<NPAR; j++)
      momentumFirstDerivatives[i][j]=0;
  return;
}
void TrackParticleFitObject::resetMomentumSecondDerivatives() const {
  for (int i=0; i<4; i++)
    for (int j=0; j<NPAR; j++)
      for (int k=0; k<NPAR; k++)
        momentumSecondDerivatives[i][j][k]=0;
  return;
}

void TrackParticleFitObject::resetNormalFirstDerivatives() const {
  for (int i=0; i<3; i++)
    for (int j=0; j<NPAR; j++)
      normalFirstDerivatives[i][j]=0;
  return;
}

void TrackParticleFitObject::resetNormalSecondDerivatives() const {
  for (int i=0; i<3; i++)
    for (int j=0; j<NPAR; j++)
      for (int k=0; k<NPAR; k++)
        normalSecondDerivatives[i][j][k]=0;
  return;
}

// B field in Tesla
double TrackParticleFitObject::bfield = 3.5;

double TrackParticleFitObject::setBfield (double bfield_) {
  return bfield = bfield_;
}

double TrackParticleFitObject::getBfield () {
  return bfield;
}
