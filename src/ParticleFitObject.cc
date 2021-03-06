/*! \file 
 *  \brief Implements class ParticleFitObject
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: ParticleFitObject.cc,v $
 * - Revision 1.7  2009/02/17 12:46:35  blist
 * - Improved version of NewtonFitterGSL, JetFitObject changed
 * -
 * - Revision 1.6  2009/02/11 15:33:49  mbeckman
 * - Bug fixes: mass initialization in ParticleFitObject, parameter handling in PhotonFitObjectPxyg
 * -
 * - Revision 1.5  2008/11/23 17:53:41  mbeckman
 * - Fixed minor bug in ParticleFitObject.cc
 * -
 * - Revision 1.4  2008/10/17 13:17:17  blist
 * - Avoid variable-size arrays
 * -
 * - Revision 1.3  2008/10/16 08:13:44  blist
 * - New versions of OPALfitter and Newtonfitter using GSL
 * -
 * - Revision 1.2  2008/09/26 09:58:11  boehmej
 * - removed ~100 semicolons after } at end of function implementation :)
 * -
 * - Revision 1.1  2008/02/12 10:19:09  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.7  2008/02/04 17:30:54  blist
 * - NewtonFitter works now!
 * -
 * - Revision 1.6  2008/01/30 21:48:03  blist
 * - Newton Fitter still doesnt work :-(
 * -
 * - Revision 1.5  2008/01/30 09:14:54  blist
 * - Preparations for NewtonFitter
 * -
 * - Revision 1.4  2007/09/17 12:50:15  blist
 * - Some parameters reordered
 * -
 * - Revision 1.3  2007/09/13 13:33:06  blist
 * - Print methods return os
 * -
 * - Revision 1.2  2007/09/13 08:09:51  blist
 * - Updated 2nd derivatives for px,py,pz,E constraints, improved header documentation
 * -
 *
 */ 
 
#include "ParticleFitObject.h"
#include "FourVector.h"

#include <iostream>
#undef NDEBUG
#include <cassert>
#include <cmath>
using std::isfinite;
using std::cout; 
using std::endl;

// #include <TMatrixDSym.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>


ParticleFitObject::ParticleFitObject()
  : mass (0), fourMomentum( FourVector(0,0,0,0) )
{
  for (int i=0; i<BaseDefs::MAXPAR; i++)
    paramCycl[i]=-1;
}

ParticleFitObject::~ParticleFitObject()
{}

ParticleFitObject::ParticleFitObject (const ParticleFitObject& rhs)
  : mass(0), fourMomentum( FourVector(0,0,0,0) )
{
  //std::cout << "copying ParticleFitObject with name" << rhs.name << std::endl;
  ParticleFitObject::assign (rhs);
}

ParticleFitObject& ParticleFitObject::operator= (const ParticleFitObject& rhs) {
  if (this != &rhs) {
    assign (rhs); // calls virtual function assign of derived class
  }
  return *this;
}

bool ParticleFitObject::setMass (double mass_) {
  if (!isfinite(mass_)) return false;
  if (mass == mass_) return true;
  invalidateCache();
  mass = std::abs(mass_);
  return true;
}

    
ParticleFitObject& ParticleFitObject::assign (const BaseFitObject& source) {
  if (const ParticleFitObject *psource = dynamic_cast<const ParticleFitObject *>(&source)) {
    if (psource != this) {
      BaseFitObject::assign (source);
      mass = psource->mass;
      for (int i =0; i < BaseDefs::MAXPAR; ++i)
        paramCycl[i] = psource->paramCycl[i];
    }
  }
  else {
    assert (0);
  }
  return *this;
}


double ParticleFitObject::getMass () const {
  return mass;
}
    
std::ostream&  ParticleFitObject::print4Vector(std::ostream& os) const {
  os << "[" << getE() << ", " << getPx() << ", " 
      << getPy() << ", "  << getPz() << "]";
  return os;
}

FourVector ParticleFitObject::getFourMomentum() const {
  if (!cachevalid) updateCache();
  return fourMomentum;
}
double ParticleFitObject::getE()   const {
  return getFourMomentum().getE();
}
double ParticleFitObject::getPx()  const {
  return getFourMomentum().getPx();
}
double ParticleFitObject::getPy()  const {
  return getFourMomentum().getPy();
}
double ParticleFitObject::getPz()  const {
  return getFourMomentum().getPz();
}
double ParticleFitObject::getP()   const {
  return getFourMomentum().getP();
}
double ParticleFitObject::getP2()  const {
  return getFourMomentum().getP2();
}
double ParticleFitObject::getPt()  const {
  return getFourMomentum().getPt();
}
double ParticleFitObject::getPt2() const {
  return getFourMomentum().getPt2();
}


std::ostream&  ParticleFitObject::print (std::ostream& os) const {

  if (!cachevalid) updateCache();

  printParams(os);
  os << " => ";
  print4Vector(os);
  os << std::endl;
  return os;
}

void ParticleFitObject::addToGlobalChi2DerVectorNum (double *y, int idim, double eps)  {
  for (int ilocal = 0; ilocal < getNPar(); ++ilocal) {
    int iglobal = getGlobalParNum(ilocal);
    y[iglobal] += num1stDerivative (ilocal, eps);
  }
}


void ParticleFitObject::addToGlobalChi2DerMatrixNum (double *M, int idim, double eps) {
  for (int ilocal1 = 0; ilocal1 < getNPar(); ++ilocal1) {
    int iglobal1 = getGlobalParNum (ilocal1);
    for (int ilocal2 = ilocal1; ilocal2 < getNPar(); ++ilocal2) {
      int iglobal2 = getGlobalParNum (ilocal2);
      M[idim*iglobal1 + iglobal2]+= num2ndDerivative (ilocal1, eps, ilocal2, eps);
    }
  }
}

void ParticleFitObject::getDerivatives (double der[], int idim) const {
  for (int ilocal = 0; ilocal < getNPar(); ++ilocal) {
    assert (ilocal < idim);
    der [4*ilocal]   = getDE  (ilocal);
    der [4*ilocal+1] = getDPx (ilocal);
    der [4*ilocal+2] = getDPy (ilocal);
    der [4*ilocal+3] = getDPz (ilocal);
  }
}

void ParticleFitObject::test1stDerivatives () {
  cout << "ParticleFitObject::test1stDerivatives, object " << getName() << "\n";
  double ycalc[100],ynum[100];
  for (int i = 0; i < 100; ++i) ycalc[i]=ynum[i]=0;
  addToGlobalChi2DerVector (ycalc, 100);
  double eps = 0.00001;
  addToGlobalChi2DerVectorNum (ynum, 100, eps);
  for (int ilocal = 0; ilocal < getNPar(); ++ilocal) {
    int iglobal = getGlobalParNum(ilocal);
    double calc = ycalc[iglobal];
    double num = ynum[iglobal];
    cout << "fo: " << getName() << " par " << ilocal << "/" 
         << iglobal << " ("<< getParamName(ilocal)
         << ") calc: " << calc << " - num: " << num << " = " << calc-num
         << endl;
  }
}

void ParticleFitObject::test2ndDerivatives () {
  cout << "ParticleFitObject::test2ndDerivatives, object " << getName() << "\n";
  const int idim=100;
  double *Mnum = new double[idim*idim];
  double *Mcalc = new double[idim*idim];
  for (int i = 0; i < idim*idim; ++i) Mnum[i]=Mcalc[i]=0;
  addToGlobalChi2DerMatrix (Mcalc, idim);
  double eps = 0.0001;
  cout << "eps=" << eps << endl;
  addToGlobalChi2DerMatrixNum (Mnum, idim, eps);
  for (int ilocal1 = 0; ilocal1 < getNPar(); ++ilocal1) {
    int iglobal1 = getGlobalParNum (ilocal1);
    for (int ilocal2 = ilocal1; ilocal2 < getNPar(); ++ilocal2) {
      int iglobal2 = getGlobalParNum (ilocal2);
      double calc = Mcalc[idim*iglobal1 + iglobal2];
      double num = Mnum[idim*iglobal1 + iglobal2];
      cout << "fo: " << getName() << " par " << ilocal1 << "/" 
           << iglobal1 << " ("<< getParamName(ilocal1)
           << "), par " << ilocal2 << "/" 
           << iglobal2 << " ("<< getParamName(ilocal2)
           << ") calc: " << calc << " - num: " << num << " = " << calc-num
           << endl;
    }
  }
  delete[] Mnum;
  delete[] Mcalc;
}

double ParticleFitObject::num1stDerivative (int ilocal, double eps) {
    double save = getParam (ilocal);
    setParam (ilocal, save+eps);
    double v1 = getChi2();
    setParam (ilocal, save-eps);
    double v2 = getChi2();
    double result = (v1-v2)/(2*eps);
    setParam (ilocal, save);
    return result;
}

double ParticleFitObject::num2ndDerivative (int ilocal1, double eeps1,
                                            int ilocal2, double eeps2) {
  double result;

  if (ilocal1 == ilocal2) {
    double save = getParam (ilocal1);
    double v0 = getChi2();
    setParam (ilocal1, save+eeps1);
    double v1 = getChi2();
    setParam (ilocal1, save-eeps1);
    double v2 = getChi2();
    result = (v1+v2-2*v0)/(eeps1*eeps1);
    setParam (ilocal1, save);
  }
  else {
    double save1 = getParam (ilocal1);
    double save2 = getParam (ilocal2);
    setParam (ilocal1, save1+eeps1);
    setParam (ilocal2, save2+eeps2);
    double v11 = getChi2();
    setParam (ilocal2, save2-eeps2);
    double v12 = getChi2();
    setParam (ilocal1, save1-eeps1);
    double v22 = getChi2();
    setParam (ilocal2, save2+eeps2);
    double v21 = getChi2();
    result = (v11+v22-v12-v21)/(4*eeps1*eeps2);
    setParam (ilocal1, save1);
    setParam (ilocal2, save2);
  }
  return result;
}


double ParticleFitObject::getChi2 () const {
  // reimplemented here to take account of cyclical variables e.g azimuthal angle phi - DJeans

  //cout << "hello from ParticleFitObject::getChi2 () " << endl;
  
  if (!covinvvalid) calculateCovInv();
  if (!covinvvalid) return -1;

  double resid[BaseDefs::MAXPAR]={0};
  for (int i=0; i<getNPar(); i++) {
    //    cout << i << " " << isParamMeasured(i) << " " << isParamFixed(i) << endl;
    if ( isParamMeasured(i) && !isParamFixed(i) ) {
      resid[i] = par[i] - mpar[i];
      //cout << "  xxx  " << i << " " << par[i] << " " << mpar[i] << " " << resid[i] << endl;
      if ( paramCycl[i]>0 ) {
	resid[i]=fmod(resid[i],paramCycl[i]);
	if(resid[i] >  paramCycl[i]/2 ) resid[i]-=paramCycl[i];
	if(resid[i] < -paramCycl[i]/2 ) resid[i]+=paramCycl[i];
      }
    }
  }

  //cout << " ParticleFitObject::getChi2  " << endl;
  //for (int i=0; i<getNPar(); i++) {
  // cout << " --- " << i << " " << resid[i] << " " << covinv[i][i] << endl;
  //}

  double chi2 = 0;
  for (int i=0; i<getNPar(); i++) {
    if ( isParamMeasured(i) && !isParamFixed(i) ) {
      for (int j=0; j<getNPar(); j++) {
	if ( isParamMeasured(j) && !isParamFixed(j) ) {
	  chi2+=resid[i]*covinv[i][j]*resid[j];
	  //	  cout << getName () << " === " << i << " " << j << " : " << 
	  // resid[i] << "*" << covinv[i][j] << "*" << resid[j] << " = " << resid[i]*covinv[i][j]*resid[j] << " , sum " << chi2 << endl;
	}
      }
    }
  }

  //  cout << "ParticleFitObject::getChi2 () chi2 = " << chi2 << endl;

  return chi2;
}
