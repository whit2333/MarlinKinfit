/*! \file 
 *  \brief Implements class BaseFitObject
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: BaseFitObject.cc,v $
 * - Revision 1.2  2009/09/02 13:10:57  blist
 * - Added errors for NewtonFitterGSL
 * -
 * - Revision 1.1  2008/02/12 10:19:08  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.8  2008/02/04 17:30:53  blist
 * - NewtonFitter works now!
 * -
 * - Revision 1.7  2008/01/29 17:17:33  blist
 * - implemented setname
 * -
 * - Revision 1.6  2007/09/17 12:50:15  blist
 * - Some parameters reordered
 * -
 * - Revision 1.5  2007/09/13 13:33:06  blist
 * - Print methods return os
 * -
 *
 */ 
 
#include "BaseFitObject.h"
#include <cassert>
#include <cstring>
#include <iostream>
#include <cmath>
using std::isfinite;

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

BaseFitObject::BaseFitObject(): name(0) {
  setName ("???");
  invalidateCache();

  for (int ilocal = 0; ilocal < BaseDefs::MAXPAR; ++ilocal) {
    globalParNum[ilocal] = -1;
    fixed[ilocal] = false;
    for (int jlocal = 0; jlocal < BaseDefs::MAXPAR; ++jlocal) 
      cov[ilocal][jlocal] = 0; 
  }

}

BaseFitObject::BaseFitObject (const BaseFitObject& rhs)
: name(0)
{
  //std::cout << "copying BaseFitObject with name" << rhs.name << std::endl;
  if (rhs.name) setName(rhs.name);
  else setName ("???");
  invalidateCache();
}
BaseFitObject& BaseFitObject::operator= (const BaseFitObject& rhs) {
  if (this != &rhs) {
    if (rhs.name) setName(rhs.name);
    else setName ("???");
  }
  return *this;
}
    
BaseFitObject::~BaseFitObject() {
  //std::cout << "destroying BaseFitObject with name" << name << std::endl;
  delete[] name;
}

//const double BaseFitObject::eps2 = 0.00001;
const double BaseFitObject::eps2 = 0.0001; // changed to 1^-4, then sqrt(eps2) corresponds to 1%

void  BaseFitObject::setName (const char * name_) {
  if (name_ == 0) return;
  size_t l = strlen(name_);
  if (name) delete[] name;
  name = new char[l+1];
  strcpy (name, name_);
}

const char * BaseFitObject::getName () const { 
  return name ? name : "???";
}

int BaseFitObject::getNMeasured() const {
  int nmeasured = 0;
  for (int i = 0; i < getNPar(); ++i) if (isParamMeasured(i) && !isParamFixed(i)) ++nmeasured;
  return nmeasured;
}
int BaseFitObject::getNUnmeasured() const {
  int nunmeasrd = 0;
  for (int i = 0; i < getNPar(); ++i) if (!isParamMeasured(i) && !isParamFixed(i)) ++nunmeasrd;
  return nunmeasrd;
}
int BaseFitObject::getNFree() const {
  int nfree = 0;
  for (int i = 0; i < getNPar(); ++i) if (!isParamFixed(i)) ++nfree;
  return nfree;
}
int BaseFitObject::getNFixed() const {
  int nfixed = 0;
  for (int i = 0; i < getNPar(); ++i) if (isParamFixed(i)) ++nfixed;
  return nfixed;
}
    
std::ostream& BaseFitObject::printParams(std::ostream& os) const {
  os << "(";
  for (int i = 0; i < getNPar(); ++i) {
    if (i>0) os << ", ";
    os << getParam(i);
    if (isParamFixed (i))  os << " fix";
    else if (getError(i)>0) os << " \261 " << getError(i);
  }
  os << ")";  
  return os;
}

bool BaseFitObject::updateParams (double p[], int idim) {
  bool result = false;
  invalidateCache();
  for (int ilocal = 0; ilocal < getNPar(); ++ilocal) {
    if ( !isParamFixed(ilocal) ) { // daniel added this
      int iglobal = getGlobalParNum (ilocal);
      assert (iglobal >= 0 && iglobal < idim);
      //      result = result || setParam (ilocal, p[iglobal]); // daniel thinks this is a BUG !!! if first param is successfully updated, no further ones are
      // result = setParam (ilocal, p[iglobal]) || result; // daniel thinks this would be OK
      // but to makes it clearer, does the following:
      bool thisresult = setParam (ilocal, p[iglobal]);
      result = result || thisresult;
      // if illegal value: read back legal value
      if ( !thisresult ) // daniel added
	p[iglobal] = getParam (ilocal);
    }
  }
  return result;
}  

void BaseFitObject::addToGlobCov(double *globCov, int idim) const {
  //  std::cout << " hi from BaseFitObject::addToGlobCov " << std::endl;
  for (int ilocal = 0; ilocal < getNPar(); ++ilocal) {
    if (!isParamFixed(ilocal) && isParamMeasured(ilocal)) {
      int iglobal = getGlobalParNum (ilocal);
      assert (iglobal >= 0 && iglobal < idim);
      int ioffs = idim*iglobal;
      for (int jlocal = 0; jlocal < getNPar(); ++jlocal) {
        if (!isParamFixed(jlocal) && isParamMeasured(jlocal)) {
          int jglobal = getGlobalParNum (jlocal);
          assert (jglobal >= 0 && jglobal < idim);
          globCov[ioffs+jglobal] += getCov(ilocal,jlocal);
        }
      }
    }
  }
}
                                
bool BaseFitObject::calculateCovInv() const {

  // DANIEL added

  //  std::cout << "hello from BaseFitObject::calculateCovInv()" << std::endl;

  int n = getNPar();

  gsl_matrix *covm = gsl_matrix_alloc (n, n);
  gsl_matrix_set_identity (covm);
  
  for (int i = 0; i < n; ++i) {
    if (isParamMeasured (i)) {
      for (int j = 0; j < n; ++j) {
        if (isParamMeasured (j)) {
	  gsl_matrix_set (covm, i, j, cov[i][j]);
	  // std::cout << "BaseFitObject::calculateCovInv getting from cov " << i << " " << j << " " << cov[i][j] << std::endl;
	}
      }
    }
  }
  gsl_error_handler_t *e = gsl_set_error_handler_off ();
  int result = gsl_linalg_cholesky_decomp (covm);
  if (result == 0) result = gsl_linalg_cholesky_invert (covm);
  gsl_set_error_handler (e);
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      covinv[i][j] = covinv[j][i] = gsl_matrix_get (covm, i, j);
    }
    covinv[i][i] = gsl_matrix_get (covm, i, i);
  }

//  std::cout << "cov matrix:" << std::endl;
//  for (int i = 0; i < n; ++i) {
//    for (int j = 0; j < n; ++j) {
//      std::cout << cov[i][j] << " " ;
//    }
//    std::cout << std::endl;
//  }
//
//  std::cout << "corr matrix:" << std::endl;
//  for (int i = 0; i < n; ++i) {
//    for (int j = 0; j < n; ++j) {
//      std::cout << cov[i][j]/sqrt(cov[i][i]*cov[j][j]) << " " ;
//    }
//    std::cout << std::endl;
//  }
//
//  std::cout << "inverse of cov matrix:" << std::endl;
//  for (int i = 0; i < n; ++i) {
//    for (int j = 0; j < n; ++j) {
//      std::cout << covinv[i][j] << " " ;
//    }
//    std::cout << std::endl;
//  }

  gsl_matrix_free(covm);
  covinvvalid = (result == 0);

  if (!covinvvalid) {
    std::cout << "ERROR, COULD NOT INVERT COV MATR!" << std::endl;

    std::cout << "COV " << std::endl;
    for (int i = 0; i < n; ++i) {
      if (isParamMeasured (i)) {
	for (int j = 0; j < n; ++j) {
	  if (isParamMeasured (j)) {
	    std::cout << cov[i][j] << " ";
	  }
	}
	std::cout << std::endl;
      }
    }


    std::cout << "CORREL " << std::endl;
    for (int i = 0; i < n; ++i) {
      if (isParamMeasured (i)) {
	for (int j = 0; j < n; ++j) {
	  if (isParamMeasured (j)) {
	    std::cout << cov[i][j]/sqrt(cov[i][i]*cov[j][j]) << " ";
	  }
	}
	std::cout << std::endl;
      }
    }

  }

  return covinvvalid;
}

bool BaseFitObject::setParam (int ilocal, double par_, 
                                    bool measured_, bool fixed_) {

  // DANIEL moved to BaseFitObject
  assert (ilocal >= 0 && ilocal < getNPar());
  if (measured[ilocal] != measured_ || fixed[ilocal] != fixed_) invalidateCache();
  measured[ilocal] = measured_;
  fixed[ilocal] = fixed_;
  return setParam (ilocal, par_);
}  

bool BaseFitObject::setParam (int ilocal, double par_ ) {
  // DANIEL moved to BaseFitObject 
  if (!isfinite(par_)) return true;
  assert (ilocal >= 0 && ilocal < getNPar());
  if (par[ilocal] == par_) return false;
  invalidateCache();
  bool result = pow( (par_-par[ilocal]) , 2 ) > eps2*cov[ilocal][ilocal]; 
  par[ilocal] = par_;
  return result;
}
  
bool BaseFitObject::setMParam (int ilocal, double mpar_ ) {
  // DANIEL moved to BaseFitObject
  if (!isfinite(mpar_)) return false;
  assert (ilocal >= 0 && ilocal < getNPar());
  if (mpar[ilocal] == mpar_) return true;
  invalidateCache();
  mpar[ilocal] = mpar_;
  return true;
}

bool BaseFitObject::setError (int ilocal, double err_) {
  // DANIEL moved to BaseFitObject 
  if (!isfinite(err_)) return false;
  assert (ilocal >= 0 && ilocal < getNPar());
  invalidateCache();
  covinvvalid = false;
  cov[ilocal][ilocal] = err_*err_;
  return true;
}

bool BaseFitObject::setCov (int ilocal, int jlocal, double cov_) {
  // DANIEL moved to BaseFitObject 
  if (!isfinite(cov_)) return false;
  assert (ilocal >= 0 && ilocal < getNPar());
  assert (jlocal >= 0 && jlocal < getNPar());
  invalidateCache();
  covinvvalid = false;
  cov[ilocal][jlocal] = cov[jlocal][ilocal] = cov_;
  return true;
}


bool BaseFitObject::fixParam (int ilocal, bool fix) {
  // DANIEL moved to BaseFitObject 
  assert (ilocal >= 0 && ilocal < getNPar());
  return fixed [ilocal] = fix;
}

bool BaseFitObject::setGlobalParNum (int ilocal, int iglobal) {
  // DANIEL moved to BaseFitObject 
  if (ilocal < 0 || ilocal >= getNPar()) return false;
  globalParNum[ilocal] = iglobal;
  return true;
}

int  BaseFitObject::getGlobalParNum(int ilocal) const {
  // DANIEL moved to BaseFitObject 
  if (ilocal < 0 || ilocal >= getNPar()) return -1;
  return globalParNum[ilocal];
}

double BaseFitObject::getParam (int ilocal) const {
  // DANIEL moved to BaseFitObject 
  assert (ilocal >= 0 && ilocal < getNPar());
  return par[ilocal];
}

double BaseFitObject::getMParam (int ilocal) const {
  // DANIEL moved to BaseFitObject 
  assert (ilocal >= 0 && ilocal < getNPar());
  return mpar[ilocal];
}

double BaseFitObject::getError (int ilocal) const {
  // DANIEL moved to BaseFitObject 
  assert (ilocal >= 0 && ilocal < getNPar());
  return std::sqrt(cov[ilocal][ilocal]);
}
double BaseFitObject::getCov (int ilocal, int jlocal) const {
  // DANIEL moved to BaseFitObject 
  assert (ilocal >= 0 && ilocal < getNPar());
  assert (jlocal >= 0 && jlocal < getNPar());
  return cov[ilocal][jlocal];
}
bool BaseFitObject::isParamMeasured (int ilocal) const {
  // DANIEL moved to BaseFitObject 
  assert (ilocal >= 0 && ilocal < getNPar());
  return measured[ilocal];
}

bool BaseFitObject::isParamFixed (int ilocal) const {
  // DANIEL moved to BaseFitObject 
  assert (ilocal >= 0 && ilocal < getNPar());
  return fixed[ilocal];
}

double BaseFitObject::getChi2() const {
  // DANIEL moved to BaseFitObject 

  if (!covinvvalid) calculateCovInv();
  if (!covinvvalid) return -1;
  double chi2 = 0;
  static double resid[BaseDefs::MAXPAR];
  static bool chi2contr[BaseDefs::MAXPAR];
  for (int i = 0; i < getNPar(); ++i) {
    resid[i] = par[i]-mpar[i];

    std::cout << " BaseFitObject::getChi2() " << i << " " << resid[i] << " " << covinv[i][i] << std::endl;

    if (chi2contr[i] = isParamMeasured(i) && !isParamFixed(i)) {
      chi2 += resid[i]*covinv[i][i]*resid[i];
      for (int j = 0; j < i; ++j) {
        if (chi2contr[j]) chi2 += 2*(resid[i])*covinv[i][j]*(resid[j]);
      }
    }
  }
  return chi2;
}

double BaseFitObject::getDChi2DParam(int ilocal) const {
  // DANIEL moved to BaseFitObject 
  assert (ilocal >= 0 && ilocal < getNPar());
  if (isParamFixed(ilocal) || !isParamMeasured(ilocal)) return 0;
  if (!covinvvalid) calculateCovInv();
  if (!covinvvalid) return 0;
  double result = 0;
  for (int jlocal = 0; jlocal < getNPar(); jlocal++) 
    if (!isParamFixed(jlocal) && isParamMeasured(jlocal))
      result += covinv[ilocal][jlocal]*(par[jlocal]-mpar[jlocal]);
  return 2*result;
}

double BaseFitObject::getD2Chi2DParam2(int ilocal, int jlocal) const {
  // DANIEL moved to BaseFitObject 
  assert (ilocal >= 0 && ilocal < getNPar() );
  assert (jlocal >= 0 && jlocal < getNPar() );
  if (isParamFixed(ilocal) || !isParamMeasured(ilocal) || 
      isParamFixed(jlocal) || !isParamMeasured(jlocal))
    return 0;
  if (!covinvvalid) calculateCovInv();
  if (!covinvvalid) return 0;
  return 2*covinv[ilocal][jlocal];
}
       

     
void BaseFitObject::addToGlobalChi2DerMatrix (double *M, int idim) const {
  // DANIEL moved to BaseFitObject 
  if (!covinvvalid) calculateCovInv();
  assert( covinvvalid );
  //  if (!covinvvalid) return;
  for (int ilocal = 0; ilocal < getNPar(); ++ilocal) {
    if (!isParamFixed(ilocal) && isParamMeasured(ilocal)) {
      int iglobal = getGlobalParNum (ilocal);
      assert (iglobal >= 0 && iglobal < idim);
      int ioffs = idim*iglobal;
      for (int jlocal = 0; jlocal < getNPar(); ++jlocal) {
	if (!isParamFixed(jlocal) && isParamMeasured(jlocal)) {
	  int jglobal = getGlobalParNum (jlocal);
	  assert (jglobal >= 0 && jglobal < idim);
 	  M[ioffs+jglobal] += getD2Chi2DParam2(ilocal, jlocal);
	}
      }
    }
  }
}


void BaseFitObject::addToGlobalChi2DerVector (double *y, int idim) const {
  // DANIEL moved to BaseFitObject 
  // this adds the dChi2/dpar piece
  assert (getNPar() <= BaseDefs::MAXPAR);
  if (!covinvvalid) calculateCovInv();
  assert (covinvvalid);
  for (int ilocal = 0; ilocal < getNPar(); ++ilocal) {
    if (!isParamFixed(ilocal) && isParamMeasured(ilocal)) {
      int iglobal = getGlobalParNum (ilocal);
      assert (iglobal>= 0 && iglobal < idim);
      y[iglobal] += getDChi2DParam(ilocal);
    }
  }
}


void BaseFitObject::addToGlobalChi2DerVector (double *y, int idim, 
					      double lambda, double der[], int metaSet ) const {
  // DANIEL moved to BaseFitObject 
  // this adds the lambda * dConst/dpar piece
  if (!cachevalid) updateCache();
  for (int ilocal=0; ilocal<getNPar(); ilocal++) {
    int iglobal = globalParNum[ilocal];
    if ( iglobal>=0 ) {
      for (int j=0; j<BaseDefs::nMetaVars[metaSet]; j++) {
	y[iglobal] += lambda * der[j] * getFirstDerivative( j , ilocal , metaSet );
      }
    }
  }
}


void BaseFitObject::addTo1stDerivatives (double M[], int idim, 
					 double der[], int kglobal, int metaSet) const {
  // DANIEL moved to BaseFitObject 
  if (!cachevalid) updateCache();
  for (int ilocal=0; ilocal<getNPar(); ilocal++) {
    int iglobal = globalParNum[ilocal];
    if (iglobal>=0) {
      for (int j=0; j<BaseDefs::nMetaVars[metaSet]; j++) {
	double x = der[j] * getFirstDerivative( j, ilocal , metaSet);
	M[idim*kglobal + iglobal] += x;
	M[idim*iglobal + kglobal] += x;
      }
    }
  }
  return;
}


    
void   BaseFitObject::addTo2ndDerivatives (double der2[], int idim, 
					   double factor[], int metaSet
					   //double efact, double pxfact, 
					   //double pyfact, double pzfact
					   ) const {

  // DANIEL moved to BaseFitObject 
  if (!cachevalid) updateCache();
  for ( int ilocal=0; ilocal<getNPar(); ilocal++) {
    int iglobal = getGlobalParNum(ilocal);
    if ( iglobal<0 ) continue;
    for ( int jlocal=ilocal; jlocal<getNPar(); jlocal++) {
      int jglobal = getGlobalParNum(jlocal);
      if ( jglobal<0 ) continue;
      double sum(0);
      for ( int imeta=0; imeta<BaseDefs::nMetaVars[metaSet]; imeta++) {
	sum+=factor[imeta]*getSecondDerivative( imeta, ilocal , jlocal , metaSet );
      }
      der2[idim*iglobal+jglobal] += sum;
      if ( iglobal!=jglobal ) der2[idim*jglobal+iglobal] += sum;
    }
  }
  return;
}

    
void   BaseFitObject::addTo2ndDerivatives (double M[], int idim,  double lambda, double der[], int metaSet) const {
  // DANIEL moved to BaseFitObject 
  double factor[BaseDefs::MAXINTERVARS];
  for (int i=0; i<BaseDefs::nMetaVars[metaSet]; i++) factor[i]=lambda*der[i];
  addTo2ndDerivatives (M, idim, factor, metaSet );
  return;
}

void BaseFitObject::addToDerivatives (double der[], int idim, 
				      double factor[], int metaSet
				      ) const {
  // DANIEL moved to BaseFitObject 
  if (!cachevalid) updateCache();
  for (int ilocal=0; ilocal<getNPar(); ilocal++) {
    int iglobal = globalParNum[ilocal];
    if ( iglobal >= 0 ) {
      double der_sum(0);
      for ( int iInter=0; iInter<BaseDefs::nMetaVars[metaSet]; iInter++) {
	der_sum += factor[iInter]*getFirstDerivative( iInter , ilocal , metaSet );
      }
      der[iglobal] += der_sum;
    }
  }
  return;
}


void BaseFitObject::initCov() const {
  // DANIEL moved to BaseFitObject 
  for (int i = 0; i < getNPar(); ++i) {
    for (int j = 0; j < getNPar(); ++j) {
      cov[i][j] = static_cast<double>(i == j);
    }
  }    
}


double BaseFitObject::getError2 (double der[], int metaSet) const {
  // DANIEL moved to BaseFitObject 
  if (!cachevalid) updateCache();
  double totError(0);
  for (int i=0; i<BaseDefs::nMetaVars[metaSet]; i++) {
    for (int j=0; j<BaseDefs::nMetaVars[metaSet]; j++) {

      double cov_i_j=0; // this will become the covariance of intermediate variables i and j
      for (int k=0; k<getNPar(); k++) {
	for (int l=0; l<getNPar(); l++) {
	  cov_i_j += getFirstDerivative( i , k , metaSet ) * cov[k][l] * getFirstDerivative( j , l , metaSet );
	}
      }
      totError+=der[i]*der[j]*cov_i_j;

    }
  }

  return totError;
}
