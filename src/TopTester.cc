#include "TopTester.h"
#include <iostream>
#include <vector>
#include <string>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/ITuple.h>
#endif

#include <CLHEP/Vector/LorentzVector.h>
#include "JetFitObject.h"
//#include "PConstraint.h"
#include "OPALFitterGSL.h"
#include "NewFitterGSL.h"
#include "TextTracer.h"
#include "NewtonFitterGSL.h"
#include "MassConstraint.h"
#include "TopEventILC.h"

using namespace marlin ;
using namespace std ;
using namespace CLHEP ;


TopTester aTopTester ;


TopTester::TopTester() : Processor("TopTester") {
  
  // modify processor description
  _description = "TopTester does a 7C fit on 6 jet events (Px, Py, Pz, E, M34 = M56 = 80.4 GeV, M134 = M256)" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "ECM" ,
                              "Center-of-Mass Energy in GeV",
                              _ecm,
                              (float)500.);
                              
  registerProcessorParameter( "ntoy" ,
                              "number of toy events",
                              _ntoy,
                              (int)200);
                              
  registerProcessorParameter( "semileptonic" ,
                              "set true if semi-leptonic ttbar events should be used",
                              _semileptonic,
                              (bool)false);
                              
  registerProcessorParameter( "leptonasjet" ,
                              "set true if lepton should be parametrised at JetFitObject",
                              _leptonasjet,
                              (bool)false);
                              
  registerProcessorParameter( "fitter" ,
                              "0 = OPALFitter, 1 = NewFitter, 2 = NewtonFitter",
                              _ifitter,
                              (int)0);
                              
  registerProcessorParameter( "traceall" ,
                              "set true if every event should be traced",
                              _traceall,
                              (bool)false);
                              
  registerProcessorParameter( "ievttrace" ,
                              "number of individual event to be traced",
                              _ievttrace,
                              (int)0);
                              
                              
  topevent = new TopEventILC();   
                           
}


void TopTester::init() { 

  // usually a good idea to
  printParameters() ;
  topevent->leptonic = _semileptonic;
  topevent->leptonasjet = _leptonasjet;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void TopTester::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

void TopTester::processEvent( LCEvent * evt ) { 

    
    message<ERROR>( log() 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      ) ;
                      
  // this gets called for every event 
  // usually the working horse ...

#ifdef MARLIN_USE_AIDA
    
  // define a histogram pointer
  static AIDA::IHistogram1D* hRecTop1Mass;    
  static AIDA::IHistogram1D* hRecTop2Mass;    
  static AIDA::IHistogram1D* hRecTop1MassNoFitOK;    
  static AIDA::IHistogram1D* hRecTop1MassNoFitAll;    
  static AIDA::IHistogram1D* hRecTop2MassNoFitOK;    
  static AIDA::IHistogram1D* hRecTop2MassNoFitAll;    
  static AIDA::IHistogram1D* hRecTopMass;    
  static AIDA::IHistogram1D* hRecTopMassNoFitOK;    
  static AIDA::IHistogram1D* hRecTopMassNoFitAll;    
  static AIDA::IHistogram1D* hRecW1Mass;    
  static AIDA::IHistogram1D* hRecW2Mass;    
  static AIDA::IHistogram1D* hRecW1MassNoFitOK;    
  static AIDA::IHistogram1D* hRecW1MassNoFitAll;    
  static AIDA::IHistogram1D* hRecW2MassNoFitOK;    
  static AIDA::IHistogram1D* hRecW2MassNoFitAll;    
  static AIDA::IHistogram1D* hRecWMass;    
  static AIDA::IHistogram1D* hRecWMassNoFitOK;    
  static AIDA::IHistogram1D* hRecWMassNoFitAll;    
  static AIDA::IHistogram1D* hFitProb;    
  static AIDA::IHistogram1D* hNIt;    
  static AIDA::IHistogram1D* hFitError;    
               
  if( isFirstEvent() ) { 
    
    hRecTopMass = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTopMass", "M_top", 200, 125., 225. ) ; 

    hRecTop1Mass = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTop1Mass", "M_top1", 200, 125., 225. ) ; 

    hRecTop2Mass = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTop2Mass", "M_top2", 200, 125., 225. ) ; 

    hRecTopMassNoFitOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTopMassNoFitOK", "prefit average M_top, fit OK", 200, 125., 225. ) ; 

    hRecTopMassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTopMassNoFitAll", "prefit average M_top, all", 200, 125., 225. ) ; 

    hRecTop1MassNoFitOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTop1MassNoFitOK", "prefit M_top1, fit OK", 200, 125., 225. ) ; 

    hRecTop1MassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTop1MassNoFitAll", "prefit M_top1, all", 200, 125., 225. ) ; 

    hRecTop2MassNoFitOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTop2MassNoFitOK", "prefit M_top2, fit OK", 200, 125., 225. ) ; 

    hRecTop2MassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTop2MassNoFitAll", "prefit M_top2, all", 200, 125., 225. ) ; 

    hRecWMass = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMass", "M_W", 200, 0., 200. ) ; 

    hRecW1Mass = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecW1Mass", "M_W1", 200, 0., 200. ) ; 

    hRecW2Mass = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecW2Mass", "M_W2", 200, 0., 200. ) ; 

    hRecW1MassNoFitOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMass1NoFitOK", "prefit M_W1, fit OK", 200, 0., 200. ) ; 

    hRecW1MassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMass1NoFitAll", "prefit M_W1, all", 200, 0., 200. ) ; 

    hRecW2MassNoFitOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMass2NoFitOK", "prefit M_W2, fit OK", 200, 0., 200. ) ; 

    hRecW2MassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMass2NoFitAll", "prefit M_W2, all", 200, 0., 200. ) ; 

    hRecWMassNoFitOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMassNoFitOK", "prefit average M_W, fit OK", 200, 0., 200. ) ; 

    hRecWMassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMassNoFitAll", "prefit average M_W, all", 200, 0., 200. ) ; 

    hFitProb = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProb", "fit probability", 100, 0., 1. ) ; 

    hNIt = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hNIt", "number of iterations", 201, -0.5, 200.5 ) ; 

    if (_ifitter == 1) {
      hFitError = 
        AIDAProcessor::histogramFactory(this)->
        createHistogram1D( "hFitError", "Error flag", 100, -0.5, 99.5 ) ; 
    }
    else {    
      hFitError = 
        AIDAProcessor::histogramFactory(this)->
        createHistogram1D( "hFitError", "Error flag", 10, -0.5, 9.5 ) ; 
    }    

  }

#endif   
      
   
  for (int ievt = 0; ievt < _ntoy; ievt++) { 
  
     message<MESSAGE>( log()  << "start to process toy event number " << ievt ) ;
     
     double startmassW1 = 0., startmassW2 = 0.;
     double startmasstop1 = 0., startmasstop2 = 0.;
     
     bool debug = false;
     if (ievt == _ievttrace || _traceall) debug = true;
     topevent->setDebug (debug);
          
     topevent->genEvent();
       
     startmassW1 = topevent->getW1Mass();
     startmassW2 = topevent->getW2Mass();
     message<DEBUG>( log()  << "start mass of W 1: " << startmassW1 ) ;
     message<DEBUG>( log()  << "start mass of W 2: " << startmassW2 ) ;
     startmasstop1 = topevent->getTop1Mass();
     startmasstop2 = topevent->getTop2Mass();
     message<DEBUG>( log()  << "start mass of top 1: " << startmasstop1 ) ;
     message<DEBUG>( log()  << "start mass of top 2: " << startmasstop2 ) ;
                     
#ifdef MARLIN_USE_AIDA
     hRecTop1MassNoFitAll->fill( startmasstop1 ) ;
     hRecTop2MassNoFitAll->fill( startmasstop2 ) ;
     hRecTopMassNoFitAll->fill( 0.5*(startmasstop1+startmasstop2) ) ;
     hRecW1MassNoFitAll->fill( startmassW1 ) ;
     hRecW2MassNoFitAll->fill( startmassW2 ) ;
     hRecWMassNoFitAll->fill( 0.5*(startmassW1+startmassW2) ) ;
#endif

     BaseFitter *pfitter;
     if (_ifitter == 1) { 
       pfitter = new NewFitterGSL();
       (dynamic_cast<NewFitterGSL*>(pfitter))->setDebug (debug);
     }  
     else if (_ifitter == 2) { 
       pfitter = new NewtonFitterGSL();
       (dynamic_cast<NewtonFitterGSL*>(pfitter))->setDebug (debug);
     }  
     else {
       // OPALFitter has no method setDebug !
       pfitter = new OPALFitterGSL();
     }
     BaseFitter &fitter = *pfitter;
  
     TextTracer tracer (std::cout);
     if (ievt == _ievttrace || _traceall) fitter.setTracer (tracer);
           
     int ierr = topevent->fitEvent(fitter);
  
     double prob = fitter.getProbability();
     double chi2 = fitter.getChi2();
     int nit = fitter.getIterations();

     message<DEBUG>( log() << "fit probability = " << prob ) ;  
     message<DEBUG>( log() << "fit chi2 = " << chi2  ) ; 
     message<DEBUG>( log() << "error code: " << ierr ) ;
                                  
     message<DEBUG>( log()  << "final mass of W 1: " << topevent->getW1Mass() ) ;
     message<DEBUG>( log()  << "final mass of W 2: " << topevent->getW2Mass() ) ;
     message<DEBUG>( log()  << "final mass of top 1: " << topevent->getTop1Mass() ) ;
     message<DEBUG>( log()  << "final mass of top 2: " << topevent->getTop2Mass() ) ;
                  
#ifdef MARLIN_USE_AIDA
     hFitError->fill( ierr ) ;
     if (ierr == 0) {
       hFitProb->fill( prob ) ;
       hNIt->fill( nit ) ;
       hRecTop1MassNoFitOK->fill( startmasstop1 ) ;
       hRecTop2MassNoFitOK->fill( startmasstop2 ) ;
       hRecTopMassNoFitOK->fill( 0.5*(startmasstop1+startmasstop2) ) ;
       hRecW1MassNoFitOK->fill( startmassW1 ) ;
       hRecW2MassNoFitOK->fill( startmassW2 ) ;
       hRecWMassNoFitOK->fill( 0.5*(startmassW1+startmassW2) ) ;
       hRecTopMass->fill( 0.5*(topevent->getTopMass(1)+topevent->getTopMass(2)) ) ;
       hRecTop1Mass->fill( topevent->getTopMass(1) ) ;
       hRecTop2Mass->fill( topevent->getTopMass(2) ) ;
       hRecW1Mass->fill( topevent->getW1Mass() ) ;
       hRecW2Mass->fill( topevent->getW2Mass() ) ;
       hRecWMass->fill( 0.5*(topevent->getW1Mass()+topevent->getW2Mass()) ) ;
     }
#endif
     if (ierr > 0) message<WARNING>( log() << "FIT ERROR = " << ierr << " in toy event " << ievt ) ;

     _nEvt ++ ;
     
   }   
}



void TopTester::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TopTester::end(){ 

 delete topevent;

}

