#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>

#include "RooDCB.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooMath.h"
#include "TMath.h"

using namespace std;

ClassImp(RooDCB);


//_____________________________________________________________________________
RooDCB::RooDCB(const char *name, const char *title,
		       RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _sigma,
		       RooAbsReal& _alphaHi, RooAbsReal& _nHi,
		       RooAbsReal& _alphaLo, RooAbsReal& _nLo) :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  m0("m0", "M0", this, _m0),
  sigma("sigma", "Sigma", this, _sigma),
  alphaHi("alphaHi", "Alpha_Hi", this, _alphaHi),
  nHi("nHi", "Order_Hi", this, _nHi),
  alphaLo("alphaLo", "Alpha_Lo", this, _alphaLo),
  nLo("nLo", "OrderLo", this, _nLo)
{
}


//_____________________________________________________________________________
RooDCB::RooDCB(const RooDCB& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
  sigma("sigma", this, other.sigma),
  alphaHi("alphaHi", this, other.alphaHi), nHi("nHi", this, other.nHi),
  alphaLo("alphaLo", this, other.alphaLo), nLo("nLo", this, other.nLo)
{
}


Double_t RooDCB::evaluate() const {
  Double_t t = (m-m0) / sigma;

  if (t < -alphaLo) {
    Double_t a = alphaLo/nLo;
    Double_t b = (1 - a * (alphaLo + t));

    return exp(-0.5*alphaLo*alphaLo) / TMath::Power(b, nLo);
  }
  else if (t > alphaHi) {
    Double_t a = alphaHi/nHi;
    Double_t b = (1 - a * (alphaHi - t));

    return exp(-0.5*alphaHi*alphaHi) / TMath::Power(b, nHi);
  }
  else {
    return exp(-0.5*t*t);
  }
}

//_____________________________________________________________________________
Int_t RooDCB::getMaxVal(const RooArgSet& vars) const 
{
  // Advertise that we know the maximum of self for given (m0,alpha,n,sigma)
  RooArgSet dummy ;

  if (matchArgs(vars,dummy,m)) {
    return 1 ;  
  }
  return 0 ;  
}

//_____________________________________________________________________________
Double_t RooDCB::maxVal(Int_t code) const
{
  assert(code==1);

  // The maximum value for given (m0,alpha,n,sigma)
  return 1.0 ;
}
