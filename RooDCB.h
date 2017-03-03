#ifndef ROO_DCB_SHAPE
#define ROO_DCB_SHAPE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooDCB : public RooAbsPdf {
public:
  RooDCB() {} ;
  RooDCB(const char *name, const char *title, RooAbsReal& _m,
	     RooAbsReal& _m0, RooAbsReal& _sigma,
	     RooAbsReal& _alphaHi, RooAbsReal& _nHi,
	     RooAbsReal& _alphaLo, RooAbsReal& _nLo);

  RooDCB(const RooDCB& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooDCB(*this,newname); }

  inline virtual ~RooDCB() { }

  // Optimized accept/reject generator support
  virtual Int_t getMaxVal(const RooArgSet& vars) const ;
  virtual Double_t maxVal(Int_t code) const ;

protected:

  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigma;
  RooRealProxy alphaHi;
  RooRealProxy nHi;
  RooRealProxy alphaLo;
  RooRealProxy nLo;

  Double_t evaluate() const;

private:

  ClassDef(RooDCB,1) // Crystal Ball lineshape PDF
};

#endif
