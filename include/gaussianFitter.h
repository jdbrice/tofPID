#ifndef GAUSSIAN_FITTER_H
#define GAUSSIAN_FITTER_H

#include "allroot.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include "RooPlot.h"

class gaussianFitter
{
public:
	gaussianFitter( TH1* dataHist );
	~gaussianFitter();


	
	/**
	 * Roo Variables
	 */
	RooRealVar *xObs, * mean, *sigma, *amp;
	RooGaussian* gaussModel;
};




#endif