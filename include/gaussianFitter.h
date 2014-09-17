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
	gaussianFitter( ){
		RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
		useROI = false;

		x = 0;
		mean = 0;
		sigma = 0;
	}
	~gaussianFitter(){

		if ( x ){
			delete x;
			x = 0;
		}
		if ( sigma ){
			delete sigma;
			sigma = 0;
		}
		if ( mean ){
			delete mean;
			mean = 0;
		}
		cout << "[gaussianFitter." << __FUNCTION__ << "]" << endl;
	}

	void makeX( double x1, double x2, string name = "x" ) {
		if ( x ){
			delete x;
			x = 0;
		}
		x = new RooRealVar( name.c_str(), name.c_str(), x1, x2 );
	}

	void makeSigma( double val, double x1, double x2 ){
		if ( sigma ){
			delete sigma;
			sigma = 0;
		}
		sigma = new RooRealVar( "sigma", "sigma", val, x1, x2 );	
	}
	void makeMean( double val, double x1, double x2 ){
		if ( mean ){
			delete mean;
			mean = 0;
		}
		mean = new RooRealVar( "mean", "mean", val, x1, x2 );	
	}

	void setROI( double x1, double x2 ){
		if ( x ){
			useROI = true;
			x->setRange( "roi", x1, x2 );
		}
	}

	void fitTo( TH1D * h, int nFits = 1, string name = "Gauss1D" ){
		using namespace RooFit;

		rdh = new RooDataHist( "data", "data", RooArgSet( *x ), h );

		if ( gaussModel ){
			delete gaussModel;
			gaussModel = 0;
		}

		gaussModel = new RooGaussian( name.c_str(), name.c_str(), *x, *mean, *sigma );

		for ( int iFit = 0; iFit < nFits; iFit ++ ){
			if ( useROI )
				gaussModel->fitTo( *rdh, Range( "roi" ), PrintLevel(-1) );
			else 
				gaussModel->fitTo( *rdh, PrintLevel(-1) );
		}

	}

	RooPlot * drawFit( string title = "") {
		using namespace RooFit;
		RooPlot * frame = x->frame( Title( title.c_str() ) );
		rdh->plotOn( frame );
		gaussModel->plotOn( frame, LineColor( kRed), Range("full") );
		frame->Draw();
		return frame;
	}

	/**
	 * Roo Variables
	 */
	bool useROI;
	RooDataHist * rdh;
	RooRealVar *x, * mean, *sigma;
	RooGaussian* gaussModel;
};




#endif