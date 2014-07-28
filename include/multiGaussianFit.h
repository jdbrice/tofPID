#ifndef MULTI_GAUSSIAN_FIT_H
#define MULTI_GAUSSIAN_FIT_H

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


class multiGaussianFit
{
public:
	multiGaussianFit( TH2 * h, uint nSpecies );
	~multiGaussianFit();

	void setX( string name, double min, double max ){
		xObs = new RooRealVar( name.c_str(), name.c_str(),
								min, max );
		xMin = min;
		xMax = max;
	}
	void setY( string name, double min, double max ){
		yObs = new RooRealVar( name.c_str(), name.c_str(),
								min, max );
		yMin = min;
		yMax = max;
	}

	void fit();

	RooRealVar * getX(){
		return xObs;
	}
	RooRealVar * getY(){
		return yObs;
	}
	RooAddPdf * getFit(){
		return model;
	}
	RooPlot * viewFitX(){
		RooPlot * frame = xObs->frame(  );
		//rdh->plotOn( frame ); //DataError( RooAbsData::SumW2 )
		TH2* h = xObs->createHistogram( "rooData", *yObs, 0, 0, 0, 0);
		rdh->fillHistogram( h, RooArgList( *xObs, *yObs) );
		
		TH2* h2 = xObs->createHistogram( "rooModel", *yObs, 0, 0, 0, 0);
		model->fillHistogram( h2, RooArgList( *xObs, *yObs) );
		h->Draw("colz");
		h2->Draw("cont2 same");
		//frame->Draw();
		return frame;
	}

protected:

	RooAddPdf * gxy;
	void setupModel();

	TH2 * input;
	uint nSpecies;
	double xMin, xMax;
	double yMin, yMax;

	RooDataHist * rdh;

	RooRealVar * xObs;
	RooRealVar * yObs;

	vector<RooRealVar * > meanX;
	vector<RooRealVar * > meanY;
	
	vector<RooRealVar * > sigX;
	vector<RooRealVar * > sigY;

	vector<RooGaussian * > gaussX;
	vector<RooGaussian * > gaussY;

	// the scale or amplitude of the gaussian
	vector<RooRealVar * > amp;

	vector<RooProdPdf *> gauss2D;
	RooAddPdf * model;

};

#endif