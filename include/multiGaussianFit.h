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
	void viewFitX( string title = "" ){
		RooPlot * frame = xObs->frame(  );
		
		TH2* h = xObs->createHistogram( "rooData", *yObs, 0, 0, 0, 0);
		rdh->fillHistogram( h, RooArgList( *xObs, *yObs) );
		
		TH2* h2 = xObs->createHistogram( "rooModel", *yObs, 0, 0, 0, 0);
		model->fillHistogram( h2, RooArgList( *xObs, *yObs) );
		
		if ( "" != title )
			h->SetTitle( title.c_str() );
		h->Draw("colz");
		h2->Draw("cont2 same");
		
	}

	void setInitialMean( double x, double y ){
		initialMeanX.push_back( x );
		initialMeanY.push_back( y );
	}
	void limitMeanX( double x1, double x2 ){
		meanMinX.push_back( x1 );
		meanMaxX.push_back( x2 );
	}
	void limitMeanY( double y1, double y2 ){
		meanMinY.push_back( y1 );
		meanMaxY.push_back( y2 );
	}

	double getMeanX( uint iS ){
		if ( iS < nSpecies )
			return meanX[ iS ]->getVal();
	}
	double getMeanY( uint iS ){
		if ( iS < nSpecies )
			return meanY[ iS ]->getVal();
	}
	double getMeanXError( uint iS ){
		if ( iS < nSpecies )
			return meanX[ iS ]->getError();
	}
	double getMeanYError( uint iS ){
		if ( iS < nSpecies )
			return meanY[ iS ]->getError();
	}

	double getSigmaX( uint iS ){
		if ( iS < nSpecies )
			return sigX[ iS ]->getVal();
	}
	double getSigmaY( uint iS ){
		if ( iS < nSpecies )
			return sigY[ iS ]->getVal();
	}
	double getSigmaXError( uint iS ){
		if ( iS < nSpecies )
			return sigX[ iS ]->getError();
	}
	double getSigmaYError( uint iS ){
		if ( iS < nSpecies )
			return sigY[ iS ]->getError();
	}

protected:

	RooAddPdf * gxy;
	void setupModel();

	TH2 * input;
	uint nSpecies;
	double xMin, xMax;
	double yMin, yMax;

	vector<double> initialMeanX;
	vector<double> initialMeanY;
	vector<double> meanMinX;
	vector<double> meanMinY;
	vector<double> meanMaxX;
	vector<double> meanMaxY;

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