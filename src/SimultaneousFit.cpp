

#include "SimultaneousFit.h"
#include "constants.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include "RooPlot.h"

#include "Reporter.h"


SimultaneousFit::SimultaneousFit( 	TH1D* tAll, TH1D* tS0, TH1D* tS1, TH1D* tS2,
									TH1D* dAll, TH1D* dS0, TH1D* dS1, TH1D* dS2,
									double p, XmlConfig * config, string nodePath 
		){


	lg = LoggerConfig::makeLogger( config, "Logger" );

	tofAll = tAll;
	tofS0 = tS0;
	tofS1 = tS1;
	tofS2 = tS2;

	dedxAll = dAll;
	dedxS0 = dS0;
	dedxS1 = dS1;
	dedxS2 = dS2;

	this->avgP = p;

	this->config = config;
	path = nodePath;

	tGen = new tofGenerator(  );

}

SimultaneousFit::~SimultaneousFit(){

}


void SimultaneousFit::fitTofAll( ){

	lg->info(__FUNCTION__) << "Fitting to Tof All tracks Only " << endl;

	Reporter rp( "tofAll.pdf" );

	// get some limits
	double c = tGen->mean( avgP, constants::kaonMass );
	double x1 = (tGen->mean( avgP, constants::piMass ) - c);
	double x2 = (tGen->mean( avgP, constants::protonMass ) - c );

	lg->info() << " p = " << avgP << endl;
	lg->info() << " center " << c << endl;
	lg->info() << " x1 " << x1 << endl;
	lg->info() << " x2 " << x2 << endl;

	// build a 3 gaussian model
	RooRealVar *tof = new RooRealVar( "tofx", "tofx", x1, x2 );
	RooDataHist * rdh = new RooDataHist( "data", "data", RooArgSet( * tof ), tofAll  );
	

	RooPlot * frame = tof->frame( );

	rdh->plotOn( frame );

	frame->Draw();

	rp.savePage();

	tofAll->Draw();
	rp.savePage();

}
















