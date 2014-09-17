
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include "RooPlot.h"
#include "RooPolynomial.h"

#include "sstream"
using namespace std;


void slice( int pBin = 0 );

void rooLinearGaussian( TH1D* h );

TH3D* h3;		// current 3d histo
TH2D* pSlice;	// 2D slice for a given P bin
TH1D* bSlice;	// 1D slice of the delta beta inv for the current p bin
TH1D* dSlice;	// 1D slice of the dedx "" "" 

void tofFit() {


	TH1D * hMean = new TH1D( "mean", "mean", 28, 0, 28 );
	TH1D * hSigma = new TH1D( "sigma", "sigma", 28, 0, 28 );

	TFile* data = new TFile( "invBetaData.root", "READ" );

	h3 = (TH3D*) data->Get( "nSig_K" );


	TCanvas * c = new TCanvas( "c", "K Fit", 800, 600 );
	string outName = "rpTofFit.pdf";
	c->Print( (outName+"[").c_str() );
	for ( int i = 1 ; i < 26; i ++ ){
		slice( i );
		

		rooLinearGaussian( i, bSlice );
		gPad->SetLogy(1);

		c->Print( (outName).c_str() );


		hMean->SetBinContent( i-1, rrvMean->getVal() );
		hMean->SetBinError( i-1, rrvMean->getError() );
		hSigma->SetBinContent( i-1, rrvSigma->getVal() );
		hSigma->SetBinError( i-1, rrvSigma->getError() );


	}

	gPad->SetLogy(0);
	hMean->Draw( "pe" );

	c->Print( (outName).c_str() );

	hSigma->Draw( "pe" );
	c->Print( (outName).c_str() );

	c->Print( (outName+"]").c_str() );
	

}

void slice( int pBin ){

	h3->GetZaxis()->SetRange( pBin, pBin );
	pSlice = (TH2D*)h3->Project3D("xy");
	bSlice = (TH1D*)h3->Project3D( "y" );
	dSlice = (TH1D*)h3->Project3D( "x" );
}


RooRealVar * rrvMean;
RooRealVar * rrvSigma;

void rooLinearGaussian( int pBin, TH1D* h ){
	using namespace RooFit;
	

	double b1 = bSlice->FindFirstBinAbove( 50, 1 );
	double b2 = bSlice->FindLastBinAbove( 50, 1 );
	double x1 = bSlice->GetBinLowEdge( b1 );
	double x2 = bSlice->GetBinLowEdge( b2 ) + bSlice->GetBinWidth( b2 );

	if ( x1 > -10 )
		x1 = -10;
	if ( x2 < 10 )
		x2 = 10;

	RooRealVar x("x", "x", x1, x2 );
	RooDataHist rdh( "data", "data", RooArgSet( x ), h );

	RooAbsPdf *model;

	if ( pBin <= 10 ){
		model = makeGauss( "gauss1_", &x );
		
		x.setRange( "cr", -8 , 8 );
	} else if ( pBin <= 20 ) {
		
		model = makeGauss( "gauss1_", &x );
		double mod = ( pBin - 10 ) * .5;
		x.setRange( "cr", -6 + mod , 6 - mod );
	} else if ( pBin <= 22 ) {	
		model = makeGauss( "gauss1_", &x );
		x.setRange( "cr", -2 , 2 );
	} else if ( pBin <= 30 ) {	
		model = makeGauss( "gauss1_", &x );
		x.setRange( "cr", -3 , 3 );
	}


	stringstream sstr;
	sstr << " pBin : " << pBin;
	RooPlot * frame = x.frame( Title( sstr.str().c_str() ) );

	// Plot p.d.f. 
	
	rdh.plotOn( frame );
	model->fitTo( rdh, Range( "cr" ), Extended() );


	model->plotOn( frame, LineColor( kRed), Range("full") );
	
	
	frame->Draw();

	delete model;
	

}



RooAddPdf * makeSpecies(	int nSpecies, RooRealVar *x ) {
	using namespace RooFit;

	if ( nSpecies == 1 ){
		cout << " one species " << endl;
		//RooRealVar x("x", "x", x1, x2);
		RooRealVar *g1 = new RooRealVar( "mean", "gaussian mean", 0, -0.005, 0.005 );
		RooRealVar *g2 = new RooRealVar( "sigma", "gaussian sigma", .01, 0, .05 );

		RooGaussian gauss( "gauss", "gauss", *x, *g1, *g2 );

		RooRealVar *h1 = new RooRealVar( "hmean", "haussian mean", 0, -0.05, 0.05 );
		RooRealVar *h2 = new RooRealVar( "hsigma", "haussian sigma", .01, 0, .05 );

		RooGaussian hauss( "hauss", "hauss", *x, *h1, *h2 );

		RooRealVar p1( "p1", "p1", 1, -1000, 1000 );
		RooRealVar p2( "p2", "p2", 1, -1000, 1000 );
		RooRealVar p3( "p3", "p3", 1, -1000, 1000 );
		RooRealVar p4( "p4", "p4", 1, -1000, 1000 );

		RooPolynomial poly( "poly", "poly", *x, RooArgList( p1, p2, p3, p4 ) );

		RooRealVar w1( "w1", "w1", 1, 0, 10000000 );
		RooRealVar w2( "w2", "w2", 1, 0, 10000000 );
		RooRealVar w3( "w2", "w2", 1, 0, 10000000 );

		RooAddPdf * model = new RooAddPdf( "model", "model", RooArgList( poly, gauss, hauss), RooArgList( w1, w2, w3 ) );
		return model;
	}	

	return  (new RooAddPdf( "model", "model" ));
}

RooGaussian * makeGauss( string name, RooRealVar *x ) {
	using namespace RooFit;
	
	rrvMean = new RooRealVar( (name+"mean").c_str(), "gaussian mean", 0, -5, 5 );
	rrvSigma = new RooRealVar( (name+"sigma").c_str(), "gaussian sigma", 1, 0.00, 5 );

	RooGaussian * gauss = new RooGaussian( name.c_str(), "gaussian", *x, *rrvMean, *rrvSigma );

	return gauss;
}

RooAbsPdf * makeDoubleGauss( string n1, string n2, RooRealVar *x ) {
	using namespace RooFit;
	
	RooRealVar *c1 = new RooRealVar( (n1+"mean").c_str(), "gaussian mean", 0, -.03, .03 );
	RooRealVar *c2 = new RooRealVar( (n1+"sigma").c_str(), "gaussian sigma", .01, 0.005, .1 );

	RooGaussian * gauss1 = new RooGaussian( n1.c_str(), "gaussian", *x, *c1, *c2 );

	RooRealVar *c3 = new RooRealVar( (n2+"mean").c_str(), "gaussian mean", -.10, -.15, -.35 );
	RooRealVar *c4 = new RooRealVar( (n2+"sigma").c_str(), "gaussian sigma", .01, 0.00, .15 );

	RooGaussian * gauss2 = new RooGaussian( n1.c_str(), "gaussian", *x, *c3, *c4 );

	RooRealVar* w1 = new RooRealVar( "w1", "w1", 1, 0, 10000000 );
	RooRealVar* w2 = new RooRealVar( "w2", "w2", 1, 0, 10000000 );

	RooAddPdf *model = new RooAddPdf( (n1+n2).c_str(), "double gaussian", RooArgList( *gauss1, *gauss2 ), RooArgList( *w1, *w2) );

	return gauss1;
}



















