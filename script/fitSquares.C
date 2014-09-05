
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


Double_t background( Double_t * x, Double_t * par);
Double_t gaussian( Double_t * x, Double_t * par);
Double_t combFit( Double_t * x, Double_t * par);
void slice( int pBin = 0 );

void rooLinearGaussian( TH1D* h );

TH3D* h3;		// current 3d histo
TH2D* pSlice;	// 2D slice for a given P bin
TH1D* bSlice;	// 1D slice of the delta beta inv for the current p bin
TH1D* dSlice;	// 1D slice of the dedx "" "" 

void fitSquares() {


	TFile* data = new TFile( "dqa.root", "READ" );

	h3 = (TH3D*) data->Get( "nSig_K_All_eta0" );

	slice( 1 );

	//bSlice->Draw();
	//gPad->SetLogy(1);
	
/*
	TF1* bFit = new TF1( "bFit", combFit, -1, 1, 5);
	bFit->SetParameters( 1, 1, 1, 1, 1);
	bFit->SetNpx( 500 );
	bSlice->Fit( bFit );
*/

	//bSlice->Draw( "pe" );
	//bSlice->Sumw2();
	//dSlice->Sumw2();
	rooLinearGaussian( bSlice );
	//dSlice->Draw();

}

void slice( int pBin ){

	h3->GetZaxis()->SetRange( 1, 1 );
	pSlice = (TH2D*)h3->Project3D("xy");
	bSlice = (TH1D*)h3->Project3D( "y" );
	dSlice = (TH1D*)h3->Project3D( "x" );
}

Double_t background( Double_t * x, Double_t * par){
	return par[0] + par[ 1 ] * x[ 0 ];
	
}
Double_t gaussian( Double_t * x, Double_t * par){
	double alpha = (x[0] - par[ 1 ])*(x[0] - par[ 1 ]);
	double beta = 2 * par[ 2 ] * par[ 2 ];
	return par[ 0 ] * TMath::Exp( - (alpha / beta) );
}
Double_t combFit( Double_t * x, Double_t * par){
	return background(x, par) +  gaussian( x, &par[2]);
}


using namespace RooFit;

void rooLinearGaussian( TH1D* h ){

	

	RooRealVar x("x", "x", -.15, .15);
	RooDataHist rdh( "data", "data", RooArgSet( x ), h );

	RooRealVar g1( "mean", "gaussian mean", 0, -1, 1 );
	RooRealVar g2( "sigma", "gaussian sigma", .01, 0, 10 );

	RooGaussian gauss( "gauss", "gauss", x, g1, g2 );

	
	
	
	RooRealVar tau("tau","tau",.01548) ;
	// Build a gaussian resolution model
	RooRealVar bias1("bias1","bias1",0) ;
	RooRealVar sigma1("sigma1","sigma1",.01) ;
	RooGaussModel gm1("gm1", "gauss model 1", x, bias1, sigma1);

	// Construct decay(t) (x) gauss1(t)
	RooDecay decay_gm1("decay_gm1", "decay", x, tau, gm1, RooDecay::DoubleSided );

	RooPlot * frame = x.frame();

	// Plot p.d.f. 
	decay_gm1.plotOn(frame) ;


	

	
	/*
	data->plotOn( frame );
	
	gauss.fitTo( *data );

	gauss.plotOn( frame );
	*/
	frame->Draw();
	//gPad->SetLogy();

}



















