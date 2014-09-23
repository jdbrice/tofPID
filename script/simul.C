


#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include "RooPlot.h"


void simul(){

	TFile* in = new TFile( "../bin/Kdata.root", "READ" );

	TH1D* h = in->Get("tof_8");

	RooRealVar * x = new RooRealVar( "x", "x", -10, 10 );

	RooDataHist * rdh = new RooDataHist( "data", "data", RooArgSet( *x ), h);

	RooRealVar m( "mean", "mean", 0, -1, 1 );
	RooRealVar s( "sigma", "sigma", .5, 0, 1 );
	RooGaussian g( "gauss", "gauss", )

	RooPlot * frame = x->frame();
	rdh->plotOn( frame );
	frame->Draw();
	

}