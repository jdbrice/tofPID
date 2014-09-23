
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


void t3g() {


	TFile* data = new TFile( "../bin/histogram/dqa.root", "READ" );

	TH1D* h1 = (TH1D*) data->Get( "tof/tof_s0_15" );
	

	
	RooRealVar x( "x", "x", -100, 100 );
	RooDataHist rdh( "rdh", "rdh", x, h1 );

	RooPlot * frame = x.frame( RooFit::Title( "hello" ) );
	rdh.plotOn( frame );

	frame->Draw();




}
















