#ifndef PID_SQUARE_H
#define PID_SQUARE_H
#endif
/*
	RooFit Includes
*/
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooProdPdf.h"
#include "RooAbsReal.h"
#include "RooPolynomial.h"
/*
	Root Includes
*/
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include "vector"
#include "sstring"
#include "string";

class pidSquare
{

protected:
	TH2* dataHist;
	std::vector<TH2*> truth;
	// 0 = pion
	// 1 = kaon
	// 2 = proton
	int centerSpecies = -1;
	double aspect = 1.0;

public:
	pidSquare( 	TH2* data, int cSpecies, double aspect,
				TH2* p0, TH2* p1, TH2* p2 ){

		this->dataHist = data;
		this->centerSpecies = cSpecies;

		truth.push_back( p0 );
		truth.push_back( p1 );
		truth.push_back( p2 );

	}
	~pidSquare();

	double yield( 	TH2* h, 	
					double xLow, double xHigh,
					double yLow, double yHigh ){

		int xBinLow = h->GetXaxis()->FindBin( xLow );
		int xBinHigh = h->GetXaxis()->FindBin( xHigh );

		int yBinLow = h->GetYaxis()->FindBin( yLow );
		int yBinHigh = h->GetYaxis()->FindBin( yHigh );

		double yield = h->Integral( xBinLow, xBinHigh, yBinLow, yBinHigh );

		return yield;

	}
	TH1D* efficiency( string name, double low = 0.5, double high = 5.5, double step = 0.25){

		double cutLow = low;
		double cutHigh = high;
		double cutStep = step;

		int n = (cutHigh - cutLow ) / cutStep;
		
		TH1D* effPlot = new TH1D( name.c_str(), "efficiency; Cut( N_{#sigma} dedx X N_{#sigma} 1/#beta )", n, cutLow, cutHigh );


		TH2* hTruth = truth[ centerSpecies ];
		double total = hTruth->Integral();
		int bin = 0;

		for ( double cut = cutLow; cut < cutHigh; cut += cutStep ){

			double inCut = yield( hTruth, -cut, cut, -(cut*aspect), cut*aspect );

			effPlot->SetBinContent( bin, (inCut / total ) );
			bin++;

		}

		return effPlot;
	}
	/*TH1D* purity( double low = 0.5, double high = 5.5, double step = 0.25){

		double cutLow = low;
		double cutHigh = high;
		double cutStep = step;
		int n = (cutHigh - cutLow ) / cutStep;
		TH1D* purePlot = new TH1D( "purity", "purity", n, cutLow, cutHigh );

		TH1*proj = truth[ centerSpecies ]->ProjectionX();
		if ( !cutOnDedx )
			proj = truth[ centerSpecies ]->ProjectionY();

		
		int bin = 0;
		for ( double cut = cutLow; cut < cutHigh; cut += cutStep ){

			int bLow = proj->GetXaxis()->FindBin( -1 * cut );
			int bHigh = proj->GetXaxis()->FindBin( 1 * cut );

			double total = 0;
			for ( int i = 0; i < truth.size(); i++ ){
				TH1* pTruth = truth[ i ]->ProjectionX();
				
				if ( !cutOnDedx )
					pTruth = truth[ i ]->ProjectionY();

				total += pTruth->Integral( bLow, bHigh );
			}

			double inCut = proj->Integral( bLow, bHigh );

			purePlot->SetBinContent( bin, (inCut / total ) );
			bin++;

		}

		return purePlot;
	}
*/
	
};