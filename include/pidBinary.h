
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

class pidBinary
{

protected:
	TH2* dataHist;
	std::vector<TH2*> truth;
	// 0 = pion
	// 1 = kaon
	// 2 = proton
	int centerSpecies = -1;
	bool cutOnDedx = true;

public:
	pidBinary( TH2* data, int cSpecies, TH2* p0, TH2* p1, TH2* p2 ){

		this->dataHist = data;
		this->centerSpecies = cSpecies;

		truth.push_back( p0 );
		truth.push_back( p1 );
		truth.push_back( p2 );

	}
	~pidBinary();

	void cutDedx( bool yes ) {
			cutOnDedx = yes;
	}
	bool cutDedx() { 
		return cutOnDedx; 
	}

	TH1D* cutData( double cutLow, double cutHigh ){

		TH1*proj = (TH1D*)dataHist->ProjectionX()->Clone("cutData" );
		if ( !cutOnDedx )
			proj = (TH1D*)dataHist->ProjectionY()->Clone("cutData" );

		int bLow = proj->GetXaxis()->FindBin( cutLow );
		int bHigh = proj->GetXaxis()->FindBin( cutHigh );

		proj->GetXaxis()->SetRange( bLow, bHigh );
		proj->SetFillColor( kRed );

		return (TH1D*)proj;

	}

	double yield( TH1 * h, double cut ){

		int bLow = h->GetXaxis()->FindBin( -cut );
		int bHigh = h->GetXaxis()->FindBin( cut );
		return h->Integral( bLow, bHigh );
	}

	double efficiency( double cut = 1.0 ){
		TH1*proj = truth[ centerSpecies ]->ProjectionX();
		if ( !cutOnDedx ){
			proj = truth[ centerSpecies ]->ProjectionY();
		}
		double total = proj->Integral();
		double inCut = yield( proj, cut );
		if ( total > 0 )
			return inCut / total;
		return 0;
	}
	TH1D* efficiency( string name, double low = 0.5, double high = 5.5, double step = 0.25){

		double cutLow = low;
		double cutHigh = high;
		double cutStep = step;
		int n = (cutHigh - cutLow ) / cutStep;
		TH1D* effPlot = new TH1D( name.c_str(), "efficiency; Cut( N_{#sigma} dedx )", n, cutLow, cutHigh );


		TH1*proj = truth[ centerSpecies ]->ProjectionX();
		if ( !cutOnDedx ){
			proj = truth[ centerSpecies ]->ProjectionY();
			effPlot->GetXaxis()->SetTitle( "Cut ( N_{#sigma} 1/#beta )" );
		}

		double total = proj->Integral();
		int bin = 0;
		for ( double cut = cutLow; cut < cutHigh; cut += cutStep ){

			
			double inCut = yield( proj, cut );

			effPlot->SetBinContent( bin, (inCut / total ) );
			bin++;

		}

		return effPlot;
	}

	double purity( double cut = 1.0 ) {
		TH1*proj = truth[ centerSpecies ]->ProjectionX();
		if ( !cutOnDedx ){
			proj = truth[ centerSpecies ]->ProjectionY();
		}

		TH1*projSum = dataHist->ProjectionX();
		if ( !cutOnDedx ){
			projSum = dataHist->ProjectionY();
		}
		double total = yield( projSum, cut );
		double inCut = yield( proj, cut );
		if ( total > 0 )
			return inCut / total;
		return 0;

	}
	TH1D* purity( string name, double low = 0.5, double high = 5.5, double step = 0.25){

		double cutLow = low;
		double cutHigh = high;
		double cutStep = step;
		int n = (cutHigh - cutLow ) / cutStep;
		TH1D* purePlot = new TH1D( name.c_str(), "purity", n, cutLow, cutHigh );

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

	
};