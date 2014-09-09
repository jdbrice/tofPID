#ifndef PID_ELLIPSE_H
#define PID_ELLIPSE_H 

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
#include "TMath.h"

#include "vector"

class pidEllipse
{

protected:
	TH2* dataHist;
	std::vector<TH2*> truth;
	// 0 = pion
	// 1 = kaon
	// 2 = proton
	int centerSpecies = -1;

	// rDedx / rTof
	double aspect = 1.0;


public:
	pidEllipse( TH2* data, int cSpecies, double aspect,
				TH2* p0, TH2* p1, TH2* p2 ){


		this->dataHist = data;
		this->aspect = aspect;
		this->centerSpecies = cSpecies;

		truth.push_back( p0 );
		truth.push_back( p1 );
		truth.push_back( p2 );

	}
	~pidEllipse(){}

	/**
	 * Calculates the yield inside the cut region
	 * @param  dHist The data histogram to use, 
	 *               this can be the sum or a truth histo
	 * @param  rDedx Radius of the dedx cut
	 * @param  rTof  Radius of the 1/beta cut
	 * @return       The absolute yield inside the cut ellipse
	 */
	double yield( TH2* dHist, double rDedx, double rTof ){

		int dedxBinHigh = dHist->GetXaxis()->FindBin( rDedx );
		int dedxBinZero = dHist->GetXaxis()->FindBin( 0.0 );

		int tofBinHigh = dHist->GetYaxis()->FindBin( rTof );
		int tofBinZero = dHist->GetYaxis()->FindBin( 0.0 );

		int dedxBinR = dedxBinHigh - dedxBinZero;
		int tofBinR = tofBinHigh - tofBinZero;
		
		double yTotal = 0;
		// loop through the bins along one axis
		for ( int iX = dedxBinR; iX >= 0; iX-- ){

			double alpha = (iX/(double)dedxBinR)*(iX/(double)dedxBinR);
			double beta = (double)tofBinR*(double)tofBinR;
			double y = TMath::Sqrt(  ( 1.0 - alpha ) *  beta );
			int yP = tofBinZero + y;
			int yN = tofBinZero - y;

			// positive iX
			yTotal += dHist->Integral( dedxBinZero+iX, dedxBinZero+iX, yN, yP );

			if ( 0 != iX )
				yTotal += dHist->Integral( dedxBinZero-iX, dedxBinZero-iX, yN, yP );

		}

		//cout << yTotal << " for (x/" << rDedx << ")^2 + (y/" << rTof << ")^2 < 1" << endl;
		//cout << dHist->Integral( dedxBinZero-dedxBinR, dedxBinZero+dedxBinR, tofBinZero-tofBinR, tofBinZero+tofBinR) << " Square " << endl;
		
		return yTotal;

	}

	double efficiency( double cut = 1.0 ){

		TH2* hTruth = truth[ centerSpecies ];
		double total = hTruth->Integral();
		double inCut = yield( hTruth, cut, cut*aspect );
		return (inCut / total);

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

			double inCut = yield( hTruth, cut, cut*aspect );

			effPlot->SetBinContent( bin, (inCut / total ) );
			bin++;

		}

		return effPlot;
	}

	double purity( double cut = 1.0 ){

		TH2* hTruth = truth[ centerSpecies ];
		double total = yield( dataHist, cut, cut*aspect );
		double inCut = yield( hTruth, cut, cut*aspect );
		if ( 0 < total )
			return inCut / total;
		return 0.0;
	}
	TH1D* purity( string name, double low = 0.5, double high = 5.5, double step = 0.25){

		double cutLow = low;
		double cutHigh = high;
		double cutStep = step;
		int n = (cutHigh - cutLow ) / cutStep;

		TH1D* purePlot = new TH1D( name.c_str(), "purity", n, cutLow, cutHigh );

		TH2* hTruth = truth[ centerSpecies ];
		
		int bin = 0;
		for ( double cut = cutLow; cut < cutHigh; cut += cutStep ){

			double total = yield( dataHist, cut, cut*aspect );
			double inCut = yield( hTruth, cut, cut*aspect );

			//cout << name << " : cut = " << cut << " total : " << total << " inCut : " << inCut << endl;
			if ( total > 0 )
				purePlot->SetBinContent( bin, (inCut / total ) );
			else
				purePlot->SetBinContent( bin, 0.0 );
			bin++;
		}

		return purePlot;
	}


protected:

	
};


#endif