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
	double yield( TH2D* dHist, double rDedx, double rTof ){

		


	}


protected:

	
};


#endif