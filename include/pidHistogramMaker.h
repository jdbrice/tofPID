#ifndef pidHistogramMaker_H
#define pidHistogramMaker_H

#include "allroot.h"

#include "HistoBook.h"
#include "constants.h"
#include "TOFrPicoDst.h"
#include "tofGenerator.h"
#include "dedxGenerator.h"
#include <vector>

// clock_t, clock, CLOCKS_PER_SEC 
#include <time.h>       

// for testing if stdout is interactive or pipe / file
#include "XmlConfig.h"
#include "Utils.h"

using namespace jdb;

#include "reporter.h"


class pidHistogramMaker{

private:

	// the canvas used to draw report hidtos
	Reporter* report;

	
	// the main chain object
	TChain * _chain;

	// the histobook that stores all of our pidHistogramMaker histograms
	HistoBook *book;

	TFile * distroData;

	// the pico dst for simpler chain usage
	TOFrPicoDst * pico;

	// config file
	XmlConfig* config;

	clock_t startTime;

	map< string, Reporter * > pReport;


	// binning information
	double pMax, pMin;
	vector<double> pBins;
	double dedxMax, dedxMin;
	vector<double> dedxBins;
	double tofMax, tofMin;
	vector<double> tofBins;

	//QA memebers
	double vOffsetX, vOffsetY;

	// have the main histograms been prepared
	bool histosReady;

	// tof Metric
	string tofMetric;
	static const string inverseBeta;
	static const string deltaBeta;
	double inverseBetaSigma;
	tofGenerator * tofGen;

	double dedxSigma;
	dedxGenerator * dedxGen;

	string centeringMethod;
	double dedxShift;
	double tofShift;
	static const string traditionalCentering;
	static const string nonlinearCentering;

	static const vector<string> species;

	// for padding plot ranges
	double tofPadding, dedxPadding;
	double tofScalePadding, dedxScalePadding;

	vector<double> averageP;

public:


	// Constructor
	pidHistogramMaker( TChain * chain, XmlConfig *config );

	// destructor
	~pidHistogramMaker();
	
	void momentumDistributions();
	void makeQA();
	void makePidHistograms();

	TGraph* inverseBetaGraph( double m, double p1, double p2, double step = .05 );
	
	bool keepEventQA();
	bool keepTrackQA( uint iHit );
	

protected:

	
	void prepareHistograms( string pType );
	string speciesName( string pType, int charge = 0);
	void speciesReport( string pType, int charge, int etaBin = -1 );
	void distributionReport( string pType );

	double nSigDedx( string pType, int iHit ); 
	double nSigmaDedx( string pType, int iHit, double avgP );

	double lh( double t, double mu, double sigma );

	double nSigInvBeta( string pType, int iHit  ){

		double betaMeasured = pico->beta[ iHit ];
		double p = pico->p[ iHit ];
		double betaExpected = eBeta( eMass( pType ), p );
		
		double deltaInvBeta = ( 1.0 / betaMeasured ) - ( 1.0 / betaExpected );

		return (deltaInvBeta / inverseBetaSigma);
	}
	double nSigmaInverseBeta( string pType, int iHit, double avgP );

	double dBeta( string pType, int iHit  ){

		//double tof = pico->tof[ iHit ];
		//double length = pico->length[ iHit ];
		double p = pico->p[ iHit ];
		double beta = pico->beta[ iHit ];
		//double m2 = p*p * ( constants::c*constants::c * tof*tof / ( length*length ) - 1  );


		double deltaB = 1 - (beta) * TMath::Sqrt( (constants::kaonMass*constants::kaonMass) / (p*p) + 1 );

		if ( "Pi" == pType )
			deltaB = 1 - (beta) * TMath::Sqrt( (constants::piMass*constants::piMass) / (p*p) + 1 );		
		if ( "P" == pType )
			deltaB = 1 - (beta) * TMath::Sqrt( (constants::protonMass*constants::protonMass) / (p*p) + 1 );		
		
		return deltaB;
	}

	double eBeta( double m, double p ){
		return TMath::Sqrt( p*p / ( p*p + m*m ) );
	}
	double eMass( string pType ){
		if ( "P" == pType )
			return constants::protonMass;
		if ( "K" == pType )
			return constants::kaonMass;
		if ( "Pi" == pType )
			return constants::piMass;
		return -10.0;	
	}

	void autoViewport( 	string pType, double p, double * tofLow, double* tofHigh, double * dedxLow, double * dedxHigh, 
						double tofPadding = 1, double dedxPadding = 1, double tofScaledPadding = 0, double dedxScaledPadding = 0 );


	/*
	*	Utility functions that should be moved soon
	*/ 
	void startTimer( ) { startTime = clock(); }
	double elapsed( ) { return ( (clock() - startTime) / (double)CLOCKS_PER_SEC ); }
};



#endif