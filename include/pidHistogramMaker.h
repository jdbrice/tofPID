#ifndef pidHistogramMaker_H
#define pidHistogramMaker_H

#include "allroot.h"

#include "histoBook.h"
#include "constants.h"
#include "TOFrPicoDst.h"
#include <vector>

// clock_t, clock, CLOCKS_PER_SEC 
#include <time.h>       

// for testing if stdout is interactive or pipe / file
#include "xmlConfig.h"
#include "utils.h"
#include "reporter.h"


class pidHistogramMaker{

private:

	// the canvas used to draw report hidtos
	reporter* report;

	
	// the main chain object
	TChain * _chain;

	// the histobook that stores all of our pidHistogramMaker histograms
	histoBook *book;

	// the pico dst for simpler chain usage
	TOFrPicoDst * pico;

	// config file
	xmlConfig* config;

	clock_t startTime;

	map< string, reporter * > pReport;


	// removing overflows
	double pMax, pMin, nSigMax, nSigMin;
	double dBetaMax, dBetaMin;

	//QA memebers
	double vOffsetX, vOffsetY;

public:


	// Constructor
	pidHistogramMaker( TChain * chain, xmlConfig *config );

	// destructor
	~pidHistogramMaker();
	
	void makeQA();
	void make();

	TGraph* inverseBeta( double m, double p1, double p2, double step = .05 );
	
	bool keepEventQA();
	bool keepTrackQA( uint iHit );
	

protected:

	void sHisto( string pType );
	string sName( string pType, int charge );
	void speciesReport( string pType, int charge, int etaBin = -1 );

	double nSigDedx( string pType, int iHit ) { 
		if ( "P" == pType )
			return pico->nSigP[ iHit ];
		if ( "K" == pType )
			return pico->nSigK[ iHit ];
		if ( "Pi" == pType )
			return pico->nSigPi[ iHit ];
		return -999.0;
	}
	double nSigInvBeta( string pType, int iHit  );
	double dBeta( string pType, int iHit  );

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

	void dklFit( string name, TH2D* h );

	/*
	*	Utility functions that should be moved soon
	*/ 
	void startTimer( ) { startTime = clock(); }
	double elapsed( ) { return ( (clock() - startTime) / (double)CLOCKS_PER_SEC ); }
};



#endif