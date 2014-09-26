#ifndef pidHistogramMaker_H
#define pidHistogramMaker_H

#include "allroot.h"

#include "HistoBook.h"
#include "constants.h"
#include "TOFrPicoDst.h"
#include "tofGenerator.h"
#include "Bichsel.h"
#include <vector>

// clock_t, clock, CLOCKS_PER_SEC 
#include <time.h>       

// for testing if stdout is interactive or pipe / file
#include "XmlConfig.h"
#include "Utils.h"
#include "Logger.h"
#include "LoggerConfig.h"

using namespace jdb;

#include "reporter.h"


class pidHistogramMaker{

private:

	Logger * lg;

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
	// pt Bins
	double ptMax, ptMin;
	vector<double> ptBins;
	// eta bins
	double etaMax, etaMin;
	vector<double> etaBins;
	// dedx bins
	double dedxMax, dedxMin;
	vector<double> dedxBins;
	// tof bins
	double tofMax, tofMin;
	vector<double> tofBins;
	// charge bins
	double chargeMax, chargeMin;
	vector<double> chargeBins;

	//QA memebers
	double vOffsetX, vOffsetY;

	// have the main histograms been prepared
	bool histosReady;

	// tof Metric
	string tofMetric;
	static const string inverseBeta;
	static const string deltaBeta;
	double tofSigma, tofPlotSigma;
	tofGenerator * tofGen;

	double dedxSigma, dedxPlotSigma;
	Bichsel * dedxGen;

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
	void makeDedxTofHistograms();

	TGraph* inverseBetaGraph( double m, double p1, double p2, double step = .05 );
	
	bool keepEventQA();
	bool keepTrackQA( uint iHit );
	
	void prepareHistograms( string pType );

protected:


	void speciesReport( string pType, int charge, int etaBin = -1 );
	void distributionReport( string pType );

	/**
	 * Calculates the nSigma dEdx value with the centroid of plc pType
	 * set to zero
	 * @param  pType The plc species for centering
	 * @param  iHit  The hit index in the pico-dst
	 * @return       Returns the # of nSigmas from the center species.
	 * Uses the "plotSigma" to allow this to be either raw values
	 * or nSigmas calculated using a fixed sigma value
	 */
	double nSigmaDedx( string pType, int iHit ); 
	double nSigmaDedx( string pType, int iHit, double avgP );

	/**
	 * Calculates the nSigma 1/beta value with centering around plc species pType
	 * @param  pType Centering species
	 * @param  iHit  The hit index in the pico-dst
	 * @return       Returns the # of nSigma from the center species. 
	 * See nSigmaDedx for notes.
	 */
	double nSigmaInverseBeta( string pType, int iHit  );
	double nSigmaInverseBeta( string pType, int iHit, double avgP );

	/**
	 * A gaussian to calculate the likelihood of being a plc species
	 * @param  t     the input measurement
	 * @param  mu    the mean of the gauss
	 * @param  sigma the sigma to use
	 * @return       Returns the unnormalized likelihood value
	 */
	double lh( double t, double mu, double sigma );

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

	string chargeString( int charge = 0 ) {
		if ( -1 >= charge )	// negative
			return "n";
		else if ( 1 <= charge ) //positive
			return "p";
		return "a";	// all
	}
	string speciesName( string centerSpecies, int charge = 0 ){
		return "dedx_tof_" + chargeString(charge) + "_" + centerSpecies;
	}
	string speciesName( string centerSpecies, int charge, int ptBin, int etaBin = 0 ){
		return "dedx_tof_" + chargeString(charge) + "_" + centerSpecies + "_" + ts(ptBin) + "_" + ts(etaBin);
	}
	string sTofName( string centerSpecies, int charge, int ptBin, int etaBin = 0, string eSpecies = "" ){
		if ( "" == eSpecies )
			return "tof_" + chargeString(charge) + "_" + centerSpecies + "_" + ts(ptBin) + "_" + ts(etaBin);
		else
			return "tof_" + chargeString(charge) + "_" + centerSpecies + "_" + ts(ptBin) + "_" + ts(etaBin) + "_" + eSpecies;
	}
	string sDedxName( string centerSpecies, int charge, int ptBin, int etaBin = 0, string eSpecies = "" ){
		if ( "" == eSpecies )
			return "dedx_" + chargeString(charge) + "_" + centerSpecies + "_" + ts(ptBin) + "_" + ts(etaBin);
		else
			return "dedx_" + chargeString(charge) + "_" + centerSpecies + "_" + ts(ptBin) + "_" + ts(etaBin) + "_" + eSpecies;
	}

	void autoViewport( 	string pType, double p, double * tofLow, double* tofHigh, double * dedxLow, double * dedxHigh, 
						double tofPadding = 1, double dedxPadding = 1, double tofScaledPadding = 0, double dedxScaledPadding = 0 );

	vector<string> otherSpecies( string center ){
		string species[] = {"Pi", "K", "P" };
		vector<string> res;
		for ( int i = 0; i < 3; i++ ){
			if ( species[ i ] != center )
				res.push_back( species[ i ] );
		}
		return res;
	}
	vector<double> enhanceTof( string center, vector<string> others, double p ){

		double cMean = tofGen->mean( p, eMass( center ) );
		
		vector<double> res;
		for ( int i = 0; i < others.size(); i++ ){
			double m = (tofGen->mean( p, eMass( others[ i ] ) ) - cMean) / tofPlotSigma;
			res.push_back( m );
		}

		return res;
	}

	vector<double> enhanceDedx( string center, vector<string> others, double p ){
		using namespace TMath;
		double cMean = Log10(dedxGen->mean( p, eMass( center ) ));
		
		vector<double> res;
		for ( int i = 0; i < others.size(); i++ ){
			double m = (Log10(dedxGen->mean( p, eMass( others[ i ] ) )) - cMean) / dedxPlotSigma;
			res.push_back( m );
		}

		return res;
	}

	/*
	*	Utility functions that should be moved soon
	*/ 
	void startTimer( ) { startTime = clock(); }
	double elapsed( ) { return ( (clock() - startTime) / (double)CLOCKS_PER_SEC ); }
};



#endif