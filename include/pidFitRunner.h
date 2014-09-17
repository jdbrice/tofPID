#ifndef PID_FIT_RUNNER
#define PID_FIT_RUNNER 

#include "reporter.h"
#include "xmlConfig.h"
#include "histoBook.h"
#include "jdbUtils.h"
using namespace jdbUtils;
// stl
#include "string"
using namespace std;
// root
#include "TH1D.h"
#include "TH2D.h"
#include "TH3.h"
#include "TFile.h"

class pidFitRunner
{
protected:
	/**
	 * Xml Configuration
	 * Ours and the one with info on how the data was produced
	 */
	xmlConfig * config;
	xmlConfig * prodConfig;

	/**
	 * histoBook for storing data
	 */
	histoBook * book;
	/**
	 * Main pdf Reporter
	 */
	reporter * report;

	/**
	 * Configuration parameter aliases
	 */
	string dataSource;
		static const string sourceData;
		static const string sourceSimulation;
		static const string sourceTruth;
	int centerSpecies;
	string tofAxis;
	string dedxAxis;
	vector<double> pBins;
	double pMin, pMax;

	/**
	 * Data file input
	 */
	TFile * dataFile;

	/**
	 * Histogram slices
	 */
	TH3 * hPid;
	TH2D* pSlice;
	TH1D* tofSlice;
	TH1D* dedxSlice;




public:
	pidFitRunner( xmlConfig * con );
	~pidFitRunner();

	/**
	 * Runs the fit jobs specified in the config file
	 * 
	 */
	void runFit();


protected:

	/**
	 * Produces the pSlice, tofSlice and dedxSlice histograms.
	 * After this point the system is agnostic of the dataSource
	 * meaning that data, simulated data, and species truth are all treated equally.
	 */
	void sliceHistograms( int pBin );

	void run1DFit( string fitTo = "tof" );

	class fitData {

	public:
		vector<double> mean;
		vector<double> meanError;
		vector<double> sigma;
		vector<double> sigmaError;
		vector<double> efficiency;
		vector<double> efficiencyError;
		vector<double> purity;
		vector<double> purityError;

		TH1D * histFrom( string name, vector<double> pBins ){

			int nBinsP = pBins.size();

			vector<double>  val = mean;
			vector<double> error = meanError;

			if ( "sigma" == name ){
				val = sigma;
				error = sigmaError;
			}

			TH1D* h = new TH1D( name.c_str(), name.c_str(), nBinsP-1, pBins.data() );
			for ( int i = 0; i < nBinsP; i++ ){
				
				if ( i >= val.size() ){
					//h->SetBinContent( i+1, 0);
					//h->SetBinError( i+1, 0);	
				} else {
					h->SetBinContent( i+1, val[ i ]);
					h->SetBinError( i+1, error[ i ]);
				}
				
				

			}

			return h;
		}

	};

	fitData fitResult;

	
};

#endif

















