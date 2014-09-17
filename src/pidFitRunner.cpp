

#include "pidFitRunner.h"
#include "gaussianFitter.h"

const string pidFitRunner::sourceData = "data";
const string pidFitRunner::sourceSimulation = "simulation";
const string pidFitRunner::sourceTruth = "truth";

pidFitRunner::pidFitRunner( xmlConfig * con ){

	/**
	 * Removes the annoyance of root
	 */
	gErrorIgnoreLevel=kSysError ;

	/**
	 * Initialize the configuration files
	 */
	config = con;
	prodConfig = new xmlConfig( config->getString( "input.productionConfig" ).c_str() );

	/**
	 * Setup the histoBook
	 */
	book = new histoBook( config->getString( "output.dataUrl", "fitData.root"), config );
	book->cd();

	/**
	 * Open the reports
	 */
	report = new reporter( config->getString( "output.reportUrl", "fitQA.pdf" ) );
	

	/**
	 * Input data type and file
	 */
	dataSource = config->getString( "input.dataSource", sourceData );
	dataFile = new TFile( config->getString( "input.url" ).c_str(), "READ" );



	/**
	 * Data only initialization
	 */
	hPid = 0;
	if ( dataSource == sourceData ){
		hPid = (TH3*)dataFile->Get( config->getString( "h3Name", "nSig_K").c_str() );	

		pMin = prodConfig->getDouble( "binning.p:min" );
		pMax = prodConfig->getDouble( "binning.p:max" );

		if ( prodConfig->nodeExists( "binning.pBins" ) && prodConfig->getDoubleVector( "binning.pBins" ).size() >= 2 ){
			pBins = prodConfig->getDoubleVector( "binning.pBins" );
				
			pMin = pBins[ 0 ];
			pMax = pBins[ pBins.size() - 1 ];
		} else {
			// build the pBins from the range and binWidth
			double pBinWidth = prodConfig->getDouble( "binning.p:binWidth", .05 );
			for ( double i = pMin; i <= pMax; i+= pBinWidth ){
				pBins.push_back( i );
			}
		}
	} else if ( sourceSimulation == dataSource || sourceTruth == dataSource ){

	}

	tofAxis = prodConfig->getString( "binning.tofAxis", "y" );
	if ( "y" == tofAxis)
		dedxAxis = "x";
	else 
		dedxAxis = "y";

}

pidFitRunner::~pidFitRunner(){
	delete prodConfig;
	delete report;
	delete book;
}


void pidFitRunner::sliceHistograms( int pBin ){
	
	if ( dataSource == sourceData && hPid ){
		
		hPid->GetZaxis()->SetRange( pBin, pBin );
		
		pSlice = (TH2D*)hPid->Project3D( "xy" );
		tofSlice = (TH1D*)hPid->Project3D( tofAxis.c_str() );
		dedxSlice = (TH1D*)hPid->Project3D( dedxAxis.c_str() );

	}

}

void pidFitRunner::runFit(){
	cout << "[pidFitRunner." << __FUNCTION__ << "]" << endl;

	if ( config->nodeExists( "fits" ) ){
		vector<string> fitNames = config->getStringVector( "fits" );

		for ( int i = 0; i < fitNames.size(); i++ ){
			run1DFit( fitNames[ i ] );
		}
	}
	//run1DFit( "tofFit" );
	//run1DFit( "dedxFit" );


}


void pidFitRunner::run1DFit( string nodePath ) {

	string np = nodePath + ".";
	

	/**
	 * Get configuration
	 */
	string fitTo = config->getString( np + "fitTo", "tof" );
	double vpThreshold = config->getDouble( np + "viewport:threshold", 1 );
	double vpMin = config->getDouble( np + "viewport:min", -10 );
	double vpMax = config->getDouble( np + "viewport:max", 10 );

	double meanVal = config->getDouble( np + "mean:value", 0 );
	double meanMin = config->getDouble( np + "mean:min", -3 );
	double meanMax = config->getDouble( np + "mean:max", 3 );

	double sigmaVal = config->getDouble( np + "sigma:value", 0 );
	double sigmaMin = config->getDouble( np + "sigma:min", -3 );
	double sigmaMax = config->getDouble( np + "sigma:max", 3 );

	double roiValue = config->getDouble( np + "roi:value", 10);
	double roiStep = config->getDouble( np + "roi:step", 0);
	/**
	 * Get configuration
	 */

	gaussianFitter * gf = new gaussianFitter();

	int iMinP = config->getInt( np + "p:min" );
	int iMaxP = config->getInt( np + "p:max" );

	taskProgress tp( fitTo + " 1D fit", (iMaxP - iMinP) );

	for ( int ip = iMinP; ip < iMaxP; ip++ ){
		tp.showProgress( ip );	

		report->newPage();
		sliceHistograms( ip + 1 );


		TH1D * h = tofSlice;
		if ( "dedx" == fitTo )
			h = dedxSlice; 

		double b1 = h->FindFirstBinAbove( vpThreshold, 1 );
		double b2 = h->FindLastBinAbove( vpThreshold, 1 );
		double x1 = h->GetBinLowEdge( b1 );
		double x2 = h->GetBinLowEdge( b2 ) + h->GetBinWidth( b2 );

		if ( x1 > vpMin )
			x1 = vpMin;
		if ( x2 < vpMax )
			x2 = vpMax;
		
		string xName = "n#sigma 1/#beta";
		if ( "dedx" == fitTo )
			xName = "n#sigma dedx";

		gf->makeX( x1, x2, xName.c_str() );
		gf->makeMean( meanVal, meanMin, meanMax);
		gf->makeSigma( sigmaVal, sigmaMin, sigmaMax );
		
		double roi = roiValue;
		roi += ( ip * roiStep );
		
		gf->setROI( -roi, roi );

		gf->fitTo( h ); 

		gf->drawFit( );
		gPad->SetLogy( 1 );
		
		report->savePage();

		/**
		 * Store the fit results
		 */
		fitResult.mean.push_back( gf->mean->getVal() );
		fitResult.meanError.push_back( gf->mean->getError() );
		fitResult.sigma.push_back( gf->sigma->getVal() );
		fitResult.sigmaError.push_back( gf->sigma->getError() );

	}

	report->newPage();

	gPad->SetLogy( 0 );
	TH1D* hMean = fitResult.histFrom( "mean", pBins );
	hMean->Draw();
	report->savePage();

	report->newPage();	
	TH1D* hSigma = fitResult.histFrom( "sigma", pBins );
	hSigma->Draw();    

	report->savePage();
	

	delete gf;
 

}













