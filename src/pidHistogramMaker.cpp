#include "constants.h"
#include "pidHistogramMaker.h"
#include "histoBook.h"

#include "TLine.h"

#include <fstream>
#include <sstream>

// provides my own string shortcuts etc.
using namespace jdbUtils;

const string pidHistogramMaker::inverseBeta = "inverseBeta";
const string pidHistogramMaker::deltaBeta = "deltaBeta";

const string pidHistogramMaker::traditionalCentering = "traditional";
const string pidHistogramMaker::nonlinearCentering = "nonlinear";

const vector<string> pidHistogramMaker::species = {"Pi","K","P"};

/**
 * Constructor - Initializes all of the pidHistogram parameters from the configuration file
 * @param chain       The chain object containing all data compatible with the TOFrPicoDST format
 * @param con         The xml configuration defining key aspects of the calibration
 *					such as number of tot bins to use, data location etc. See repo Readme
 *					for a sample configuration.
 */
pidHistogramMaker::pidHistogramMaker( TChain* chain, xmlConfig* con )  {
	cout << "[pidHistogramMaker.pidHistogramMaker] " << endl;
	
	gErrorIgnoreLevel=kSysError ;

	config = con;


	// set the histogram info verbosity to show nothing
	gStyle->SetOptStat( config->getInt( "statBox.show", 0 ) );
	gStyle->SetStatX( config->getDouble( "statBox.pos:x", 0.85 ) );
	gStyle->SetStatY( config->getDouble( "statBox.pos:y", 0.9 ) );
	gStyle->SetStatH( config->getDouble( "statBox.pos:h", 0.2 ) );
	gStyle->SetStatW( config->getDouble( "statBox.pos:w", 0.2 ) );
	
	// create the histogram book
	book = new histoBook( ( config->getString( "output.base" ) + config->getString( "output.root" ) ), config, config->getString( "input.root" ) );

	if ( true == config->nodeExists( "pType" )  ){
		vector<string> parts = config->getStringVector( "pType" );
		
		for ( int i = 0; i < parts.size(); i++ ){
			for ( int charge = -1; charge <= 1; charge ++ ){

				string n = speciesName( parts[ i ], charge );
				pReport[ n ] = new reporter( config->getString( "output.base" ) + n + config->getString( "output.report" ) );		
			}
		}
	}

	// create a report builder 
	report = new reporter( config->getString( "output.base" ) + config->getString( "output.report" ) );


	_chain = chain;
	pico = new TOFrPicoDst( _chain );


	vOffsetX = config->getDouble( "cut.vOffset:x", 0 );
	vOffsetY = config->getDouble( "cut.vOffset:y", 0 );

	inverseBetaSigma = config->getDouble( "centering.sigma:tof", 0.012 );
	dedxSigma = config->getDouble( "centering.sigma:dedx", 0.06        );      
	
	// for centering only
	tofGen = new tofGenerator( inverseBetaSigma );
	dedxGen = new dedxGenerator(  );

	centeringMethod = config->getString( "centering.mode", traditionalCentering );
	tofShift = config->getDouble( "centering.globalShift:tof", 0.0 );
	dedxShift = config->getDouble( "centering.globalShift:dedx", 0.0 );

	// for putting nice ranges on plots
	tofPadding = config->getDouble( "binning.tof:padding", 5 );
	dedxPadding = config->getDouble( "binning.dedx:padding", 5 );
	tofScalePadding = config->getDouble( "binning.tof:paddingScale", .05 );
	dedxScalePadding = config->getDouble( "binning.dedx:paddingScale", .05 );
}

/**
 *	Destructor - Deletes the histoBook ensuring it is saved.
 */
pidHistogramMaker::~pidHistogramMaker() {
	
	delete book;
	delete report;

	if ( true == config->nodeExists( "pType" )  ){
		vector<string> parts = config->getStringVector( "pType" );
		for ( int i = 0; i < parts.size(); i++ ){
			delete pReport[ speciesName( parts[ i ], -1 ) ];
			delete pReport[ speciesName( parts[ i ], 0 ) ];
			delete pReport[ speciesName( parts[ i ], 1 ) ];
		}
	}
	
	cout << "[pidHistogramMaker.~pidHistogramMaker] " << endl;
}

void pidHistogramMaker::makeQA() {


	

	startTimer();

	gStyle->SetOptStat( 11 );

	if ( !_chain ){
		cout << "[pidHistogramMaker." << __FUNCTION__ << "] ERROR: Invalid chain " << endl;
		return;
	}

	Int_t nevents = (Int_t)_chain->GetEntries();
	cout << "[pidHistogramMaker." << __FUNCTION__ << "] Loaded: " << nevents << " events " << endl;

	book->cd( "" );

	//book->make1D( "length", "length", 1000, 0, 500 );

	book->makeAll( "histo" );
	
	
	// loop over all events
	for(Int_t i=0; i<nevents; i++) {
    	_chain->GetEntry(i);

    	progressBar( i, nevents, 60 );

    	int nT0 = pico->nTZero;
    	int nTofHits = pico->nTofHits;
    	

    	double vz = pico->vertexZ;
    	double vx = pico->vertexX;
    	double vy = pico->vertexY;
    	double vr = TMath::Sqrt( vx*vx + vy*vy );
    	double vrOff = TMath::Sqrt( (vx - vOffsetX)*(vx - vOffsetX) + (vy - vOffsetY)*(vy - vOffsetY) );
    	
    	// nHits QA
    	book->fill( "nTZero", nT0 );
    	book->fill( "nTofHits", nTofHits );
    	//book->fill( "nHits", nHits );

    	// Vertex distributions before any cuts
    	book->fill( "preVz", vz );

    	book->fill( "preOffsetVr", vr );
    	book->fill( "preOffsetVxy", vx, vy);

    	// with the x,y brought back to 0
    	book->fill( "preVr", vrOff );
    	book->fill( "preVxy", vx - vOffsetX, vy - vOffsetY);

		if ( TMath::Abs( vz ) < config->getDouble( "cut.vZ", 30 ) ){
			book->fill( "postVz", vz );
			book->fill( "postVrCutVz", vrOff );
		}
		if (  vrOff < config->getDouble( "cut.vR", 10 ) ){
			book->fill( "postVr", vrOff );
			book->fill( "postVxy", vx - vOffsetX, vy - vOffsetY );
		}  
		if ( TMath::Abs( vz ) < config->getDouble( "cut.vZ", 30 ) && vrOff < config->getDouble( "cut.vR", 10 ) ){
			book->fill( "postVzCutVzVr", vz );	
			book->fill( "postVrCutVzVr", vrOff );	
		}  	
		if ( nT0 > config->getDouble( "cut.nT0", 0 ) ){
			book->fill( "postNTZero", nT0 );
		}
    	if ( nTofHits > config->getDouble( "cut.nTof", 0 ) ){
    		book->fill( "postNTofHits", nTofHits );
    	}

    	

    	if ( !keepEventQA() )
    		continue;


    	double tStart = pico->tStart;
    	Int_t refMult = pico->refMult;

    
    	for ( int iHit = 0; iHit < nTofHits; iHit++ ){


    		book->fill( "nHits", pico->nHits[ iHit ] );
    		book->fill( "nHitsFit", pico->nHitsFit[ iHit ] );
    		book->fill( "nHitsPossible", pico->nHitsPossible[ iHit ] );


    		if ( !keepTrackQA( iHit ) )
    			continue;

    		//if ( abs( pico->nSigPi[ iHit ] ) > 2 )
    		//	continue;
    		

    		double le = pico->leTime[ iHit ];
    		double length = pico->length[ iHit ];
    		double tof = pico->tof[ iHit ];
    		double beta = pico->beta[ iHit ];

    		double p = pico->p[ iHit ];
    		double dedx = pico->dedx[ iHit ];
    		
    		double m2 = p*p * ( constants::c*constants::c * tof*tof / ( length*length ) - 1  );

    		// mass cut?
    		//if ( TMath::Abs(m2 - constants::kaonMass*constants::kaonMass) > .1 )
    		//	continue;

    		double deltaB = 1 - beta * TMath::Sqrt( 1 - m2 / ( p*p )  );
    		//double deltaBPi = 1 - beta * TMath::Sqrt( 1 - (.139*.139) / ( p*p ) );

    		double deltaBPi = 1 - (beta) * TMath::Sqrt( (constants::piMass*constants::piMass) / (p*p) + 1 );
    		double deltaBK = 1 - (beta) * TMath::Sqrt( (constants::kaonMass*constants::kaonMass) / (p*p) + 1 );

    		double sigIBeta = 0.013;
    		double bnSigK = ((1.0/beta) - (1.0 / eBeta( constants::kaonMass, p ) )) / sigIBeta ;
    		double bnSigP = ((1.0/beta) - (1.0 / eBeta( constants::protonMass, p ) )) / sigIBeta ;
    		double bnSigPi = ((1.0/beta) - (1.0 / eBeta( constants::piMass, p ) )) / sigIBeta ;
    		
    		book->fill( "m2", m2 );
    		book->fill( "m2p", p, m2 );
    		book->fill( "m2dedx", TMath::Log(dedx), m2 );
    		
    		book->fill( "dedxP", p, dedx );
    		book->fill( "iBeta", p, (1.0/beta) );

    		book->fill( "dedxVsBeta", TMath::Exp( beta ), TMath::Log( dedx ) );

    		
    		book->fill( "deltaB", deltaB );
    		book->fill( "deltaBPi", deltaBPi );
    		book->fill( "deltaBSigPi", deltaB, pico->nSigPi[iHit] );

    		//book->fill( "deltaBVsP", p, deltaBKaon );

    		book->fill( "bnSigK", bnSigK );
    		book->fill( "bnSigKVsP", p, bnSigK );
    		book->fill( "bnSigPVsP", p, bnSigP );
    		book->fill( "bnSigPiVsP", p, bnSigPi );

    		book->fill( "nSigBetaDedxK", pico->nSigK[iHit], bnSigK );
    		book->fill( "nSigBetaDedxP", pico->nSigP[iHit], bnSigP );
    		book->fill( "nSigBetaDedxPi", pico->nSigPi[iHit], bnSigPi );

    		if ( p < 1.15 && p > 1.049 ){
    			book->fill( "deltaBSigPi", deltaBPi, pico->nSigPi[ iHit ] );
    			book->fill( "deltaBSigK", deltaBK, pico->nSigK[ iHit ] );
    		}

  

    	}
    	
	} // end loop on events

	report->newPage( 2, 2);
	book->style( "nTZero" )->draw();
	report->next();
	book->style( "nTofHits" )->draw();
	report->next();
	book->style( "postNTZero" )->draw();
	report->next();
	book->style( "postNTofHits" )->draw();
	report->savePage();

	report->newPage( 2, 2);
	book->style( "preVz" )->draw();
	report->next();
	book->style( "postVz" )->draw();
	report->next();
	book->style( "postVzCutVzVr" )->draw();
	report->savePage();

	report->newPage( 2, 2);
	book->style( "preVr" )->draw();
	report->next();
	book->style( "postVrCutVz" )->draw();
	report->next();
	book->style( "postVr" )->draw();
	report->next();
	book->style( "postVrCutVzVr" )->draw();
	report->savePage();

	report->newPage( 2, 2);
	gStyle->SetOptStat( 111 );
	book->style( "preOffsetVxy" )->draw();
	report->next();
	book->style( "preVxy" )->draw();
	report->next();
	book->style( "postVxy" )->draw();
	report->savePage();
	gStyle->SetOptStat( 11 );


	report->newPage( 2, 2 );
	book->style( "nHits" )->draw();
	report->next();
	book->style( "nHitsFit" )->draw();
	report->next();
	book->style( "nHitsPossible" )->draw();
	report->savePage();

	report->newPage(1, 2);
	book->style("iBeta")->draw();

	TGraph * g1 = inverseBetaGraph( constants::piMass, 0.15, 3, .05 );
	TGraph * g2 = inverseBetaGraph( constants::kaonMass, 0.15, 3, .05 );
	TGraph * g3 = inverseBetaGraph( constants::protonMass, 0.15, 3, .05 );
	TGraph * g4 = inverseBetaGraph( constants::eMass, 0.15, 3, .05 );
	g1->Draw( "same" );
	g2->Draw( "same" );
	g3->Draw( "same" );
	g4->Draw( "same" );
	report->savePage();

	report->newPage( 1, 2);
	book->style("dedxP")->draw();
	report->next();
	book->style("dedxVsBeta")->draw();
	report->savePage();

	report->newPage( 2, 2);
	book->style("m2")->draw();
	report->next();
	book->style("m2p")->draw();
	report->next();
	book->style("m2dedx")->draw();
	report->savePage();



	cout << "[pidHistogramMaker." << __FUNCTION__ << "] completed in " << elapsed() << " seconds " << endl;
}

bool pidHistogramMaker::keepEventQA(){


	double vz = pico->vertexZ;
	double vx = pico->vertexX;
	double vy = pico->vertexY;
	double vr = TMath::Sqrt( vx*vx + vy*vy );
	double vrOff = TMath::Sqrt( (vx - vOffsetX)*(vx - vOffsetX) + (vy - vOffsetY)*(vy - vOffsetY) );

	if ( TMath::Abs( vz ) > config->getDouble( "cut.vZ", 30 ) )
    	return false;
    if (  vrOff > config->getDouble( "cut.vR", 10 ) )
    	return false;

    int nT0 = pico->nTZero;
    int nTofHits = pico->nTofHits;

    if ( nT0 < config->getDouble( "cut.nT0", 0 ) )
    	return false;
    if ( nTofHits < config->getDouble( "cut.nTof", 0 ) )
    	return false;


    return true;
}

bool pidHistogramMaker::keepTrackQA( uint iHit ){

	double p = pico->p[ iHit ];
	double nHits = pico->nHits[ iHit ];
	double eta = pico->eta[ iHit ];

	if ( TMath::Abs( eta ) > .20 )
		return false;

	return true;
}


TGraph * pidHistogramMaker::inverseBetaGraph( double m, double p1, double p2, double step ){

	if ( p1 <= 0 )
		p1 = 0.01;

	int n = (p2 - p1) / step ;
	double * x = new double[ n + 2 ];
	double * y = new double[ n + 2 ];

	int i = 0;
	for ( double p = p1; p <= p2; p += step ){

		x[ i ] = p ;
		y[ i ] =  1.0 / eBeta( m, p );
		i++;

	}
	TGraph* g = new TGraph( i, x, y );
	g->SetLineColor( kRed );
	g->SetLineWidth( 2 );
	return  g ;

}

void pidHistogramMaker::momentumDistributions() {

	book->make( "histograms.momentum" );

	Int_t nevents = (Int_t)_chain->GetEntries();
	cout << "[pidHistogramMaker." << __FUNCTION__ << "] Loaded: " << nevents << " events " << endl;

	book->cd( "" );

	vector<string> parts = config->getStringVector( "pType" );

	for ( int i = 0; i < parts.size(); i++ ){
		prepareHistograms( parts[ i ] );
	}

	taskProgress tp( "Making Pid Histograms", nevents );
	// loop over all events to produce the histograms
	for(Int_t i=0; i<nevents; i++) {
    	_chain->GetEntry(i);

    	// report progress
    	tp.showProgress( i );

    	// use the QA event cuts
    	if ( !keepEventQA() )
    		continue;
    	
    	// Loop over all tof hits
    	int nTofHits = pico->nTofHits;
    	for ( int iHit = 0; iHit < nTofHits; iHit++ ){

    		// Use the QA Track cuts
    		if ( !keepTrackQA( iHit ) )
	    		continue;
  		
    		double pt = pico->pt[ iHit ];
    		double p = pico->pt[ iHit ];

			string pType = parts[ 0 ];
			
			// check the limits so we dont process more than we need to
			if ( pt > pMax || pt < pMin )
				continue;

			// get and fill the histogram for this pType and charge
			string name = "nSig_" + speciesName( pType, charge );
			TH3* h3 = ((TH3*)book->get( name ));
			
			int ptBin = h3->GetZaxis()->FindBin( pt );
			book->
		}
	}


}


void pidHistogramMaker::makePidHistograms() {

	if ( !_chain ){
		cout << "[pidHistogramMaker." << __FUNCTION__ << "] ERROR: Invalid chain " << endl;
		return;
	}

	Int_t nevents = (Int_t)_chain->GetEntries();
	cout << "[pidHistogramMaker." << __FUNCTION__ << "] Loaded: " << nevents << " events " << endl;

	book->cd( "" );

	vector<string> parts = config->getStringVector( "pType" );

	for ( int i = 0; i < parts.size(); i++ ){
		prepareHistograms( parts[ i ] );
	}

	taskProgress tp( "Making Pid Histograms", nevents );
	// loop over all events to produce the histograms
	for(Int_t i=0; i<nevents; i++) {
    	_chain->GetEntry(i);

    	// report progress
    	tp.showProgress( i );

    	// use the QA event cuts
    	if ( !keepEventQA() )
    		continue;
    	
    	// Loop over all tof hits
    	int nTofHits = pico->nTofHits;
    	for ( int iHit = 0; iHit < nTofHits; iHit++ ){

    		// Use the QA Track cuts
    		if ( !keepTrackQA( iHit ) )
	    		continue;
  		
    		double p = pico->pt[ iHit ];

    		// compute centered distributions for each particle type
    		for ( int i = 0; i < parts.size(); i++ ){
    			string pType = parts[ i ];
    			
    			// collect variables 
				int charge = pico->charge[ iHit ];
				double eta = pico->eta[ iHit ];
				double dedx = nSigDedx( pType, iHit ) - dedxShift;
				
				// Configuration allows you to choose between these two metrics
				double invBeta = nSigInvBeta( pType, iHit) - tofShift;
				double deltaB = dBeta( pType, iHit );
				double tof = invBeta;
				if ( tofMetric == deltaBeta )
					tof = deltaB;
				
				// check the limits so we dont process more than we need to
				if ( p > pMax || p < pMin )
					continue;
				if ( dedx >= dedxMax || dedx <= dedxMin )
					continue;
				if ( deltaB >= tofMax || deltaB <= tofMin )
					continue;

				

				// get and fill the histogram for this pType and charge
				string name = "nSig_" + speciesName( pType, charge );
				TH3* h3 = ((TH3*)book->get( name ));

				// get pBin;
				double avgP = 0;
				if ( h3 ){
					int pBin = h3->GetZaxis()->FindBin( p );
					if ( pBin > 0 && pBin < pBins.size() )
						avgP = (pBins[ pBin ] + pBins[ pBin - 1 ] ) / 2.0;
				}

				/**
				 * Switch centering methods
				 */
				if (  nonlinearCentering == centeringMethod ){
					dedx = nSigmaDedx( pType, iHit, avgP ) - dedxShift;
					tof = nSigmaInverseBeta( pType, iHit, avgP ) - tofShift;
				}

				if ( h3 ){
					h3->Fill( dedx, tof, p );
				}

				// always fill the charge agnostic histogram
				name = "nSig_" + speciesName( pType, 0 );
				h3 = ((TH3*)book->get( name ));
				if ( h3 ){
					h3->Fill( dedx, tof, p );
				}
			

			} // loop pTypes
    		
    	}
    	
	} // end loop on events

	// make particle type reports
	for ( int i = 0; i < parts.size(); i++ ){
		//speciesReport( parts[ i ], -1 );
		//speciesReport( parts[ i ], 0 );
		distributionReport( parts[ i ] );
		//speciesReport( parts[ i ], 1 );
	}

}

string pidHistogramMaker::speciesName( string pType, int charge ){

	if ( -1 == charge )
		return pType + "_Negative";
	if ( 1 == charge )
		return pType + "_Positive";
	if ( 0 == charge )
		return pType + "";
	return "";
}

void pidHistogramMaker::prepareHistograms( string pType ) {
	cout << "[pidHistogramMaker." << __FUNCTION__ << "] Center Species : " << pType << endl;
	/**
	 * Make the dedx binning 
	 */
	double dedxBinWidth = config->getDouble( "binning.dedx:binWidth" );
	dedxMin = config->getDouble( "binning.dedx:min" );
	dedxMax = config->getDouble( "binning.dedx:max" );

	for ( double i = dedxMin; i <= dedxMax; i += dedxBinWidth ){
		dedxBins.push_back( i );
	}
	
	cout << "\tDedx bins created" << endl;

	/**
	 * Make the Tof Binning
	 * could be for either 1/beta of delta 1/beta 
	 */
	tofMetric = config->getString( "binning.tofMetric", "inverseBeta" );
	double tofBinWidth = config->getDouble( "binning.tof:binWidth" );
	tofMin = config->getDouble( "binning.tof:min" );
	tofMax = config->getDouble( "binning.tof:max" );

	for ( double i = tofMin; i <= tofMax; i += tofBinWidth ){
		tofBins.push_back( i );
	}
	cout << "\tTof bins created" << endl;


	/**
	 * Make the momentum binning
	 */
	pMin = config->getDouble( "binning.p:min" );
	pMax = config->getDouble( "binning.p:max" );

	if ( config->nodeExists( "binning.pBins" ) && config->getDoubleVector( "binning.pBins" ).size() >= 2 ){
		pBins = config->getDoubleVector( "binning.pBins" );
			
		pMin = pBins[ 0 ];
		pMax = pBins[ pBins.size() - 1 ];
	} else {
		// build the pBins from the range and binWidth
		double pBinWidth = config->getDouble( "binning.p:binWidth", .05 );
		for ( double i = pMin; i <= pMax; i+= pBinWidth ){
			pBins.push_back( i );
		}
	}
	cout << "\tMomentum bins created" << endl;


	string title = "";

	if ( tofMetric == inverseBeta )
		title = "; n#sigma dedx; n#sigma #beta^{-1} ";
	else 
		title = "; n#sigma dedx; #Delta #beta^{-1} / #beta^{-1} ";
		
	
	
	// create a combined, plus, and minus
	for ( int charge = -1; charge <= 1; charge ++ ){

		string name = "nSig_" + speciesName( pType, charge ) ;

		TH3D * h3 = new TH3D( name.c_str(), title.c_str(), 
				dedxBins.size()-1, dedxBins.data(), 
				tofBins.size()-1, tofBins.data(),
				pBins.size()-1, pBins.data() );

		book->add( name, h3 );
	}
	


}


void pidHistogramMaker::speciesReport( string pType, int charge, int etaBin ){
	
	string name = speciesName( pType, charge );

	cout << "\tSpecies Report : " << name << endl;

	uint nBinsP = pBins.size();
	

	TH3 * h3 = book->get3D( "nSig_" + name );

	taskProgress tp( pType + " report", nBinsP );

	for ( uint i = 0; i < nBinsP; i ++ ){

		tp.showProgress( i );

		// momentum value used for finding nice range
		double p = pBins[ i ];


		pReport[ name ]->newPage( 2, 2 );
		pReport[ name ]->cd( 1, 2 );
		h3->GetZaxis()->SetRange( i, i );
		
		TH2* proj;
		TH1* proj1D;

		double pLow = h3->GetZaxis()->GetBinLowEdge( i );
		double pHi = h3->GetZaxis()->GetBinLowEdge( i + 1 );
		if ( 0 == i ){
			pLow = h3->GetZaxis()->GetBinLowEdge( 1 );
			pHi = h3->GetZaxis()->GetBinLowEdge( pBins.size()-1 );
		}
		

		string order = "xy";
		string tofAxis = config->getString( "binning.tofAxis", "x" );
		if ( tofAxis == "y" )
			order = "yx";
		

		proj = (TH2*)h3->Project3D( order.c_str() );
		
		if ( i > 0 && tofMetric == inverseBeta ){
			double tofLow, tofHigh, dedxLow, dedxHigh;
			autoViewport( pType, p, &tofLow, &tofHigh, &dedxLow, &dedxHigh, tofPadding, dedxPadding, tofScalePadding, dedxScalePadding );
		
			proj->GetYaxis()->SetRangeUser( tofLow, tofHigh );
			proj->GetXaxis()->SetRangeUser( dedxLow, dedxHigh );
		}
		
		string hTitle = (pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 ) );
		string hName = (pType + "_" + ts(i) + "_" + ts(etaBin) );
		proj->SetTitle( hTitle.c_str()  );
		gPad->SetLogz( 1 );
		proj->Draw( "colz" );

	
		if ( tofAxis == "x" ){
			pReport[ name ]->cd( 2, 2 );
			proj1D = proj->ProjectionX();
			proj1D->SetTitle( ( "dedx : " + pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 )  ).c_str()  );
			proj1D->SetFillColor( kBlue );
			gPad->SetLogx( 1 );
			proj1D->Draw( "hbar" );

			pReport[ name ]->cd( 1, 1 );
			proj1D = proj->ProjectionY();
			proj1D->SetTitle( ( "1/#beta : " + pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 )  ).c_str()  );
			gPad->SetLogy( 1 );
			proj1D->SetFillColor( kBlue );
			proj1D->Draw( "" );
		} else {
			pReport[ name ]->cd( 2, 2 );
			proj1D = proj->ProjectionY();
			proj1D->SetTitle( ( "#beta^{-1} : " + pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 )  ).c_str()  );
			proj1D->SetFillColor( kBlue );
			gPad->SetLogx( 1 );
			proj1D->Draw( "hbar" );

			pReport[ name ]->cd( 1, 1 );
			proj1D = proj->ProjectionX();
			proj1D->SetTitle( ( "dedx : " + pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 )  ).c_str()  );
			gPad->SetLogy( 1 );
			proj1D->SetFillColor( kBlue );
			proj1D->Draw( "" );
		}

		


		pReport[ name ]->savePage();

	}

	

} 

void pidHistogramMaker::distributionReport( string pType ){


	histoBook * dBook = new histoBook( pType + "data.root", config );

	string name = speciesName( pType, 0 );

	cout << "\tSpecies Report : " << name << endl;

	uint nBinsP = pBins.size();
	

	TH3 * h3 = book->get3D( "nSig_" + name );

	taskProgress tp( pType + " distribution report", nBinsP );

	for ( uint i = 0; i < nBinsP; i ++ ){

		tp.showProgress( i );

		// momentum value used for finding nice range
		double p = pBins[ i ];

		// start a new page on the report file
		pReport[ name ]->newPage( 2, 2 );

		// get the h3 and set the current p bin range
		h3->GetZaxis()->SetRange( i, i );	

		// get the 2D projection in dedx X tof space
		TH2* proj = (TH2*)h3->Project3D( "xy" );

		// get information on plot ranges
		double tofLow, tofHigh, dedxLow, dedxHigh;
		autoViewport( pType, p, &tofLow, &tofHigh, &dedxLow, &dedxHigh, tofPadding, dedxPadding, tofScalePadding, dedxScalePadding );
		
		// Make the all tof tracks histogram
		TH1* hTof = proj->ProjectionX();
		dBook->add( "tof_"+ts(i), (TH1*)hTof->Clone( ("tof_"+ts(i)).c_str()  ) );
		dBook->style( "tof_"+ts(i) )->set( "style.tof" ) 	
			 ->set("domain", tofLow, tofHigh )->draw();

		double yMax = hTof->GetMaximum();
		double piMean = (tofGen->mean( p, eMass( "Pi" ) ) - tofGen->mean( p, eMass( pType ) ) ) / inverseBetaSigma;
		TLine * l1 = new TLine( piMean, 0, piMean, yMax );
		l1->Draw();

		// Make the enhanced histogram for the center species
		pReport[ name ]->cd( 2, 1 );
		TH2* proj2 = (TH2*)proj->Clone( "proj2");
		int y1 = proj2->GetYaxis()->FindBin( -1 );
		int y2 = proj2->GetYaxis()->FindBin( 1 );

		// sets the rnge to effect the cut in dedx space
		proj2->GetYaxis()->SetRange( y1, y2 );
		TH1* hTofEnhanced = proj2->ProjectionX( "_px" );

		dBook->add( "tof_" + pType + ts(i), (TH1*)hTofEnhanced->Clone( ("tof_" + pType + ts(i)).c_str()  ) );
		dBook->style( "tof_" + pType + ts(i) )->set( "style.tof" )
			 ->set("domain", tofLow, tofHigh )->draw();
		

		pReport[ name ]->savePage();

	}

	delete dBook;


}

void pidHistogramMaker::autoViewport( 	string pType, 
										double p, double * tofLow, double* tofHigh, double * dedxLow, double * dedxHigh, 
										double tofPadding, double dedxPadding, double tofScaledPadding, double dedxScaledPadding  ){


	double tofCenter = tofGen->mean( p, eMass( pType ) );
	double dedxCenter = TMath::Log10( dedxGen->mean( p, eMass( pType ) ) );

	vector<double> tofMean;
	vector<double> dedxMean;

	for ( int i = 0; i < species.size(); i ++ ){
		tofMean.push_back(  tofGen->mean( p, eMass( species[ i ] ) ) );
		dedxMean.push_back( TMath::Log10( dedxGen->mean( p, eMass( species[ i ] ) ) ) );
	}

	double tHigh = (tofMean[ 0 ] - tofCenter);
	double tLow = (tofMean[ 0 ] - tofCenter);
	double dHigh = (dedxMean[ 0 ] - dedxCenter);
	double dLow = (dedxMean[ 0 ] - dedxCenter);
	for ( int i = 0; i < species.size(); i++ ){
		if ( (tofMean[ i ] - tofCenter) > tHigh )
			tHigh = (tofMean[ i ] - tofCenter);
		if ( (tofMean[ i ] - tofCenter) < tLow )
			tLow = (tofMean[ i ] - tofCenter);

		if ( (dedxMean[ i ] - dedxCenter) > dHigh )
			dHigh = (dedxMean[ i ] - dedxCenter);
		if ( (dedxMean[ i ] - dedxCenter) < dLow )
			dLow = (dedxMean[ i ] - dedxCenter);
	}
	
	


	*tofHigh = ( tHigh  / inverseBetaSigma) + tofPadding;
	*tofLow = ( tLow / inverseBetaSigma) - tofPadding;

	*dedxHigh = ( dHigh / dedxSigma) + dedxPadding;
	*dedxLow = ( dLow / dedxSigma) - dedxPadding;

	double tofRange = *tofHigh - *tofLow;
	double dedxRange = *dedxHigh - *dedxLow;

	double scaledPaddingTof = tofRange * tofScaledPadding;
	double scaledPaddingDedx = dedxRange * dedxScaledPadding;

	*tofLow -= scaledPaddingTof;
	*tofHigh += scaledPaddingTof;

	*dedxLow -= scaledPaddingDedx;
	*dedxHigh += scaledPaddingDedx;

	if ( *tofLow < tofMin )
		*tofLow = tofMin;
	if ( *tofHigh > tofMax )
		*tofHigh = tofMax;

	if ( *dedxLow < dedxMin )
		*dedxLow = dedxMin;
	if ( *dedxHigh > dedxMax )
		*dedxHigh = dedxMax;

	return;

}


double pidHistogramMaker::nSigDedx( string pType, int iHit ){ 
			
	double p = pico->p[ iHit ];
	double mean = TMath::Log10( dedxGen->mean( p, eMass( pType ) ) * 1000 );
	double dedx = TMath::Log10(pico->dedx[ iHit ]);
	double nSig = (( dedx - mean ) / dedxSigma );

	return nSig;

	return -999.0;
}

double pidHistogramMaker::nSigmaDedx( string pType, int iHit, double avgP ){

	double p = pico->p[ iHit ];
	double dedx = TMath::Log10( pico->dedx[ iHit ] );
	//cout << "\tdedx = " << dedx << endl;

	// mean for this species
	double mu = TMath::Log10( dedxGen->mean( p, eMass( pType ) ) * 1000 );
	double muAvg = TMath::Log10( dedxGen->mean( avgP, eMass( pType ) ) * 1000 );
	//cout << "\tmu = " << mu << endl;

	vector< string > species = { "K", "P", "Pi" };

	double n1 = 0, n2 = 0;
	double d1 = 0, d2 = 0;
	for ( int i = 0; i < species.size(); i++ ){

		double iMu = TMath::Log10( dedxGen->mean( p, eMass( species[ i ] ) ) * 1000 );
		double iMuAvg = TMath::Log10( dedxGen->mean( avgP, eMass( species[ i ] ) ) * 1000 );
		// may improve
		double sigma = dedxSigma; 

		double iL = lh( dedx, iMu, sigma );
		//cout << " iL " + species[ i ] + " " << iL << endl;
		double w = dedx + iMuAvg - iMu;
		

		n1 += (iL * w);
		d1 += iL;

		double iLc = lh( mu, iMu, sigma );
		n2 += (iLc * w);
		d2 += iLc;
	}

	//cout << "\tn1 = " << n1 << endl;
	//cout << "\td1 = " << d1 << endl;
	//cout << "\tn2 = " << n2 << endl;
	//cout << "\td2 = " << d2 << endl;

	double p1 = (n1/d1);
	double p2 = (n2/d2);
	//cout << "\tp1 = " << p1 << endl;
	//cout << "\tp2 = " << p2 << endl;
	double nSig = (n1/d1) - muAvg;

	return (nSig / dedxSigma);

}

double pidHistogramMaker::nSigmaInverseBeta( string pType, int iHit, double avgP ){

	double p = pico->p[ iHit ];
	double tof = 1.0 /  pico->beta[ iHit ];
	

	// mean for this species
	double mu =  tofGen->mean( p, eMass( pType ) );
	double muAvg =  tofGen->mean( avgP, eMass( pType ) );
	//cout << "\tmu = " << mu << endl;

	vector< string > species = { "K", "P", "Pi" };

	double n1 = 0, n2 = 0;
	double d1 = 0, d2 = 0;

	for ( int i = 0; i < species.size(); i++ ){

		double iMu =  tofGen->mean( p, eMass( species[ i ] ) ) ;
		double iMuAvg =  tofGen->mean( avgP, eMass( species[ i ] ) ) ;
		
		double sigma = inverseBetaSigma; 

		double iL = lh( tof, iMu, sigma );
		
		double w = tof + iMuAvg - iMu;
		
		n1 += (iL * w);
		d1 += iL;

		double iLc = lh( mu, iMu, sigma );
		n2 += (iLc * w);
		d2 += iLc;
	}



	double p1 = (n1/d1);
	double p2 = (n2/d2);
	
	double nSig = (n1/d1) - muAvg;

	return (nSig /  inverseBetaSigma);

}


double pidHistogramMaker::lh( double x, double mu, double sigma ){

	double a = sigma * TMath::Sqrt( 2 * TMath::Pi() );
	double b = ( x - mu );
	double c = 2 * sigma*sigma;
	double d = (1/a) * TMath::Exp( -b*b / c );

	return d;
}






