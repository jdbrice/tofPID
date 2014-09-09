
#include "constants.h"
#include "pidHistogramMaker.h"
#include "histoBook.h"
#include "dklMinimizer.h"
#include <fstream>
#include <sstream>

// provides my own string shortcuts etc.
using namespace jdbUtils;

const string pidHistogramMaker::inverseBeta = "inverseBeta";
const string pidHistogramMaker::deltaBeta = "deltaBeta";

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

	inverseBetaSigma = config->getDouble( "binning.invBetaSigma", 0.012 );
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




void pidHistogramMaker::makePidHistograms() {

	startTimer();

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

	// loop over all events to produce the histograms
	for(Int_t i=0; i<nevents; i++) {
    	_chain->GetEntry(i);

    	// report progress
    	progressBar( i, nevents, 60 );

    	// use the QA event cuts
    	if ( !keepEventQA() )
    		continue;
    	
    	// Loop over all tof hits
    	int nTofHits = pico->nTofHits;
    	for ( int iHit = 0; iHit < nTofHits; iHit++ ){

    		// Use the QA Track cuts
    		if ( !keepTrackQA( iHit ) )
	    		continue;
  		
    		double p = pico->p[ iHit ];

    		// compute centered distributions for each particle type
    		for ( int i = 0; i < parts.size(); i++ ){
    			string pType = parts[ i ];
    			
    			// collect variables 
				int charge = pico->charge[ iHit ];
				double eta = pico->eta[ iHit ];
				double dedx = nSigDedx( pType, iHit );
				
				// Configuration allows you to choose between these two metrics
				double invBeta = nSigInvBeta( pType, iHit);
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
		speciesReport( parts[ i ], 0 );
		//speciesReport( parts[ i ], 1 );
	}


	cout << "[pidHistogramMaker." << __FUNCTION__ << "] completed in " << elapsed() << " seconds " << endl;
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
		title = "; n#sigma dedx; #Delta #beta^{-1} / #beta^{-1} ";
	else 
		title = "; n#sigma dedx; n#sigma#beta^{-1} ";
	
	
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

	for ( uint i = 0; i < nBinsP; i ++ ){

		progressBar( i, nBinsP, 60 );

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

		proj = (TH2*)h3->Project3D( "xy" );
		string hTitle = (pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 ) );
		string hName = (pType + "_" + ts(i) + "_" + ts(etaBin) );
		proj->SetTitle( hTitle.c_str()  );
		gPad->SetLogz( 1 );
		proj->Draw( "colz" );

		
		pReport[ name ]->cd( 2, 2 );
		proj1D = h3->Project3D( "x" );
		proj1D->SetTitle( ( "dedx : " + pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 )  ).c_str()  );
		proj1D->SetFillColor( kBlue );
		gPad->SetLogx( 1 );
		proj1D->Draw( "hbar" );

		
		
		pReport[ name ]->cd( 1, 1 );
		proj1D = h3->Project3D( "y" );
		proj1D->SetTitle( ( "1/#beta : " + pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 )  ).c_str()  );
		gPad->SetLogy( 1 );
		proj1D->SetFillColor( kBlue );
		proj1D->Draw( "" );


		pReport[ name ]->savePage();
		//pReport[ name ]->saveImage( "report/png/" + hName + ".png" );

	}

	

}










