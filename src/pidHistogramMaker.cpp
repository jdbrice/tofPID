
#include "constants.h"
#include "pidHistogramMaker.h"
#include "histoBook.h"
#include <fstream>
#include <sstream>

// provides my own string shortcuts etc.
using namespace jdbUtils;


/**
 * Constructor - Initializes all of the pidHistogram parameters from the configuration file
 * @param chain       The chain object containing all data compatible with the TOFrPicoDST format
 * @param con         The xml configuration defining key aspects of the calibration
 *					such as number of tot bins to use, data location etc. See repo Readme
 *					for a sample configuration.
 */
pidHistogramMaker::pidHistogramMaker( TChain* chain, xmlConfig* con )  {
	cout << "[pidHistogramMaker.pidHistogramMaker] " << endl;
	
	gErrorIgnoreLevel=kError;

	config = con;


	// set the histogram info verbosity to show nothing
	gStyle->SetOptStat( 0 );
	
	// create the histogram book
	book = new histoBook( ( config->getString( "output.base" ) + config->getString( "output.root" ) ), config );
	


	vector<string> parts = config->getStringVector( "pType" );
	for ( int i = 0; i < parts.size(); i++ ){
		pReport[ parts[ i ] ] = new reporter( config->getString( "output.base" ) + parts[ i ] + config->getString( "output.report" ) );
	}

	// create a report builder 
	report = new reporter( config->getString( "output.base" ) + config->getString( "output.report" ) );


	_chain = chain;
	pico = new TOFrPicoDst( _chain );
}

/**
 *	Destructor - Deletes the histoBook ensuring it is saved.
 */
pidHistogramMaker::~pidHistogramMaker() {
	
	delete book;
	delete report;
	vector<string> parts = config->getStringVector( "pType" );
	for ( int i = 0; i < parts.size(); i++ ){
		delete pReport[ parts[ i ] ];
	}
	
	cout << "[pidHistogramMaker.~pidHistogramMaker] " << endl;
}



void pidHistogramMaker::loopEvents() {

	startTimer();

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

    	double vz = pico->vertexZ;
    	if ( TMath::Abs( vz ) > config->getDouble( "cut.vZ", 30 ) )
    		continue;

    	double tStart = pico->tStart;
    	Int_t refMult = pico->refMult;

    	

    	int nTofHits = pico->nTofHits;
    	for ( int iHit = 0; iHit < nTofHits; iHit++ ){

    		//if ( abs( pico->nSigPi[ iHit ] ) > 2 )
    		//	continue;

    		double le = pico->leTime[ iHit ];
    		double length = pico->length[ iHit ];
    		double tof = pico->tof[ iHit ];
    		double beta = pico->beta[ iHit ];

    		double p = pico->p[ iHit ];
    		double dedx = pico->dedx[ iHit ];
    		
    		double m2 = p*p * ( constants::c*constants::c * tof*tof / ( length*length ) - 1  );
    		double deltaB = 1 - beta * TMath::Sqrt( 1 - m2 / ( p*p )  );
    		double deltaBPi = 1 - beta * TMath::Sqrt( 1 - (.139*.139) / ( p*p ) );

    		double deltaBKaon = 1 - (beta) * TMath::Sqrt( 1 - (.009*.009) / ( p*p ) );

    		double sigIBeta = 0.013;
    		double bnSigK = ((1.0/beta) - (1.0 / eBeta( constants::kaonMass, p ) )) / sigIBeta ;
    		double bnSigP = ((1.0/beta) - (1.0 / eBeta( constants::protonMass, p ) )) / sigIBeta ;
    		double bnSigPi = ((1.0/beta) - (1.0 / eBeta( constants::piMass, p ) )) / sigIBeta ;
    		
    		book->fill( "m2p", p, m2 );
    		book->fill( "m2dedx", TMath::Log(dedx), m2 );
    		book->fill( "m2", m2 );
    		book->fill( "iBeta", p, (1.0/beta) );
    		book->fill( "dedxVsBeta", TMath::Exp( beta ), TMath::Log( dedx ) );
    		book->fill( "dedxP", p, dedx );
    		book->fill( "deltaB", deltaB );
    		book->fill( "deltaBPi", deltaBPi );
    		book->fill( "deltaBSigPi", deltaB, pico->nSigPi[iHit] );

    		book->fill( "deltaBVsP", p, deltaBKaon );

    		book->fill( "bnSigK", bnSigK );
    		book->fill( "bnSigKVsP", p, bnSigK );
    		book->fill( "bnSigPVsP", p, bnSigP );
    		book->fill( "bnSigPiVsP", p, bnSigPi );

    		book->fill( "nSigBetaDedxK", pico->nSigK[iHit], bnSigK );
    		book->fill( "nSigBetaDedxP", pico->nSigP[iHit], bnSigP );
    		book->fill( "nSigBetaDedxPi", pico->nSigPi[iHit], bnSigPi );

    	}
    	
	} // end loop on events

	report->newPage();
	book->style("iBeta")->set( "style.log2D" )->draw();

	TGraph * g1 = inverseBeta( constants::piMass, 0.15, 3, .05 );
	TGraph * g2 = inverseBeta( constants::kaonMass, 0.15, 3, .05 );
	TGraph * g3 = inverseBeta( constants::protonMass, 0.15, 3, .05 );
	g1->Draw( "same" );
	g2->Draw( "same" );
	g3->Draw( "same" );
	report->savePage();

	report->newPage();
	book->style("bnSigK")->draw();
	report->savePage();

	report->newPage();
	book->style("bnSigKVsP")->set( "style.log2D" )->draw();
	report->savePage();

	report->newPage();
	book->style("bnSigPVsP")->set( "style.log2D" )->draw();
	report->savePage();

	report->newPage();
	book->style("bnSigPiVsP")->set( "style.log2D" )->draw();
	report->savePage();

	report->newPage();
	book->style("nSigBetaDedxK")->set( "style.log2D" )->draw();
	report->savePage();

	report->newPage();
	book->style("nSigBetaDedxP")->set( "style.log2D" )->draw();
	report->savePage();

	report->newPage();
	book->style("nSigBetaDedxPi")->set( "style.log2D" )->draw();
	report->savePage();

	


	/*report->newPage();
	book->style("iBeta")->set( "style.log2D" )->draw();
	report->savePage();

	report->newPage();
	book->style("m2dedx")->set( "style.log2D" )->draw();
	report->savePage();

	report->newPage();
	book->style("dedxVsBeta")->set( "style.log2D" )->draw();
	report->savePage();

	report->newPage();
	book->style("dedxP")->set( "style.log2D" )->draw();
	report->savePage();

	report->newPage();
	book->style("deltaBSigPi")->set( "style.log2D" )->draw();
	report->savePage();

	report->newPage();
	book->style("m2")->set( "style.log1D" )->draw();
	report->savePage();*/


	cout << "[pidHistogramMaker." << __FUNCTION__ << "] completed in " << elapsed() << " seconds " << endl;
}



TGraph * pidHistogramMaker::inverseBeta( double m, double p1, double p2, double step ){

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


void pidHistogramMaker::make() {

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
		sHisto( parts[ i ] );
	}

	book->make( "nSigBetaDedxK" );


	
	// loop over all events
	for(Int_t i=0; i<nevents; i++) {
    	_chain->GetEntry(i);

    	double vz = pico->vertexZ;
    	if ( TMath::Abs( vz ) > config->getDouble( "cut.vZ", 30 ) )
    		continue;
    	

    	int nTofHits = pico->nTofHits;
    	for ( int iHit = 0; iHit < nTofHits; iHit++ ){
  		
    		double p = pico->p[ iHit ];

    		for ( int i = 0; i < parts.size(); i++ ){
    			string pType = parts[ i ];
    			
				int charge = pico->charge[ iHit ];
				double dedx = nSigDedx( pType, iHit );
				double invBeta = nSigInvBeta( pType, iHit);
				
				if ( "Pi" == pType ){
					//cout << "nSigPi : " << pico->nSigPi[ iHit ] << endl ;
					//cout << "nSigPi : " << dedx << endl ;
				}

				//if ( dedx < 0 )
				//	cout << "dedx " << endl;
				//if ( invBeta < 0 )
				//	cout << "beta " << endl;


				string name = "nSig_" + sName( pType, charge );
				TH3* h3 = ((TH3*)book->get( name ));
				if ( h3 )
					h3->Fill( dedx, invBeta, p );

				// always fill the both charges histo
				name = "nSig_" + sName( pType, 0 );
				h3 = ((TH3*)book->get( name ));
				if ( h3 )
					h3->Fill( dedx, invBeta, p );

				if ( "K" == pType ){
					book->fill( "nSigBetaDedxK", dedx, invBeta );

				}
			}
    		
    	}
    	
	} // end loop on events

	// make particle type reports
	for ( int i = 0; i < parts.size(); i++ ){
		speciesReport( parts[ i ] );
	}


	cout << "[pidHistogramMaker." << __FUNCTION__ << "] completed in " << elapsed() << " seconds " << endl;
}

string pidHistogramMaker::sName( string pType, int charge ){

	if ( -1 == charge )
		return pType + "_Negative";
	if ( 1 == charge )
		return pType + "_Positive";
	if ( 0 == charge )
		return pType + "_All";
	return "";
}

void pidHistogramMaker::sHisto( string pType ) {


	int nSigBins = config->getDouble( "binning.nSig:nBins" );
	double nSigMin = config->getDouble( "binning.nSig:min" );
	double nSigMax = config->getDouble( "binning.nSig:max" );

	int nPBins = config->getDouble( "binning.p:nBins" );
	double pMin = config->getDouble( "binning.p:min" );
	double pMax = config->getDouble( "binning.p:max" );





	// make the nsigma bins
	double step = (nSigMax - nSigMin ) / nSigBins;
	vector<double>nSigBinEdges;
	for ( double i = nSigMin; i < nSigMax; i += step ){
		nSigBinEdges.push_back( i );
	}
	nSigBinEdges.push_back( nSigMax );


	//vector<double> ppBins = { 0, .5, 1.0, 1.5, 2.0, 2.25, 3, 4 };

	vector<double> pBins = config->getDoubleVector( "binning.pBins" );

	// create a combined, plus, and minus
	for ( int charge = -1; charge <= 1; charge ++ ){
		string name = "nSig_" + sName( pType, charge );
		TH3D * h3 = new TH3D( name.c_str(), name.c_str(), 
				nSigBinEdges.size()-1, nSigBinEdges.data(), 
				nSigBinEdges.size()-1, nSigBinEdges.data(),
				pBins.size()-1, pBins.data() );

		book->add( name, h3 );
	}


}


double pidHistogramMaker::nSigInvBeta( string pType, int iHit  ){

	double b = pico->beta[ iHit ];
	double p = pico->p[ iHit ];
	double expectedBeta = eBeta( eMass( pType ), p );
	double invBetaSig = config->getDouble( "binning.invBetaSig" );

	double deltaInvBeta = ( 1.0 / b ) - ( 1.0 / expectedBeta );

	return (deltaInvBeta / invBetaSig);
}

void pidHistogramMaker::speciesReport( string pType ){


	vector<double>pBins = config->getDoubleVector( "binning.pBins" );

	TH3 * h3 = book->get3D( "nSig_" + sName( pType, 0 ) );
	for ( int i = 0; i < pBins.size(); i ++ ){

		pReport[ pType ]->newPage();
		h3->GetZaxis()->SetRange( i, i );
		TH2* proj = (TH2*)h3->Project3D( "xy" );


		double pLow = h3->GetZaxis()->GetBinLowEdge( i );
		double pHi = h3->GetZaxis()->GetBinLowEdge( i + 1 );
		if ( 0 == i ){
			pLow = h3->GetZaxis()->GetBinLowEdge( 1 );
			pHi = h3->GetZaxis()->GetBinLowEdge( pBins.size()-1 );
		}

		proj->SetTitle( (pType + " : " + ts( pLow ) + "#leq" + "P #leq" + ts( pHi )).c_str()  );
		proj->Draw( "colz" );

		pReport[ pType ]->savePage();

	}

}









