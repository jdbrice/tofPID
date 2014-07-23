
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

	book->make( "histo.m2p" );
	book->make( "histo.m2" );
	book->make( "histo.m2dedx" );
	book->make( "histo.iBeta" );
	book->make( "histo.dedxVsBeta" );
	book->make( "histo.dedxP" );
	book->make( "histo.deltaB" );
	book->make( "histo.deltaBPi" );
	book->make( "histo.deltaBSigPi" );
	
	// loop over all events
	for(Int_t i=0; i<nevents; i++) {
    	_chain->GetEntry(i);

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
    		
    		book->fill( "m2p", p, m2 );
    		book->fill( "m2dedx", TMath::Log(dedx), m2 );
    		book->fill( "m2", m2 );
    		book->fill( "iBeta", p, (1.0/beta) );
    		book->fill( "dedxVsBeta", TMath::Exp( beta ), TMath::Log( dedx ) );
    		book->fill( "dedxP", p, dedx );
    		book->fill( "deltaB", deltaB );
    		book->fill( "deltaBPi", deltaBPi );
    		book->fill( "deltaBSigPi", deltaB, pico->nSigPi[iHit] );

    	}
    	
	} // end loop on events

	report->newPage();
	book->style("m2p")->set( "style.log2D" )->draw();
	report->savePage();

	report->newPage();
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
	report->savePage();


	cout << "[pidHistogramMaker." << __FUNCTION__ << "] completed in " << elapsed() << " seconds " << endl;
}
