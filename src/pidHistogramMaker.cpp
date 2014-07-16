
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
	book = new histoBook( ( config->getString( "output.base" ) + config->getString( "output.root" ) ) );
	
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

	book->make( config, "histo.m2p" );
	book->make( config, "histo.m2" );
	
	// loop over all events
	for(Int_t i=0; i<nevents; i++) {
    	_chain->GetEntry(i);

    	double tStart = pico->tStart;
    	Int_t refMult = pico->refMult;

    	int nTofHits = pico->nTofHits;
    	for ( int iHit = 0; iHit < nTofHits; iHit++ ){
    		double le = pico->leTime[ iHit ];
    		double length = pico->length[ iHit ];
    		double tof = pico->tof[ iHit ];

    		double p = pico->p[ iHit ];

    		
    		double m2 = p*p * ( constants::c*constants::c * tof*tof / ( length*length ) - 1  );
    		book->fill( "m2p", p, m2 );
    		book->fill( "m2", m2 );
    		

    	}
    	
	} // end loop on events

	report->newPage();
	book->style("m2p")->set( config, "style.mass2" )->draw();
	report->savePage();

	report->newPage();
	gPad->SetLogy( 1 );
	book->style("m2")->set( "draw", config->getString("histo.m2:draw"))->draw();
	report->savePage();


	cout << "[pidHistogramMaker." << __FUNCTION__ << "] completed in " << elapsed() << " seconds " << endl;
}
