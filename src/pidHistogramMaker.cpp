#include "constants.h"
#include "pidHistogramMaker.h"


#include "TLine.h"

#include <fstream>
#include <sstream>


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
pidHistogramMaker::pidHistogramMaker( TChain* chain, XmlConfig* con )  {
	

	lg = LoggerConfig::makeLogger( con, "Logger" );

	lg->info(__FUNCTION__) << endl;
	
	gErrorIgnoreLevel=kSysError ;

	config = con;


	// set the histogram info verbosity to show nothing
	gStyle->SetOptStat( config->getInt( "statBox.show", 0 ) );
	gStyle->SetStatX( config->getDouble( "statBox.pos:x", 0.85 ) );
	gStyle->SetStatY( config->getDouble( "statBox.pos:y", 0.9 ) );
	gStyle->SetStatH( config->getDouble( "statBox.pos:h", 0.2 ) );
	gStyle->SetStatW( config->getDouble( "statBox.pos:w", 0.2 ) );
	

	distroData = new TFile( "distros.root", "RECREATE" );
	// create the histogram book
	book = new HistoBook( ( config->getString( "output.base" ) + config->getString( "output.root" ) ), config, config->getString( "input.root" ) );

	if ( true == config->nodeExists( "pType" )  ){
		vector<string> parts = config->getStringVector( "pType" );
		
		for ( int i = 0; i < parts.size(); i++ ){
			for ( int charge = -1; charge <= 1; charge ++ ){

				string n = speciesName( parts[ i ], charge );
				pReport[ n ] = new Reporter( config->getString( "output.base" ) + n + config->getString( "output.report" ) );		
			}
		}
	}

	// create a report builder 
	report = new Reporter( config->getString( "output.base" ) + config->getString( "output.report" ) );


	_chain = chain;
	pico = new TOFrPicoDst( _chain );


	vOffsetX = config->getDouble( "cut.vOffset:x", 0 );
	vOffsetY = config->getDouble( "cut.vOffset:y", 0 );

	tofSigma = config->getDouble( "centering.sigma:tof", 0.012 );
	dedxSigma = config->getDouble( "centering.sigma:dedx", 0.06 );

	tofPlotSigma = config->getDouble( "centering.plotSigma:tof", 0.012 );
	dedxPlotSigma = config->getDouble( "centering.plotSigma:dedx", 0.06 );
	
	// for centering only
	tofGen = new tofGenerator( tofSigma );
	dedxGen = new Bichsel( 	config->getString( "bichsel.table", "dedxBichsel.root"),
							config->getInt( "bichsel.method", 0 ) );

	centeringMethod = config->getString( "centering.mode", traditionalCentering );
	tofShift = config->getDouble( "centering.globalShift:tof", 0.0 );
	dedxShift = config->getDouble( "centering.globalShift:dedx", 0.0 );

	// for putting nice ranges on plots
	tofPadding = config->getDouble( "binning.tof:padding", 5 );
	dedxPadding = config->getDouble( "binning.dedx:padding", 5 );
	tofScalePadding = config->getDouble( "binning.tof:paddingScale", .05 );
	dedxScalePadding = config->getDouble( "binning.dedx:paddingScale", .05 );

	histosReady = false;

}

/**
 *	Destructor - Deletes the histoBook ensuring it is saved.
 */
pidHistogramMaker::~pidHistogramMaker() {
	lg->info(__FUNCTION__) << endl;
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
	
	lg->info(__FUNCTION__) << "Freed Memory" << endl;
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


    	//double tStart = pico->tStart;
    	//Int_t refMult = pico->refMult;

    
    	for ( int iHit = 0; iHit < nTofHits; iHit++ ){


    		book->fill( "nHits", pico->nHits[ iHit ] );
    		book->fill( "nHitsFit", pico->nHitsFit[ iHit ] );
    		book->fill( "nHitsPossible", pico->nHitsPossible[ iHit ] );


    		if ( !keepTrackQA( iHit ) )
    			continue;

    		//if ( abs( pico->nSigPi[ iHit ] ) > 2 )
    		//	continue;
    		

    		//double le = pico->leTime[ iHit ];
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
	//double vr = TMath::Sqrt( vx*vx + vy*vy );
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

	//double p = pico->p[ iHit ];
	//double nHits = pico->nHits[ iHit ];
	double eta = pico->eta[ iHit ];

	//if ( TMath::Abs( eta ) > .20 )
	//	return false;

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

	vector<string> plcs = config->getStringVector( "centering.species" );

	if ( false == histosReady ){
		for ( int i = 0; i < plcs.size(); i++ ){
			prepareHistograms( plcs[ i ] );
		}
	}

	book->makeAll( "histograms" );

	for ( int i = 0; i < etaBins.size() - 1; i ++ ){
		book->clone( "momentum", "momentum_"+ts(i) );
		book->clone( "pVsPt", "pVsPt_"+ts(i) );
	}


	

	Int_t nevents = (Int_t)_chain->GetEntries();
	cout << "[pidHistogramMaker." << __FUNCTION__ << "] Loaded: " << nevents << " events " << endl;

	

	taskProgress tp( "Filling Momentum distributions", nevents );
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

    		double eta = pico->eta[ iHit ];
    		double pt = pico->pt[ iHit ];
    		double p = pico->p[ iHit ];

    		double dedx = pico->dedx[ iHit ];
    		
    		int etaBin = HistoBook::findBin( etaBins, eta );
    		if ( etaBin < 0 )
    			continue;

			string pType = plcs[ 0 ];
			
			book->fill( "eta", eta );
			book->fill( "momentum" , p );
			book->fill( "momentum_" + ts(etaBin) , p );
			book->fill( "momentumTransverse", pt );
    		
    		// Use the QA Track cuts
    		if ( !keepTrackQA( iHit ) )
	    		continue;
	    	
  			
			book->fill( "pVsPt", pt, p);
			book->fill( "pVsPt_"+ts(etaBin), pt, p);

			book->fill( "dedxVsP", p, dedx );

		}
	}

	// loop through the pt bins
	for ( int ip = 0; ip < ptBins.size()-1; ip++ ){
		TH2* h2 = book->get2D( "pVsPt" );
		if ( h2 ){
			h2->GetXaxis()->SetRange( ip+1, ip+1 );
			double avgP = h2->GetMean( 2 );
			//cout << " <p[ " << ip << " ]> = " << avgP << endl; 
			averageP.push_back( avgP );
		}
	}
	

}


void pidHistogramMaker::makeDedxTofHistograms() {

	using namespace TMath;

	if ( !_chain ){
		lg->error( __FUNCTION__) << " Invalid chain " << endl;
		return;
	}

	Int_t nEvents = (Int_t)_chain->GetEntries();
	lg->info(__FUNCTION__) << "Loaded: " << nEvents << endl;

	

	vector<string> plcs = config->getStringVector( "centering.species" );

	if ( false == histosReady ){
		for ( int i = 0; i < plcs.size(); i++ ){
			prepareHistograms( plcs[ i ] );
		}
	}

	taskProgress tp( "Making Dedx vs. Tof Histograms", nEvents );
	
	book->cd( "dedx_tof" );

	// loop over all events to produce the histograms
	for(Int_t i=0; i<nEvents; i++) {
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

    		// compute centered distributions for each particle type
    		for ( int i = 0; i < plcs.size(); i++ ){
    			string pType = plcs[ i ];
    			
    			// collect variables 
				int charge = pico->charge[ iHit ];
				double eta = pico->eta[ iHit ];


				double dedx = nSigmaDedx( pType, iHit ) - dedxShift;
				double tof = nSigmaInverseBeta( pType, iHit) - tofShift;
				
				// check the limits so we dont process more than we need to
				if ( pt > ptMax || pt < ptMin )
					continue;
				if ( eta > etaMax || eta < etaMin )
					continue;

				/**
				 * Get the bin for the current kinematic values
				 */
				int ptBin = HistoBook::findBin( ptBins, pt );
				int etaBin = HistoBook::findBin( etaBins, Abs( eta ) );
				int chargeBin = HistoBook::findBin( chargeBins, charge );
				

				//double avgP = averageP[ ptBin ];
				double avgP = (ptBins[ ptBin ] + ptBins[ ptBin+1 ] ) / 2.0;

				/**
				 * Switch centering methods
				 */
				if (  nonlinearCentering == centeringMethod ){
					dedx = nSigmaDedx( pType, iHit, avgP ) - dedxShift;
					tof = nSigmaInverseBeta( pType, iHit, avgP ) - tofShift;
				}


				
				if ( chargeMin <= 0 && chargeMax >= 0 ){ // incude 0 charge
					string name = speciesName( pType, 0, ptBin, etaBin );	
					TH2* h2 = book->get2D( name );

					if ( h2 ){
						h2->Fill( dedx, tof );
					} else {
						// this one chould always exist
						lg->error(__FUNCTION__) << "Could not fill for " << name << endl;
					}
				}
				
				// fill the histogram for the charge if desired
				if ( chargeMin <= charge && chargeMax >= charge ){ // incude charge?
					string name = speciesName( pType, charge, ptBin, etaBin );	
					TH2* h2 = book->get2D( name );

					if ( h2 ){
						h2->Fill( dedx, tof );
					} else {
						// this one chould always exist
						lg->error(__FUNCTION__) << "Could not fill for " << name << endl;
					}
				}

			} // loop pTypes
    		
    	}
    	
	} // end loop on events


	// make particle type reports
	for ( int i = 0; i < plcs.size(); i++ ){
		distributionReport( plcs[ i ] );
	}

}

void pidHistogramMaker::prepareHistograms( string pType ) {
	
	lg->info(__FUNCTION__) << "Making Histograms with centering spceies: " << pType << endl;

	/**
	 * Make the dedx binning 
	 */
	double dedxBinWidth = config->getDouble( "binning.dedx:binWidth" );
	dedxMin = config->getDouble( "binning.dedx:min" );
	dedxMax = config->getDouble( "binning.dedx:max" );
	dedxBins = HistoBook::makeFixedWidthBins( dedxBinWidth, dedxMin, dedxMax );
	lg->info(__FUNCTION__) << "Dedx bins created ( " << dedxMin << ", " << dedxMax << " )" << endl;
	/**
	 * Make the Tof Binning
	 * could be for either 1/beta of delta 1/beta 
	 */
	tofMetric = config->getString( "binning.tofMetric", "inverseBeta" );
	double tofBinWidth = config->getDouble( "binning.tof:binWidth" );
	tofMin = config->getDouble( "binning.tof:min" );
	tofMax = config->getDouble( "binning.tof:max" );
	tofBins = HistoBook::makeFixedWidthBins( tofBinWidth, tofMin, tofMax );
	lg->info(__FUNCTION__) << "Tof bins created ( " << tofMin << ", " << tofMax << " )" << endl;

	/**
	 * Make the momentum transverse binning
	 * Can choose between fixed and variable width bins
	 */
	if ( config->nodeExists( "binning.ptBins" ) && config->getDoubleVector( "binning.ptBins" ).size() >= 2 ){
		ptBins = config->getDoubleVector( "binning.ptBins" );
	} else {
		// build the ptBins from the range and binWidth
		ptBins = HistoBook::makeFixedWidthBins( 
			config->getDouble( "binning.pt:binWidth", .05 ), 
			config->getDouble( "binning.pt:min", 0.2 ), 
			config->getDouble( "binning.pt:min", 0.2 ) 
		);
	}
	
	
	ptMin = ptBins[ 0 ];
	ptMax = ptBins[ ptBins.size() - 1 ];
	lg->info(__FUNCTION__) << "pT bins created ( " << ptMin << ", " << ptMax << " )" << endl;

	/**
	 * Make the eta binning
	 * Can choose between fixed and variable width bins
	 */
	if ( config->isVector( "binning.etaBins" ) ){
		etaBins = config->getDoubleVector( "binning.etaBins" );
	} else {
		etaBins = HistoBook::makeFixedWidthBins(
			config->getDouble( "binning.eta:binWidth", 0.2 ),
			config->getDouble( "binning.eta:min", 0.0 ),
			config->getDouble( "binning.eta:max", 1.0 )
		);
	}
	etaMin = etaBins[ 0 ];
	etaMax = etaBins[ etaBins.size() - 1 ];
	lg->info(__FUNCTION__) << "eta bins created ( " << etaMin << ", " << etaMax << " )" << endl;

	// this one could be length 0
	// in general not allowed for bins so use nodeExists instead
	if ( config->nodeExists( "binning.chargeBins" ) ){
		chargeBins = config->getDoubleVector( "binning.chargeBins" );
	} else {
		lg->info(__FUNCTION__) << "Vector?" <<  endl;
		chargeBins = HistoBook::makeFixedWidthBins(
			config->getDouble( "binning.charge:binWidth", 1 ),
			config->getDouble( "binning.charge:min", -1.0 ),
			config->getDouble( "binning.charge:max", 1.0 )
		);
	}
	chargeMin = chargeBins[ 0 ];
	chargeMax = chargeBins[ chargeBins.size() - 1 ];
	lg->info(__FUNCTION__) << "charge bins created ( " << chargeMin << ", " << chargeMax << " )" << endl;


	book->cd( "dedx_tof" );	

	// Loop through pt, then eta then charge
	for ( int ptBin = 0; ptBin < ptBins.size()-1; ptBin++ ){

		double p = ptBins[ ptBin ] + ptBins[ ptBin + 1 ];
		p /= 2.0;

		double tofLow, tofHigh, dedxLow, dedxHigh;
		autoViewport( pType, p, &tofLow, &tofHigh, &dedxLow, &dedxHigh, tofPadding, dedxPadding, tofScalePadding, dedxScalePadding );
		
		tofBins.clear();
		dedxBins.clear();
		tofBins = HistoBook::makeFixedWidthBins( tofBinWidth, tofLow, tofHigh );
		dedxBins = HistoBook::makeFixedWidthBins( dedxBinWidth, dedxLow, dedxHigh );

		for ( int etaBin = 0; etaBin < etaBins.size()-1; etaBin++ ){
			for ( int chargeBin = 0; chargeBin < chargeBins.size(); chargeBin++ ){

				// the name of the histogram
				string hName = speciesName( pType, chargeBins[ chargeBin ], ptBin, etaBin );

				string title = "dedx vs. tof; dedx; 1/#beta";

				// make it and keep it in the HistoBook
				book->make2D( hName, title, 
					dedxBins.size()-1, dedxBins.data(),
					tofBins.size()-1, tofBins.data() );

			}// loop on charge
		}// loop on eta bins
	} // loop on ptBins
	
	histosReady = true;

}


void pidHistogramMaker::speciesReport( string pType, int charge, int etaBin ){
	
	string name = speciesName( pType, charge );

	cout << "\tSpecies Report : " << name << endl;

	uint nBinsP = ptBins.size();
	

	TH3 * h3 = book->get3D( "nSig_" + name );

	taskProgress tp( pType + " report", nBinsP );

	for ( uint i = 0; i < nBinsP; i ++ ){

		tp.showProgress( i );

		// momentum value used for finding nice range
		double p = ptBins[ i ];


		pReport[ name ]->newPage( 2, 2 );
		pReport[ name ]->cd( 1, 2 );
		h3->GetZaxis()->SetRange( i, i );
		
		TH2* proj;
		TH1* proj1D;

		double pLow = h3->GetZaxis()->GetBinLowEdge( i );
		double pHi = h3->GetZaxis()->GetBinLowEdge( i + 1 );
		if ( 0 == i ){
			pLow = h3->GetZaxis()->GetBinLowEdge( 1 );
			pHi = h3->GetZaxis()->GetBinLowEdge( ptBins.size()-1 );
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

	uint nBinsPt = ptBins.size() - 1;
	string rName = speciesName( pType, 0 );

	taskProgress tp( pType + " distribution report", nBinsPt );

	book->cd( "tof" );
	for ( uint i = 0; i < nBinsPt; i ++ ){

		tp.showProgress( i );

		// momentum value used for finding nice range
		double p = ptBins[ i ];
		double p2 = ptBins[ i + 1 ];
		double avgP = 0.2;
		avgP = (ptBins[ i ] + ptBins[ i + 1])/2.0;
		

		string name = speciesName( pType, 0, i, 0 );

		book->cd( "dedx_tof" );
		TH2 * pTof = book->get2D( name );
		book->cd( "scratch" );
		TH2 * pDedx = (TH2*)pTof->Clone( "pDedx__" );

		// start a new page on the report file
		pReport[ rName ]->newPage( 2, 2 );


		// get information on plot ranges
		double tofLow, tofHigh, dedxLow, dedxHigh;
		autoViewport( pType, p, &tofLow, &tofHigh, &dedxLow, &dedxHigh, tofPadding, dedxPadding, tofScalePadding, dedxScalePadding );
		
		if ( true ) {	// show the tof proj

			string title = "#beta^{-1} : " + ts(ptBins[ i ], 4) + " < pT < " + ts(ptBins[i+1], 4);
			vector<string> others = otherSpecies( pType );
			vector< double > tofMean = enhanceTof( pType, others, avgP );
			vector< double > dedxMean = enhanceDedx( pType, others, avgP );
			pReport[ rName ]->cd( 1, 1 );
			//hdt->GetXaxis()->SetRangeUser( -.06, .06 );
			// Make the all tof tracks histogram
			string hName = sTofName( pType, 0, i );
			book->cd( "scratch" );
			TH1* hTof = (TH1D*)pTof->ProjectionY( "_py" );
			book->cd( "tof" );
			book->add( hName, (TH1*)hTof->Clone( hName.c_str() )  );
			book->style( hName )->set( "style.tof" )
			 	->set( "title", title )->draw();

			TLine * l1 = new TLine( tofMean[ 0 ], hTof->GetMinimum(), tofMean[ 0 ], hTof->GetMaximum() );
			l1->Draw();
			TLine * l2 = new TLine( tofMean[ 1 ], hTof->GetMinimum(), tofMean[ 1 ], hTof->GetMaximum() );
			l2->Draw();

			pReport[ rName ]->cd( 2, 1 );
			pTof->GetXaxis()->SetRangeUser( -.06, .06 );
			// Make the all tof tracks histogram
			hName = sTofName( pType, 0, i, 0, pType );
			book->cd( "scratch" );
			hTof = (TH1D*)pTof->ProjectionY( "_py" );
			book->cd( "tof" );
			book->add( hName, (TH1*)hTof->Clone( hName.c_str() )  );
			book->style( hName )->set( "style.tof" )
			 	->set( "title", title + " " + pType + " enhanced" )->draw();

			for ( int j = 0; j < dedxMean.size(); j++ ){

				pReport[ rName ]->cd( j+1, 2 );
				pTof->GetXaxis()->SetRangeUser( dedxMean[j]-0.06, dedxMean[j]+0.06 );
				// Make the all tof tracks histogram
				hName = sTofName( pType, 0, i, 0, others[ j ] );
				book->cd( "scratch" );
				hTof = (TH1D*)pTof->ProjectionY( "_py" );
				book->cd( "tof" ); 
				book->add( hName, (TH1*)hTof->Clone( hName.c_str() )  );
				book->style( hName )->set( "style.tof" )
				 	->set( "title", title + " " + others[ j ] + " enhanced" )->draw();
			}

		}

		pReport[ rName ]->savePage();
		pReport[ rName ]->newPage( 2, 2 );


		if ( true ) {	// show the dedx proj

			string title = "dEdx : " + ts(ptBins[ i ], 4) + " < pT < " + ts(ptBins[i+1], 4);
			pTof->GetXaxis()->SetRange( 1, pTof->GetXaxis()->GetNbins() );
			pTof->GetYaxis()->SetRange( 1, pTof->GetYaxis()->GetNbins() );

			vector<string> others = otherSpecies( pType );
			vector< double > tofMean = enhanceTof( pType, others, avgP );
			vector< double > dedxMean = enhanceDedx( pType, others, avgP );
			pReport[ rName ]->cd( 1, 1 );
			
			// Make the all dedx tracks histogram
			string hName = sDedxName( pType, 0, i );
			book->cd( "scratch" );
			TH1* hDedx = (TH1D*)pTof->ProjectionX( "_px" );
			book->cd( "dedx" );
			book->add( hName, (TH1*)hDedx->Clone( hName.c_str() )  );
			book->style( hName )->set( "style.dedx" )
			 	->set( "title", title )->draw();

			TLine * l1 = new TLine( dedxMean[ 0 ], hDedx->GetMinimum(), dedxMean[ 0 ], hDedx->GetMaximum() );
			l1->Draw();
			TLine * l2 = new TLine( dedxMean[ 1 ], hDedx->GetMinimum(), dedxMean[ 1 ], hDedx->GetMaximum() );
			l2->Draw();

			pReport[ rName ]->cd( 2, 1 );
			pTof->GetYaxis()->SetRangeUser( -.012, .012 );
			// Make the all tof tracks histogram
			hName = sDedxName( pType, 0, i, 0, pType );
			book->cd( "scratch" );
			hDedx = (TH1D*)pTof->ProjectionX( "_px" );
			book->cd( "dedx" );
			book->add( hName, (TH1*)hDedx->Clone( hName.c_str() )  );
			book->style( hName )->set( "style.dedx" )
			 	->set( "title", title + " " + pType + " enhanced" )->draw();
			
			for ( int j = 0; j < dedxMean.size(); j++ ){

				pReport[ rName ]->cd( j+1, 2 );
				pTof->GetYaxis()->SetRangeUser( tofMean[j]-0.012, tofMean[j]+0.012 );
				// Make the all tof tracks histogram
				hName = sDedxName( pType, 0, i, 0, others[ j ] );
				book->cd( "scratch" );
				hDedx = (TH1D*)pTof->ProjectionX( "_px" );
				book->cd( "dedx" ); 
				book->add( hName, (TH1*)hDedx->Clone( hName.c_str() )  );
				book->style( hName )->set( "style.dedx" )
				 	->set( "title", title + " " + others[ j ] + " enhanced" )->draw();
			}

		}

		pReport[ rName ]->savePage();
		
	}



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
	
	


	*tofHigh = ( tHigh  / tofPlotSigma) + tofPadding;
	*tofLow = ( tLow / tofPlotSigma) - tofPadding;

	*dedxHigh = ( dHigh / dedxPlotSigma) + dedxPadding;
	*dedxLow = ( dLow / dedxPlotSigma) - dedxPadding;

	double tofRange = *tofHigh - *tofLow;
	double dedxRange = *dedxHigh - *dedxLow;

	double scaledPaddingTof = tofRange * tofScaledPadding;
	double scaledPaddingDedx = dedxRange * dedxScaledPadding;

	*tofLow -= scaledPaddingTof;
	*tofHigh += scaledPaddingTof;

	*dedxLow -= scaledPaddingDedx;
	*dedxHigh += scaledPaddingDedx;

/*
	if ( *tofLow < tofMin )
		*tofLow = tofMin;
	if ( *tofHigh > tofMax )
		*tofHigh = tofMax;

	if ( *dedxLow < dedxMin )
		*dedxLow = dedxMin;
	if ( *dedxHigh > dedxMax )
		*dedxHigh = dedxMax;*/

	return;

}


/**
 *  N Sigma Calculations ************************************
 */

/**
 * Calculates the difference from the expected value of log10( dedx )
 * @param  pType particle type
 * @param  iHit  hit index in the nTuple
 * @return       n Sigma from expectation
 */
double pidHistogramMaker::nSigmaDedx( string pType, int iHit ){ 
			
	double p = pico->p[ iHit ];

	double mean = TMath::Log10( dedxGen->mean( p, eMass( pType ) ) * 1000 );
	double dedx = TMath::Log10( pico->dedx[ iHit ]);
	double nSig = (( dedx - mean ) / dedxPlotSigma );

	return nSig;

	return -999.0;
}

/**
 * Calculates the difference from the expected value of log10( dedx )
 * Uses the non-linear recentering scheme
 * @param  pType particle type
 * @param  iHit  hit index in the nTuple
 * @return       n Sigma from expectation
 */
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
		double w = dedx + iMuAvg - iMu;
		

		n1 += (iL * w);
		d1 += iL;

		double iLc = lh( mu, iMu, sigma );
		n2 += (iLc * w);
		d2 += iLc;
	}

	double nSig = (n1/d1) - (muAvg );

	return (nSig / dedxPlotSigma);

}


/**
 * Calculates the difference from the expected value of 1/beta
 * @param  pType particle type
 * @param  iHit  hit index in the nTuple
 * @return       n Sigma from expectation
 */
double pidHistogramMaker::nSigmaInverseBeta( string pType, int iHit  ){

	double betaMeasured = pico->beta[ iHit ];
	double p = pico->p[ iHit ];
	double betaExpected = eBeta( eMass( pType ), p );
	double deltaInvBeta = ( 1.0 / betaMeasured ) - ( 1.0 / betaExpected );

	return (deltaInvBeta / tofPlotSigma);
}

/**
 * Calculates the difference from the expected value of 1/beta
 * Uses the nonlinear recentering scheme
 * @param  pType particle type
 * @param  iHit  hit index in the nTuple
 * @return       n Sigma from expectation
 */
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
		
		double sigma = tofSigma; 

		double iL = lh( tof, iMu, sigma );
		
		double w = tof + iMuAvg - iMu;
		
		n1 += (iL * w);
		d1 += iL;

		double iLc = lh( mu, iMu, sigma );
		n2 += (iLc * w);
		d2 += iLc;
	}

	
	double nSig = (n1/d1) - ( muAvg );

	return (nSig /  tofPlotSigma);

}

/**
 * Likelihood function
 * A gauss around the expected value with expected sigma
 * @param  x     measured value
 * @param  mu    expected mean
 * @param  sigma expected sigma
 * @return       the unnormalized likelihood
 */
double pidHistogramMaker::lh( double x, double mu, double sigma ){

	double a = sigma * TMath::Sqrt( 2 * TMath::Pi() );
	double b = ( x - mu );
	double c = 2 * sigma*sigma;
	double d = (1/a) * TMath::Exp( -b*b / c );

	return d;
}






