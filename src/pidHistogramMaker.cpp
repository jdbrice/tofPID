
#include "constants.h"
#include "pidHistogramMaker.h"
#include "histoBook.h"
#include "dklMinimizer.h"
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
	
	vector<double> etaBins = config->getDoubleVector( "binning.eta" );
	int nEta = etaBins.size();

	if ( true == config->nodeExists( "pType" )  ){
		vector<string> parts = config->getStringVector( "pType" );
		
		for ( int i = 0; i < parts.size(); i++ ){
			for ( int charge = -1; charge <= 1; charge ++ ){

				string n = sName( parts[ i ], charge );
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
			delete pReport[ sName( parts[ i ], -1 ) ];
			delete pReport[ sName( parts[ i ], 0 ) ];
			delete pReport[ sName( parts[ i ], 1 ) ];
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

	report->newPage();
	book->style("iBeta")->draw();

	TGraph * g1 = inverseBeta( constants::piMass, 0.15, 3, .05 );
	TGraph * g2 = inverseBeta( constants::kaonMass, 0.15, 3, .05 );
	TGraph * g3 = inverseBeta( constants::protonMass, 0.15, 3, .05 );
	TGraph * g4 = inverseBeta( constants::eMass, 0.15, 3, .05 );
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

	vector<double> etaBins = config->getDoubleVector( "binning.eta" );
	int nEta = etaBins.size();

	for ( int i = 0; i < parts.size(); i++ ){
		sHisto( parts[ i ] );
	}

	// loop over all events
	for(Int_t i=0; i<nevents; i++) {
    	_chain->GetEntry(i);

    	progressBar( i, nevents, 60 );

    	if ( !keepEventQA() )
    		continue;
    	

    	int nTofHits = pico->nTofHits;
    	for ( int iHit = 0; iHit < nTofHits; iHit++ ){

    		if ( !keepTrackQA( iHit ) )
	    		continue;
  		
    		double p = pico->p[ iHit ];

    		for ( int i = 0; i < parts.size(); i++ ){
    			string pType = parts[ i ];
    			
				int charge = pico->charge[ iHit ];
				double eta = pico->eta[ iHit ];
				double dedx = nSigDedx( pType, iHit );
				double invBeta = nSigInvBeta( pType, iHit);
				double deltaB = dBeta( pType, iHit );
				
				
				if ( p > pMax || p < pMin )
					continue;
				if ( dedx >= nSigMax || dedx <= nSigMin )
					continue;
				if ( deltaB >= dBetaMax || deltaB <= dBetaMin )
					continue;

				//double angle = 0 * ( 3.1415926 / 180 );
				
				//double rX = dedx;
				//double rY = invBeta;
				//double rX = dedx * TMath::Cos( angle ) - invBeta * TMath::Sin( angle );
				//double rY = dedx * TMath::Sin( angle ) + invBeta * TMath::Cos( angle );

				//double pR = TMath::Hypot( dedx+5, invBeta+5 );
				//double pT = TMath::ATan2( invBeta+5, dedx+5 );
				

				string name = "nSig_" + sName( pType, charge );
				
				TH3* h3 = ((TH3*)book->get( name ));
				
				if ( h3 ){
					h3->Fill( dedx, deltaB, p );
				}

				// always fill the both charges histo
				name = "nSig_" + sName( pType, 0 );
				h3 = ((TH3*)book->get( name ));
				if ( h3 ){
					h3->Fill( dedx, deltaB, p );
				}

				for ( int iEta = 0; iEta < nEta; iEta++ ){
					double etaHigh = etaBins[ iEta ];
					double etaLow = 0.0;


					if ( iEta >= 1 )
						etaLow = etaBins[ iEta - 1];

					if ( TMath::Abs( eta ) >= etaLow && TMath::Abs( eta ) < etaHigh ){
						
						name = "nSig_" + sName( pType, 0 ) + "_eta" + ts( iEta);
						
						h3 = book->get3D( name );
						if ( h3 ){
							h3->Fill( dedx, deltaB, p );
						}

					}

				} // loop over the eta bins

			} // loop pTypes
    		
    	}
    	
	} // end loop on events

	// make particle type reports
	for ( int i = 0; i < parts.size(); i++ ){
		//speciesReport( parts[ i ], -1 );
		speciesReport( parts[ i ], 0 );
		//speciesReport( parts[ i ], 0, 0 );
		//speciesReport( parts[ i ], 1 );
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
	nSigMin = config->getDouble( "binning.nSig:min" );
	nSigMax = config->getDouble( "binning.nSig:max" );

	int dBetaBins = config->getDouble( "binning.deltaBeta:nBins" );
	dBetaMin = config->getDouble( "binning.deltaBeta:min" );
	dBetaMax = config->getDouble( "binning.deltaBeta:max" );

	vector<double> etaBins = config->getDoubleVector( "binning.eta" );
	
	
	pMin = config->getDouble( "binning.p:min" );
	pMax = config->getDouble( "binning.p:max" );

	// make the nsigma bins
	double step = (nSigMax - nSigMin ) / nSigBins;
	vector<double>nSigBinEdges;
	for ( double i = nSigMin; i < nSigMax; i += step ){
		nSigBinEdges.push_back( i );
	}
	nSigBinEdges.push_back( nSigMax );

	// make the delta beta bins
	double stepB = (dBetaMax - dBetaMin ) / dBetaBins;
	vector<double>dBetaBinEdges;
	for ( double i = dBetaMin; i < dBetaMax; i += stepB ){
		dBetaBinEdges.push_back( i );
	}
	dBetaBinEdges.push_back( dBetaMax );



	vector<double> pBins = config->getDoubleVector( "binning.pBins" );

	string title = "; n#sigma dedx; #Delta #beta^{-1} / #beta^{-1} ";
	
	
	// create a combined, plus, and minus
	for ( int charge = -1; charge <= 1; charge ++ ){
		for ( int iEta = 0; iEta < etaBins.size(); iEta++ ){

			double eta = etaBins[ iEta ];

			string name = "nSig_" + sName( pType, charge ) + "_eta" + ts( iEta );

			TH3D * h3 = new TH3D( name.c_str(), title.c_str(), 
					nSigBinEdges.size()-1, nSigBinEdges.data(), 
					dBetaBinEdges.size()-1, dBetaBinEdges.data(),
					pBins.size()-1, pBins.data() );

			book->add( name, h3 );
		}


		string name = "nSig_" + sName( pType, charge ) ;

		TH3D * h3 = new TH3D( name.c_str(), title.c_str(), 
				nSigBinEdges.size()-1, nSigBinEdges.data(), 
				dBetaBinEdges.size()-1, dBetaBinEdges.data(),
				pBins.size()-1, pBins.data() );

		book->add( name, h3 );

	}
	


}
double pidHistogramMaker::dBeta( string pType, int iHit ){

	double tof = pico->tof[ iHit ];
	double length = pico->length[ iHit ];
	double p = pico->p[ iHit ];
	double beta = pico->beta[ iHit ];
	double m2 = p*p * ( constants::c*constants::c * tof*tof / ( length*length ) - 1  );


	double deltaB = 1 - (beta) * TMath::Sqrt( (constants::kaonMass*constants::kaonMass) / (p*p) + 1 );

	if ( "Pi" == pType )
		deltaB = 1 - (beta) * TMath::Sqrt( (constants::piMass*constants::piMass) / (p*p) + 1 );		
	if ( "P" == pType )
		deltaB = 1 - (beta) * TMath::Sqrt( (constants::protonMass*constants::protonMass) / (p*p) + 1 );		
	
	return deltaB;
}

double pidHistogramMaker::nSigInvBeta( string pType, int iHit  ){

	double b = pico->beta[ iHit ];
	double p = pico->p[ iHit ];
	double expectedBeta = eBeta( eMass( pType ), p );
	double invBetaSig = config->getDouble( "binning.invBetaSig" );

	double deltaInvBeta = ( 1.0 / b ) - ( 1.0 / expectedBeta );

	return (deltaInvBeta / invBetaSig);
}

void pidHistogramMaker::speciesReport( string pType, int charge, int etaBin ){
	
	
	vector<double> etaBins = config->getDoubleVector( "binning.eta" );
	string name = sName( pType, charge );

	cout << "\tSpecies Report : " << name << endl;

	string etaRange = "";
	if ( 0 == etaBin)
		etaRange = " && 0 < |#eta| < " + ts( etaBins[ etaBin ], 3 );
	else if ( 0 < etaBin )
		if ( 0 == etaBin)
		etaRange = " && " + ts( etaBins[ etaBin - 1 ], 3 ) + " < |#eta| < " + ts( etaBins[ etaBin ], 3 );
	
	vector<double>pBins = config->getDoubleVector( "binning.pBins" );
	bool fitGauss = config->getBool( "pReport.fit1DGauss", false );
	double fitX1 = config->getDouble( "pReport.fit1DGauss:x1", nSigMin );
	double fitX2 = config->getDouble( "pReport.fit1DGauss:x2", nSigMax );
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
		string hTitle = (pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 ) + etaRange );
		string hName = (pType + "_" + ts(i) + "_" + ts(etaBin) );
		proj->SetTitle( hTitle.c_str()  );
		gPad->SetLogz( 1 );
		proj->Draw( "colz" );
		

		
		if ( fitGauss ){
			pReport[ name ]->cd( 2, 1 );
			proj1D = (TH1D*)h3->Project3D( "x" )->Clone( "fit");
			proj1D->SetTitle( ("dedx : " + pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 ) + etaRange ).c_str()  );
			gPad->SetLogy( 1 );
			proj1D->Draw( "h" );
			TF1 * f1 = new TF1( "g1", "gaus", fitX1, fitX2 );
			f1->SetRange( fitX1, fitX2 );
			f1->SetParameter( 1, 0 );
			proj1D->Fit( f1, "QR", "", fitX1, fitX2 );
		}

		
		pReport[ name ]->cd( 2, 2 );
		proj1D = h3->Project3D( "x" );
		proj1D->SetTitle( ( "dedx : " + pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 ) + etaRange ).c_str()  );
		proj1D->SetFillColor( kBlue );
		gPad->SetLogx( 1 );
		proj1D->Draw( "hbar" );

		//can->Print( "test.pdf" );
		

		pReport[ name ]->cd( 1, 1 );
		proj1D = h3->Project3D( "y" );
		proj1D->SetTitle( ( "1/#beta : " + pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 ) + etaRange ).c_str()  );
		gPad->SetLogy( 1 );
		proj1D->SetFillColor( kBlue );
		proj1D->Draw( "" );

		if ( fitGauss ){
			TF1 * f2 = new TF1( "g2", "gaus", fitX1, fitX2);
			f2->SetRange( fitX1, fitX2 );
			//f2->SetParameter( 1, 0.0 );
			proj1D->Fit( f2, "QR", "", fitX1, fitX2 );
		}
		

		pReport[ name ]->savePage();
		//pReport[ name ]->saveImage( "report/png/" + hName + ".png" );

		if ( pType == config->getString( "dklFit.pType" ) && charge == 0 && i >= config->getInt( "dklFit.pBin:min", 1 ) && i <= config->getInt( "dklFit.pBin:max", 1 )){
			dklFit( name, (TH2D*)proj );
		}

	}

	//can->Print( "test.pdf]");

}


void pidHistogramMaker::dklFit( string pName, TH2D * h ) {

	dklMinimizer *dkl = new dklMinimizer( h, config->getInt( "dklFit.nSpecies", 1 ) );

	dkl->run( config->getInt( "dklFit.nIterations", 1 ) );

	pReport[ pName ]->newPage( 3, 2 );

	pReport[ pName ]->cd( 1, 1);
	gPad->SetLogz(1);
	dkl->viewInput()->Draw("colz");

	pReport[ pName ]->cd( 2, 1);
	gPad->SetLogz(1);
	TH2D* ap = dkl->viewApproximation();
	double min = ap->GetMinimum();
	double max = ap->GetMaximum();
	ap->Draw("colz");

	pReport[ pName ]->cd( 1, 2);
	gPad->SetLogz(1);

	TH2D* s1 = dkl->viewSpecies( 0 );
	s1->GetZaxis()->SetRangeUser( min, max );
	s1->Draw("colz");

	if ( config->getInt( "dklFit.nSpecies" ) >= 2 ){
		pReport[ pName ]->cd( 2, 2);
		gPad->SetLogz(1);
		TH2D* s2 = dkl->viewSpecies( 1 );
		s2->GetZaxis()->SetRangeUser( min, max );
		s2->Draw("colz");
	}

	if ( config->getInt( "dklFit.nSpecies" ) >= 3 ){
		pReport[ pName ]->cd( 3, 2);
		gPad->SetLogz(1);
		TH2D* s3 = dkl->viewSpecies( 2 );
		s3->GetZaxis()->SetRangeUser( min, max );
		s3->Draw("colz");
	}

	pReport[ pName ]->savePage();



}









