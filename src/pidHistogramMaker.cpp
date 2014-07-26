
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
	
	gErrorIgnoreLevel=kError;

	config = con;


	// set the histogram info verbosity to show nothing
	gStyle->SetOptStat( 0 );
	
	// create the histogram book
	book = new histoBook( ( config->getString( "output.base" ) + config->getString( "output.root" ) ), config, config->getString( "input.root" ) );
	


	vector<string> parts = config->getStringVector( "pType" );
	for ( int i = 0; i < parts.size(); i++ ){
		for ( int charge = -1; charge <= 1; charge ++ ){
			string n = sName( parts[ i ], charge );
				pReport[ n ] = new reporter( config->getString( "output.base" ) + n + config->getString( "output.report" ) );		
		}
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
		delete pReport[ sName( parts[ i ], -1 ) ];
		delete pReport[ sName( parts[ i ], 0 ) ];
		delete pReport[ sName( parts[ i ], 1 ) ];
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
	/*
	TH2D* h = new TH2D( "dan", "dan", 50, -10, 40, 50, -10, 40 );
	
   	for (Int_t i = 0; i < 500000; i++) {
    	double px = gRandom->Gaus( 25, 3 );
    	double py = gRandom->Gaus( 25, 3 );
     	h->Fill( px, py );
   	}
   	for (Int_t i = 0; i < 500000; i++) {
    	double px = gRandom->Gaus( 15, 3 );
    	double py = gRandom->Gaus( 15, 3 );
     	h->Fill( px, py );
   	}
   	for (Int_t i = 0; i < 500000; i++) {
    	double px = gRandom->Gaus( 5, 3 );
    	double py = gRandom->Gaus( 5, 3 );
     	h->Fill( px, py );
   	}


	dklMinimizer * dkl = new dklMinimizer( h, 3 );

	//dkl->printAll();

	dkl->run( 50000 );

	//dkl->printAll();
	cout << "iy : " << dkl->inputYield() << endl;
	cout << "ay : " << dkl->approximationYield() << endl;	
	cout << "s0y : " << dkl->speciesYield( 0 ) << endl;	
	cout << "s1y : " << dkl->speciesYield( 1 ) << endl;	
	cout << "s2y : " << dkl->speciesYield( 2 ) << endl;	

	report->newPage(1, 2);
	gPad->SetLogz(1);
	dkl->viewInput( ) ->Draw("colz");
	report->cd( 1, 2 );
	gPad->SetLogz(1);
	dkl->viewApproximation( )->Draw("colz");
	report->savePage();

	report->newPage(1, 2);
	gPad->SetLogz(1);
	dkl->viewSpecies( 0 ) ->Draw("colz");
	report->cd( 1, 2 );
	gPad->SetLogz(1);
	dkl->viewSpecies( 1 ) ->Draw("colz");
	report->savePage();

	report->newPage(1, 2);
	gPad->SetLogz(1);
	dkl->viewSpecies( 2 ) ->Draw("colz");

	report->savePage();
	
	return;
	*/



	if ( !_chain ){
		cout << "[pidHistogramMaker." << __FUNCTION__ << "] ERROR: Invalid chain " << endl;
		return;
	}

	Int_t nevents = (Int_t)_chain->GetEntries();
	cout << "[pidHistogramMaker." << __FUNCTION__ << "] Loaded: " << nevents << " events " << endl;

	book->cd( "" );

	vector<string> parts = config->getStringVector( "pType" );

	if ( !book->get( "nSig_" + sName( parts[ 0 ], 0 ) ) ){
		for ( int i = 0; i < parts.size(); i++ ){
			sHisto( parts[ i ] );
		}

		book->make( "nSigBetaDedxK" );

	
	
		// loop over all events
		for(Int_t i=0; i<nevents; i++) {
	    	_chain->GetEntry(i);

	    	progressBar( i, nevents, 60 );

	    	double vz = pico->vertexZ;
	    	int nTof = pico->nTofHits;
	    	int nT0 = pico->nTZero;
	    	if ( TMath::Abs( vz ) > config->getDouble( "cut.vZ", 30 ) )
	    		continue;
	    	if ( nTof < 10 || nT0 < 10 )
	    		continue;
	    	

	    	int nTofHits = pico->nTofHits;
	    	for ( int iHit = 0; iHit < nTofHits; iHit++ ){
	  		
	    		double p = pico->p[ iHit ];

	    		for ( int i = 0; i < parts.size(); i++ ){
	    			string pType = parts[ i ];
	    			
					int charge = pico->charge[ iHit ];
					double dedx = nSigDedx( pType, iHit );
					double invBeta = nSigInvBeta( pType, iHit);
					
					
					if ( p > pMax || p < pMin )
						continue;
					if ( dedx > nSigMax || dedx < nSigMin )
						continue;
					if ( invBeta > nSigMax || invBeta < nSigMin )
						continue;

					double angle = 0 * ( 3.1415926 / 180 );
					
					//double rX = dedx;
					//double rY = invBeta;
					double rX = dedx * TMath::Cos( angle ) - invBeta * TMath::Sin( angle );
					double rY = dedx * TMath::Sin( angle ) + invBeta * TMath::Cos( angle );

					double pR = TMath::Hypot( dedx, invBeta );
					double pT = TMath::ATan2( invBeta, dedx );

					string name = "nSig_" + sName( pType, charge );
					TH3* h3 = ((TH3*)book->get( name ));
					if ( h3 )
						h3->Fill( dedx, invBeta, p );

					// always fill the both charges histo
					name = "nSig_" + sName( pType, 0 );
					h3 = ((TH3*)book->get( name ));
					if ( h3 )
						h3->Fill( rX, rY, p );

					if ( "K" == pType ){
						book->fill( "nSigBetaDedxK", dedx, invBeta );

					}
				}
	    		
	    	}
	    	
		} // end loop on events
	}

	// make particle type reports
	for ( int i = 0; i < parts.size(); i++ ){
		//speciesReport( parts[ i ], -1 );
		speciesReport( parts[ i ], 0 );
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

	//int nPBins = config->getDouble( "binning.p:nBins" );
	
	pMin = config->getDouble( "binning.p:min" );
	pMax = config->getDouble( "binning.p:max" );





	// make the nsigma bins
	double step = (nSigMax - nSigMin ) / nSigBins;
	vector<double>nSigBinEdges;
	for ( double i = nSigMin; i < nSigMax; i += step ){
		nSigBinEdges.push_back( i );
	}
	nSigBinEdges.push_back( nSigMax );


	//vector<double> ppBins = { 0, .5, 1.0, 1.5, 2.0, 2.25, 3, 4 };

	vector<double> pBins = config->getDoubleVector( "binning.pBins" );

	string title = "; n#sigma dedx; n#sigma 1/#beta ";
	// create a combined, plus, and minus
	for ( int charge = -1; charge <= 1; charge ++ ){
		string name = "nSig_" + sName( pType, charge );

		TH3D * h3 = new TH3D( name.c_str(), title.c_str(), 
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

void pidHistogramMaker::speciesReport( string pType, int charge ){

	string name = sName( pType, charge );
	//TCanvas * can = new TCanvas( "can", "can", 800, 800 );
	//can->Print( "test.pdf[");
	vector<double>pBins = config->getDoubleVector( "binning.pBins" );
	bool fitGauss = config->getBool( "pReport.fit1DGauss", false );
	double fitX1 = config->getDouble( "pReport.fit1DGauss:x1", nSigMin );
	double fitX2 = config->getDouble( "pReport.fit1DGauss:x2", nSigMax );
	
	TH3 * h3 = book->get3D( "nSig_" + name );
	for ( int i = 0; i < pBins.size(); i ++ ){

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
		proj->SetTitle( (pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 ) ).c_str()  );
		gPad->SetLogz( 1 );
		proj->Draw( "colz" );


		
		if ( fitGauss ){
			pReport[ name ]->cd( 2, 1 );
			proj1D = (TH1D*)h3->Project3D( "x" )->Clone( "fit");
			proj1D->SetTitle( ("dedx : " + pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 ) ).c_str()  );
			gPad->SetLogy( 1 );
			proj1D->Draw( "h" );
			TF1 * f1 = new TF1( "g1", "gaus", fitX1, fitX2 );
			f1->SetRange( fitX1, fitX2 );
			f1->SetParameter( 1, 0 );
			proj1D->Fit( f1, "QR", "", fitX1, fitX2 );
		}

		
		pReport[ name ]->cd( 2, 2 );
		proj1D = h3->Project3D( "x" );
		proj1D->SetTitle( ( "dedx : " + pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 ) ).c_str()  );
		proj1D->SetFillColor( kBlue );
		gPad->SetLogx( 1 );
		proj1D->Draw( "hbar" );

		//can->Print( "test.pdf" );
		

		pReport[ name ]->cd( 1, 1 );
		proj1D = h3->Project3D( "y" );
		proj1D->SetTitle( ( "1/#beta : " + pType + " : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 ) ).c_str()  );
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

		if ( pType == config->getString( "dklFit.pType" ) && charge == 0 && i >= config->getInt( "dklFit.pBin:min", 1 ) && i <= config->getInt( "dklFit.pBin:max", 1 )){
			//dklFit( name, (TH2D*)proj );
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









