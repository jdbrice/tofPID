

#include "pidFitter.h"
#include "dklMinimizer.h"
#define sDelete(x) {delete x; x = NULL;}

// provides my own string shortcuts etc.
using namespace jdbUtils;

pidFitter::pidFitter( xmlConfig * con ){
	cout << "[pidFitter.pidFitter] " << endl;
	
	gErrorIgnoreLevel=kError;

	config = con;

	// set the histogram info verbosity to show nothing
	gStyle->SetOptStat( 0 );
	
	// create the histogram book
	book = new histoBook( ( config->getString( "output.base" ) + config->getString( "output.root" ) ), config, config->getString( "input.root:file" ) );	


	vector<string> parts = config->getStringVector( "pType" );
	for ( int i = 0; i < parts.size(); i++ ){
		for ( int charge = -1; charge <= 1; charge ++ ){
			string n = sName( parts[ i ], charge );
				pReport[ n ] = new reporter( config->getString( "output.base" ) + n + ".pdf" );		
		}
	}

	// create a report builder 
	report = new reporter( config->getString( "output.base" ) + config->getString( "output.report" ) );

}


pidFitter::~pidFitter() {
	cout << "[pidFitter." << __FUNCTION__ << "]" << endl;
	delete book;
	delete report;

	vector<string> parts = config->getStringVector( "pType" );
	for ( int i = 0; i < parts.size(); i++ ){
		delete pReport[ sName( parts[ i ], -1 ) ];
		delete pReport[ sName( parts[ i ], 0 ) ];
		delete pReport[ sName( parts[ i ], 1 ) ];
	}	
}

string pidFitter::sName( string pType, int charge ){

	if ( -1 == charge )
		return pType + "_Negative";
	if ( 1 == charge )
		return pType + "_Positive";
	if ( 0 == charge )
		return pType + "_All";
	return "";
}


void pidFitter::runDkl( TH2D* h, reporter * rp, uint nS, uint nIt ){

	if ( nS <= 0 )
		return;
	if ( nIt <= 0 )
		return;
	if ( NULL == h || NULL == rp )
		return;

	double axisMin = config->getDouble( "binning.nSig:min", 0 );
	double axisMax = config->getDouble( "binning.nSig:max", 0 );

	dklMinimizer *dkl = new dklMinimizer( h, nS );
	dkl->run( nIt );

	rp->newPage( 1, 2 );

	rp->cd( 1, 1);
	gPad->SetLogz(1);
	dkl->viewInput( axisMin, axisMax, axisMin, axisMax)->Draw("colz");

	rp->cd( 2, 1);
	gPad->SetLogz(1);
	TH2D* ap = dkl->viewApproximation( axisMin, axisMax, axisMin, axisMax );
	double min = ap->GetMinimum();
	double max = ap->GetMaximum();
	ap->Draw("colz");
	rp->savePage();

	rp->newPage( 1, 3 );
	rp->cd( 1, 1);
	gPad->SetLogz(1);

	TH2D* s1, *s2, *s3;
	s1 = dkl->viewSpecies( 0, axisMin, axisMax, axisMin, axisMax );
	s1->GetZaxis()->SetRangeUser( min, max );
	s1->Draw("colz");

	if ( nS >= 2 ){
		rp->cd( 1, 2);
		gPad->SetLogz(1);
		s2 = dkl->viewSpecies( 1, axisMin, axisMax, axisMin, axisMax );
		s2->GetZaxis()->SetRangeUser( min, max );
		s2->Draw("colz");
	}

	if ( nS >= 3 ){
		rp->cd( 1, 3);
		gPad->SetLogz(1);
		s3 = dkl->viewSpecies( 2, axisMin, axisMax, axisMin, axisMax );
		s3->GetZaxis()->SetRangeUser( min, max );
		s3->Draw("colz");
	}

	rp->savePage();

	sDelete( s1 );
	sDelete( s2 );
	sDelete( s3 );
}


void pidFitter::runFit(){

	
	processSpecies( "K", 0, report );

}


void pidFitter::processSpecies( string species, int charge, reporter * rp ){
	cout << "[pidFitter." << __FUNCTION__ << "]" << endl;

	string hName = "nSig_" + sName( species, charge );
	// get the pt Binning
	vector<double>pBins = config->getDoubleVector( "binning.pBins" );

	TH3* h3 = book->get3D( hName );
	for ( int i = 1; i < pBins.size(); i ++ ){

		// get the Pt range for title etc.
		double pLow = h3->GetZaxis()->GetBinLowEdge( i );
		double pHi = h3->GetZaxis()->GetBinLowEdge( i + 1 );

		// get the config entry for this p bin
		string nodePath = species + "_Fit.p" + ts(i);
		if ( !config->nodeExists( nodePath ) )
			continue;

		// look at one Pt bin at a time
		h3->GetZaxis()->SetRange( i, i );
		
		// Get the 2D projection we want
		TH2D* proj;
		proj = (TH2D*)h3->Project3D( "xy" );
		string name = proj->GetName();
		TH2D* cut = (TH2D*) proj->Clone( (name + "cut").c_str() );
		rp->newPage( 1, 2 );
		proj->SetTitle( ( ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 ) ).c_str()  );
		gPad->SetLogz( 1 );
		proj->Draw( "colz" );

		// process the square cuts
		makeSquareCuts( cut, nodePath + ".squareCut" );

		// draw the distribution after square cuts
		rp->cd( 1, 2 );
		gPad->SetLogz( 1 );
		cut->SetTitle( ( "After Square Cuts : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 ) ).c_str()  );
		cut->Draw( "colz" );
		rp->savePage( );

		// now fit using the dkl algorithm
		runDkl( cut, rp, config->getInt( nodePath + ".dkl:nSpecies", 1 ), config->getInt( nodePath + ".dkl:nRuns", 10 ) );


		
	}
}


double pidFitter::squareCut( TH2D* h, string axis, string cut, double value ){

	if ( "x" != axis && "y" != axis && "X" != axis && "Y" != axis )
		return 0;
	if ( ">" != cut && "<" != cut )
		return 0;
	if ( NULL == h)
		return 0;

	int nX = h->GetNbinsX();
	int nY = h->GetNbinsY();

	double nRemoved = 0;

	for ( int bX = 1; bX <= nX; bX++ ){
		for ( int bY = 1; bY <= nY; bY++ ){

			double xEdge = h->GetXaxis()->GetBinLowEdge( bX );
			double yEdge = h->GetXaxis()->GetBinLowEdge( bY );

			if ( "x" == axis ){
				if ( xEdge > value && ">" == cut ){
					nRemoved += h->GetBinContent( bX, bY );
					h->SetBinContent( bX, bY, 0);
					h->SetBinError( bX, bY, 0);
				} else if ( xEdge < value && "<" == cut ){
					nRemoved += h->GetBinContent( bX, bY );
					h->SetBinContent( bX, bY, 0);
					h->SetBinError( bX, bY, 0);
				}
			}

			if ( "y" == axis ){
				if ( yEdge > value && ">" == cut ){
					nRemoved += h->GetBinContent( bX, bY );
					h->SetBinContent( bX, bY, 0);
					h->SetBinError( bX, bY, 0);
				} else if ( yEdge < value && "<" == cut ){
					nRemoved += h->GetBinContent( bX, bY );
					h->SetBinContent( bX, bY, 0);
					h->SetBinError( bX, bY, 0);
				}
			}

		} // loop bY
	} // loop bX

	//cout << "Square Cut : " << h->GetName() << " Removed : "<< nRemoved << endl;
}

void pidFitter::makeSquareCuts( TH2D* h, string nodePath ){

	// four possible square cuts
	for ( int i = 1; i <= 4; i++ ){

		string axis = config->getString( nodePath + ":axis"+ ts(i), "" );
		string cut = config->getString( nodePath + ":cut"	+ ts(i), "" );
		double val = config->getDouble( nodePath + ":value"+ ts(i), 0 );

		squareCut( h, axis, cut, val );

	}
}




