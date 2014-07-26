

#include "pidFitter.h"


// provides my own string shortcuts etc.
using namespace jdbUtils;

pidFitter::pidFitter( xmlConfig * con ){
	cout << "[pidFitter.pidFitter] " << endl;
	
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
				pReport[ n ] = new reporter( config->getString( "output.base" ) + n + ".pdf" );		
		}
	}

	// create a report builder 
	report = new reporter( config->getString( "output.base" ) + config->getString( "output.report" ) );

}


pidFitter::~pidFitter() {
	
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

	dklMinimizer *dkl = new dklMinimizer( h, nS );
	dkl->run( nIt );

	rp->newPage( 1, 2 );

	rp->cd( 1, 1);
	gPad->SetLogz(1);
	dkl->viewInput()->Draw("colz");

	rp->cd( 2, 1);
	gPad->SetLogz(1);
	TH2D* ap = dkl->viewApproximation();
	double min = ap->GetMinimum();
	double max = ap->GetMaximum();
	ap->Draw("colz");
	rp->savePage();

	rp->newPage( 1, 3 );
	rp->cd( 1, 1);
	gPad->SetLogz(1);

	TH2D* s1 = dkl->viewSpecies( 0 );
	s1->GetZaxis()->SetRangeUser( min, max );
	s1->Draw("colz");

	if ( config->getInt( "dklFit.nSpecies" ) >= 2 ){
		rp->cd( 1, 2);
		gPad->SetLogz(1);
		TH2D* s2 = dkl->viewSpecies( 1 );
		s2->GetZaxis()->SetRangeUser( min, max );
		s2->Draw("colz");
	}

	if ( config->getInt( "dklFit.nSpecies" ) >= 3 ){
		rp->cd( 1, 3);
		gPad->SetLogz(1);
		TH2D* s3 = dkl->viewSpecies( 2 );
		s3->GetZaxis()->SetRangeUser( min, max );
		s3->Draw("colz");
	}

	rp->savePage();
}


void pidFitter::runFit(){

	cout << "K_All" << book->get( "nSig_K_All" ) << endl;


}







