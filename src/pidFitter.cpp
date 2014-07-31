

#include "pidFitter.h"
#include "dklMinimizer.h"
#include "multiGaussianFit.h"


// provides my own string shortcuts etc.
using namespace jdbUtils;

pidFitter::pidFitter( xmlConfig * con ){
	cout << "[pidFitter.pidFitter] " << endl;
	
	gErrorIgnoreLevel=kSysError;

	config = con;

	// set the histogram info verbosity to show nothing
	gStyle->SetOptStat( 0 );
	
	// create the histogram book
	
	book = new histoBook( ( config->getString( "output.base" ) + config->getString( "output.root" ) ), config, config->getString( "input.root:file" ) );	
	lutBook = new histoBook( (  config->getString( "output.lut" ) ), config );	


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
	delete lutBook;
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


void pidFitter::runDkl( TH2D* h, reporter * rp, string optPath ){

	gStyle->SetOptStat( 1 );

	int nS = config->getInt( optPath + ".dkl:nSpecies", 1 );
	int nIt = config->getInt( optPath + ".dkl:nRuns", 1 );
	double angle = config->getDouble( optPath + ".dkl:angle", 0 );

	if ( nS <= 0 )
		return;
	if ( nIt <= 0 )
		return;
	if ( NULL == h || NULL == rp )
		return;

	double axisMin = config->getDouble( "binning.nSig:min", 0 );
	double axisMax = config->getDouble( "binning.nSig:max", 0 );

	dklMinimizer *dkl = new dklMinimizer( h, nS, angle );
	dkl->run( nIt );

	rp->newPage( 1, 2 );

	rp->cd( 1, 1);
	gPad->SetLogz(1);
	dkl->viewInput( )->Draw("colz");

	rp->cd( 2, 1);
	gPad->SetLogz(1);
	//TH2D* ap = dkl->viewApproximation( axisMin, axisMax, axisMin, axisMax );
	TH2D* ap = dkl->viewApproximation( );
	double min = 0.01;//ap->GetMinimum() * 0.1;
	double max = ap->GetMaximum();
	ap->Draw("colz");
	rp->savePage();

	rp->newPage( 1, 3 );
	rp->cd( 1, 1);
	gPad->SetLogz(1);

	TH2D* s1 = NULL, *s2 = NULL, *s3 = NULL;
	s1 = dkl->viewSpecies( 0 );
	s1->GetZaxis()->SetRangeUser( min, max );
	s1->Draw("colz");

	if ( nS >= 2 ){
		rp->cd( 1, 2);
		gPad->SetLogz(1);
		s2 = dkl->viewSpecies( 1 );
		s2->GetZaxis()->SetRangeUser( min, max );
		s2->Draw("colz");
	}

	if ( nS >= 3 ){
		rp->cd( 1, 3);
		gPad->SetLogz(1);
		s3 = dkl->viewSpecies( 2 );
		s3->GetZaxis()->SetRangeUser( min, max );
		s3->Draw("colz");
	}

	rp->savePage();

	rp->newPage();
	uint ciS = dklMinimizer::speciesClosestTo( 0, 0, s1, s2, s3 );
	//TH2 * cSpecies = dklMinimizer::viewSpeciesClosestTo( 0, 0, s1, s2, s3 );
	//TH2D* ps1 = (TH2D*)cSpecies->Clone( "probabilityS1" );
	TH2D * ps1 = (TH2D*) dkl->speciesProbabilityMap( ciS );
	cout << " Input Yield: " << dkl->inputYield( ) << endl;;
	cout << " Species Yield: " << dkl->speciesYield( ciS ) << endl;;
	cout << " Total Approx Yield: " << dkl->approximationYield( ) << endl;;

	makeSquareCuts( ps1, "K_Fit.dklPostFitCut" );
	//ps1->Divide( ap );

	ps1->Draw( "colz" );
	ps1->GetZaxis()->SetRangeUser( 0, 1 );

	rp->savePage();

	//sDelete( s1 );
	//sDelete( s2 );
	//sDelete( s3 );
}


void pidFitter::runFit(){

	lutBook->cd("");
	lutBook->makeAll( "hist" );

	book->cd("");
	processSpecies( "Pi", 0, report );
	processSpecies( "K", 0, report );

}


void pidFitter::processSpecies( string species, int charge, reporter * rp ){
	cout << "[pidFitter." << __FUNCTION__ << "]" << endl;

	string hName = "nSig_" + sName( species, charge );
	// get the pt Binning
	vector<double>pBins = config->getDoubleVector( "binning.pBins" );

	string useNode = "";

	TH3* h3 = book->get3D( hName );
	int nFits = config->getInt( species + "_Fit:nFits", 1 );
	cout << "Number of Fit Categories: " << nFits << endl;
	for ( int iFit = 1; iFit <= nFits; iFit++ ){

		string optPath = species + "_Fit.opt" + config->getString( species + "_Fit.fit"+ts(iFit)+":options");
		int fBin = config->getInt( species + "_Fit.fit"+ts(iFit)+":min", 1 );
		int lBin = config->getInt( species + "_Fit.fit"+ts(iFit)+":max", pBins.size() );	

		cout << "Fitting P bins ( " << fBin << " --> " << lBin << " ) " << endl;
		for ( int i = fBin; i <= lBin; i ++ ){

			// get the Pt range for title etc.
			double pLow = h3->GetZaxis()->GetBinLowEdge( i );
			double pHi = h3->GetZaxis()->GetBinLowEdge( i + 1 );	

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
			//nProj->Draw("colz");

			// process the square cuts
			makeSquareCuts( cut, optPath + ".squareCut" );

			// draw the distribution after square cuts
			rp->cd( 1, 2 );
			gPad->SetLogz( 1 );
			cut->SetTitle( ( "After 1D Cuts : " + ts( pLow, 4 ) + " #leq " + " P #leq" + ts( pHi, 4 ) ).c_str()  );
			cut->Draw( "colz" );
			rp->savePage( );

			// fit using the dkl algorithm
			if ( config->nodeExists( optPath + ".dkl" ) )
				runDkl( cut, rp, optPath );
			// runs the 2d gaussian fit
			else if ( config->nodeExists( optPath + ".mgf" ) )
				runMultiGauss( cut, rp, species, optPath, i );
			
			
		}

	} 
	
	
}

void pidFitter::runMultiGauss( TH2D* h, reporter* rp, string pType, string nodePath, uint pBin ){

	h->Sumw2();
	uint nS = config->getInt( nodePath + ".mgf:nSpecies", 1 );

	double axisMin = config->getDouble( "binning.nSig:min", 0 );
	double axisMax = config->getDouble( "binning.nSig:max", 0 );

	multiGaussianFit * mgf = new multiGaussianFit( h, nS );
	mgf->setX( "#Delta#beta^{-1}/#beta^{-1}",  axisMin, axisMax );
	mgf->setY( "n#sigma dedx",  axisMin, axisMax );

	for ( int i = 1; i < nS + 1; i ++ ){

		double mX = config->getDouble( nodePath + ".mgf.initialMean:x" + ts( i ), 0);
		double mY = config->getDouble( nodePath + ".mgf.initialMean:y" + ts( i ), 0);
		mgf->setInitialMean( mX, mY );

		double minX = config->getDouble( nodePath + ".mgf.meanLimits:min" + ts( i ), -999);
		double maxX = config->getDouble( nodePath + ".mgf.meanLimits:max" + ts( i ), -999);
		if ( minX > -998 && maxX > -998 )
			mgf->limitMeanX( minX, maxX );		

	}

	
	mgf->fit();

	lutBook->cd("");
	lutBook->get( "xMean" + pType )->SetBinContent( pBin, mgf->getMeanX( 0 ) );
	lutBook->get( "xMean" + pType )->SetBinError( pBin, mgf->getMeanXError( 0 ) );
	lutBook->get( "yMean" + pType )->SetBinContent( pBin, mgf->getMeanY( 0 ) );
	lutBook->get( "yMean" + pType )->SetBinError( pBin, mgf->getMeanYError( 0 ) );

	// sigmas
	lutBook->get( "xSigma" + pType )->SetBinContent( pBin, mgf->getSigmaX( 0 ) );
	lutBook->get( "xSigma" + pType )->SetBinError( pBin, mgf->getSigmaXError( 0 ) );
	lutBook->get( "ySigma" + pType )->SetBinContent( pBin, mgf->getSigmaY( 0 ) );
	lutBook->get( "ySigma" + pType )->SetBinError( pBin, mgf->getSigmaYError( 0 ) );

	
	rp->newPage();
	mgf->viewFitX( "Fit");
	gPad->SetLogz(1);
	rp->savePage();

	delete mgf;

}


double pidFitter::squareCut( TH2D* h, string axis, string cut, double value ){

	if ( "x" != axis && "y" != axis && "X" != axis && "Y" != axis )
		return 0;
	if ( ">" != cut && "<" != cut )
		return 0;
	if ( NULL == h)
		return 0;

	//cout << axis  << cut << value << endl;

	int nX = h->GetNbinsX();
	int nY = h->GetNbinsY();

	double nRemoved = 0;

	for ( int bX = 1; bX <= nX; bX++ ){
		for ( int bY = 1; bY <= nY; bY++ ){

			double xEdge = h->GetXaxis()->GetBinLowEdge( bX );
			double yEdge = h->GetYaxis()->GetBinLowEdge( bY );

			if ( "x" == axis ){
				if ( xEdge > value && ">" == cut ){
					nRemoved += h->GetBinContent( bX, bY );
					h->SetBinContent( bX, bY, 0.0);
					h->SetBinError( bX, bY, 0);
				} else if ( xEdge < value && "<" == cut ){
					nRemoved += h->GetBinContent( bX, bY );
					h->SetBinContent( bX, bY, 0.0);
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




