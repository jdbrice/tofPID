
#include "multiGaussianFit.h"
#include "utils.h"

using namespace jdbUtils;


multiGaussianFit::multiGaussianFit( TH2* h, uint nS ){
	cout << "[multiGaussianFit." << __FUNCTION__ << "]" << endl;
	input = h;
	nSpecies = nS;

	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
}

multiGaussianFit::~multiGaussianFit() {
	cout << "[multiGaussianFit." << __FUNCTION__ << "]" << endl;

}

void multiGaussianFit::fit() {
	cout << "[multiGaussianFit." << __FUNCTION__ << "]" << endl;
	setupModel();

	model->fitTo( (*rdh), RooFit::PrintLevel(-1) );
	

	//RooGaussian gx("gx","gx",*xObs,RooFit::RooConst(-2),RooFit::RooConst(3)) ;
  	//RooGaussian gy("gy","gy",*yObs,RooFit::RooConst(+2),RooFit::RooConst(2)) ;

  	// Create gxy = gx(x)*gy(y)
  	//gxy = new RooAddPdf("gxy","gxy",RooArgSet(gx,gy)) ;

  	//gxy->fitTo( *rdh ); 

  	//cout << "y Obs " << yObs << endl;
  	//RooPlot * frame = yObs->frame();
  	//gxy.plotOn( frame );
  	//frame->Draw();  
  	//delete frame;

}

void multiGaussianFit::setupModel() {
	cout << "[multiGaussianFit." << __FUNCTION__ << "]" << endl;
	rdh = new RooDataHist( "multiGaussianData", "multiGaussianData", 
				RooArgSet( *xObs, *yObs), input );

	RooArgList *ralGauss = new RooArgList();
	RooArgList *ralAmp = new RooArgList();

	for ( uint iS = 0; iS < nSpecies; iS++ ){

		// create the mean for X and Y observables
		double mX = 0;
		if ( initialMeanX.size() > iS )
			mX = initialMeanX[ iS ];
		double mMinX = xMin, mMaxX = xMax;
		if ( meanMinX.size() > iS && meanMaxX.size() > iS ){
			mMinX = meanMinX[ iS ];
			mMaxX = meanMaxX[ iS ];
		}
		cout << "Mean X " << mX << " in ( " << mMinX << " -> " << mMaxX << " ) " << endl;
		RooRealVar *mx = new RooRealVar( 	("meanX" + ts(iS) ).c_str(), ("meanX" + ts(iS) ).c_str(),
											mX, mMinX, mMaxX );
		meanX.push_back( mx );

		double mY = 0;
		if ( initialMeanY.size() > iS )
			mY = initialMeanY[ iS ];
		RooRealVar *my = new RooRealVar( 	("meanY" + ts(iS) ).c_str(), ("meanY" + ts(iS) ).c_str(),
											mY, yMin, yMax );
		meanY.push_back( my );
		
		// create the sigma for X and Y gaussians
		RooRealVar *sx = new RooRealVar( 	("sigX" + ts(iS) ).c_str(), ("sigX" + ts(iS) ).c_str(),
											.005, 0.0001, .02 );
		sigX.push_back( sx );

		RooRealVar *sy = new RooRealVar( 	("sigY" + ts(iS) ).c_str(), ("sigY" + ts(iS) ).c_str(),
											5, 1, 20 );
		sigY.push_back( sy );
		
		// create the gaussians
		RooGaussian * gx = new RooGaussian( ("gaussX" + ts(iS) ).c_str(), ("gaussX" + ts(iS) ).c_str(),
											*xObs, *mx, *sx );
		gaussX.push_back( gx );

		RooGaussian * gy = new RooGaussian( ("gaussY" + ts(iS) ).c_str(), ("gaussY" + ts(iS) ).c_str(),
											*yObs, *my, *sy );
		gaussY.push_back( gy );

		RooProdPdf * g2D = new RooProdPdf( ("gauss2D" + ts(iS) ).c_str(), ("gauss2D" + ts(iS) ).c_str(),
											RooArgList( *gx, *gy ) );

		// the amplitude of each gaussian
		RooRealVar *a = new RooRealVar( 	("amp" + ts(iS) ).c_str(), ("amp" + ts(iS) ).c_str(),
											0.5, 0, 1 );
		amp.push_back( a );

		// setup the arg lists
		ralGauss->add( *g2D );
		ralAmp->add( *a );

	}

	// now create the final model containing all our 2D gaussians
	model = new RooAddPdf( "multi2DGaussian", "multi2DGaussian", *ralGauss, *ralAmp );

}











