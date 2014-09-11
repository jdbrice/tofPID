

#include "dklMinimizer.h"
#include <time.h> 
#include "jdbUtils.h"
using namespace jdbUtils;

dklMinimizer::dklMinimizer( TH2D* d, uint nSpecies, double rotate ) {
	cout << "[dklMinimizer." << __FUNCTION__ << "] " << endl;

	this->nSpecies = nSpecies;
	data = d;
	iCurrent = 0;

	// store the number of bins
	nBinsX = data->GetNbinsX();
	nBinsY = data->GetNbinsY();

	// now trim the matrix to make computation faster
	T = histogramToMatrix( data, false /* force trimming */ );

	if ( rotate != 0 ){
		rot = rotate;
		cout << "Rotating Matrix " << endl;
		TMatrixD* temp = T;
		T = rotateMatrix( (*T), rotate );
		delete temp;
	}


	// do the trimming manually to get the User Space min /max
	int x1 = 0, x2 = 0, y1 = 0, y2 = 0;
	trimHistogram( d, &x1, &x2, &y1, &y2 );

	minX = d->GetXaxis()->GetBinLowEdge( x1 );
	minY = d->GetYaxis()->GetBinLowEdge( y1 );

	maxX = d->GetXaxis()->GetBinLowEdge( x2 ) + d->GetXaxis()->GetBinWidth( x2 );
	maxY = d->GetYaxis()->GetBinLowEdge( y2 ) + d->GetYaxis()->GetBinWidth( y2 );

	nRows = T->GetNrows();
	nCols = T->GetNcols();

	A = new TMatrixD( nRows, nCols );
	U = new TMatrixD( nRows, nSpecies );
	V = new TMatrixD( nSpecies, nCols );
	

	// initial conditions
	Double_t seed = (double) clock();
	U->Randomize( 1.0, 0.0, seed );
	V->Randomize( 1.0, 0.0, seed );
	

	updateA();

}


dklMinimizer::~dklMinimizer() {

	delete A;
	delete T;
	delete U;
	delete V;

}

TMatrixD * dklMinimizer::histogramToMatrix( TH2* h, bool trim) {
	cout << "[dklMinimizer." << __FUNCTION__ << "] " << endl;


	if ( !trim ){

		uint rows = h->GetNbinsY();
		uint cols = h->GetNbinsX();

		TMatrixD * m = new TMatrixD( rows, cols );

		for ( int x = 0; x < cols; x ++ ){
			for ( int y = 0; y < rows; y ++ ){
				(*m)[ y ][ x ] = h->GetBinContent( x+1, y+1 );
			}
		}

		return m;
	} else {

		int x1 = 0, x2 = 0, y1 = 0, y2 = 0;
		trimHistogram( h, &x1, &x2, &y1, &y2 );

		uint rows = y2 - y1 + 1;
		uint cols = x2 - x1 + 1;

		TMatrixD * m = new TMatrixD( rows, cols );

		for ( int x = x1; x <= x2; x++ ){
			for ( int y = y1; y <= y2; y++ ){
				(*m)[ y - y1 ][ x - x1 ] = h->GetBinContent( x+1, y+1 );
			}
		}
		return m;
	}

}

void dklMinimizer::trimHistogram( TH2* d, int* xBin1, int* xBin2, int* yBin1, int* yBin2, double thresh ){

	try{

		if ( NULL == xBin1 || NULL == xBin2 || NULL == yBin1 || NULL == yBin2 )
			throw 0;

		// first do X axis
		int nBinsX = d->GetNbinsX();
		int nBinsY = d->GetNbinsY();

		int x1 = 1;
		for ( int x = 1; x < nBinsX; x++ ){
			
			bool empty = true;
			for ( int y = 1; y < nBinsY; y ++ ){
				double v = d->GetBinContent( x, y );
				if ( v > thresh ){
					empty = false ;
					break;
				}
			}

			if ( empty == true )
				x1 = x;
			else 
				break;
		}

		// find x2
		int x2 = 1;
		for ( int x = nBinsX; x > 1; x-- ){
			
			bool empty = true;
			for ( int y = 1; y < nBinsY; y ++ ){
				double v = d->GetBinContent( x, y );
				if ( v > thresh ){
					empty = false;
					break;
				}
			}

			if ( empty == true )
				x2 = x;
			else 
				break;
		}

		//find y1
		int y1 = 1;
		for ( int y = 1; y < nBinsY; y++ ){
			
			bool empty = true;
			for ( int x = 1; x < nBinsX; x++ ){
				double v = d->GetBinContent( x, y );
				if ( v > thresh ){
					empty = false ;
					break;
				}
			}

			if ( empty == true )
				y1 = y;
			else 
				break;
		}

		// find y2
		int y2 = 1;
		for ( int y = nBinsY; y > 1; y-- ){
			
			bool empty = true;
			for ( int x = 1; x < nBinsX; x++ ){
				double v = d->GetBinContent( x, y );
				if ( v > thresh ){
					empty = false ;
					break;
				}
			}

			if ( empty == true )
				y2 = y;
			else 
				break;
		}

		//cout << "effective Bin Range X : ( " << x1 << ", " << x2 << " ) Y : ( " << y1 << ", " << y2 << " ) " << endl;
		// copy the values to the output pointers
		
		(*xBin1) = x1;
		(*xBin2) = x2;
		(*yBin1) = y1;
		(*yBin2) = y2;


	} catch ( ... ){
		cout << "could not trim histogram " << endl;
	}

	return;
}

TH2 * dklMinimizer::matrixToHistogram( TMatrixD * m, string name, double x1, double x2, double y1, double y2 ) {
	cout << "[dklMinimizer." << __FUNCTION__ << "] " << endl;

	uint rows = m->GetNrows();
	uint cols = m->GetNcols();

	double xBin1 = x1, xBin2 = x2;
	double yBin1 = y1, yBin2 = y2;

	if ( xBin1 - xBin2  > 0 ||
		 xBin2 - xBin1  <= numeric_limits<double>::epsilon()*3 ){
		xBin1 = 0;
		xBin2 = cols;
	}
	if ( yBin1 - yBin2  > 0 ||
		 yBin2 - yBin1  <= numeric_limits<double>::epsilon()*3 ){
		yBin1 = 0;
		yBin1 = rows; 
	}

	TH2D * h = new TH2D( name.c_str(), name.c_str(), cols, xBin1, xBin2, rows, yBin1, yBin2 );

	for ( int x = 0; x < cols; x ++ ){
		for ( int y = 0; y < rows; y ++ ){
			double v = (*m)[ y ][ x ];
			h->SetBinContent( x+1, y+1, v );
		}
	}
	return h;
}

TH2 * dklMinimizer::viewSpeciesClosestTo( double x, double y, TH2 * s1, TH2 * s2, TH2 * s3, TH2 * s4 ){

	double r1 = distanceFrom( x, y, s1);
	double r2 = distanceFrom( x, y, s2);
	double r3 = distanceFrom( x, y, s3);
	double r4 = distanceFrom( x, y, s4);

	if ( r1 == -1 )
		r1 = std::numeric_limits<double>::max();
	if ( r2 == -1 )
		r2 = std::numeric_limits<double>::max();
	if ( r3 == -1 )
		r3 = std::numeric_limits<double>::max();
	if ( r4 == -1 )
		r4 = std::numeric_limits<double>::max();

	double min = std::min( {r1, r2, r3, r4} );

	if ( r1 == min)
		return s1;
	if ( r2 == min)
		return s2;
	if ( r3 == min)
		return s3;

	cout << "def" << endl;
	return s4;


}

uint dklMinimizer::speciesClosestTo( double x, double y, TH2 * s1, TH2 * s2, TH2 * s3, TH2 * s4 ){

	double r1 = distanceFrom( x, y, s1);
	double r2 = distanceFrom( x, y, s2);
	double r3 = distanceFrom( x, y, s3);
	double r4 = distanceFrom( x, y, s4);

	if ( r1 == -1 )
		r1 = std::numeric_limits<double>::max();
	if ( r2 == -1 )
		r2 = std::numeric_limits<double>::max();
	if ( r3 == -1 )
		r3 = std::numeric_limits<double>::max();
	if ( r4 == -1 )
		r4 = std::numeric_limits<double>::max();

	double min = std::min( {r1, r2, r3, r4} );

	if ( r1 == min)
		return 0;
	if ( r2 == min)
		return 1;
	if ( r3 == min)
		return 2;

	cout << "def" << endl;
	return 3;

}

TMatrixD dklMinimizer::probabilityMap( TMatrixD s, TMatrixD total ){

	uint rows = total.GetNrows();
	uint cols = total.GetNcols();

	TMatrixD pMap( rows, cols );
	for ( uint r = 0; r < rows; r++ ){
		for ( uint c = 0; c < cols; c++ ){
			
			double tV = total[ r ][ c ];
			double sV = s[ r ][ c ];

			if ( 0 == tV  )
				pMap[ r ][ c ] = 0;
			else
				pMap[ r ][ c ]= sV / tV;
		}
	}

	return pMap;
}




inline void dklMinimizer::updateU( ) {

	for ( int i = 0; i < nRows; i++ ){
		for ( int p = 0; p < nSpecies; p++ ){

			double num = 0;
			double den = 0;
			double orig = (*U)[ i ][ p ];

			// loop through a row of V
			for ( int a = 0; a < nCols; a++ ){
				double vT = (*T)[ i ][ a ];
				double vA = (*A)[ i ][ a ];
				double res = (vT / vA );
				
				if ( 0 == vA && 0 == vT )
					res = 0;
				else if ( 0 == vA )
					res = 0;

				num += (*V)[ p ][ a ] * res;
				den += (*V)[ p][ a ];
			}

			if ( den == 0 ){
				(*U)[i][p] = 0;
			} else {
				(*U)[ i ][ p ] = orig * ( num / den );
			}
		}
	}
}

inline void dklMinimizer::updateV( ) {

	for ( int p = 0; p < nSpecies; p++ ){
		for ( int j = 0; j < nCols; j++ ){
			
			double num = 0;
			double den = 0;
			double orig = (*V)[ p ][ j ];

			// loop through a col of U
			for ( int a = 0; a < nRows; a++ ){
				double vT = (*T)[ a ][ j ];
				double vA = (*A)[ a ][ j ];
				double res = (vT / vA ) ;
				
				if ( 0 == vA && 0 == vT )
					res = 0;
				else if ( 0 == vA )
					res = 0;

				num += (*U)[ a ][ p ] * res;
				den += (*U)[ a ][ p ];
			}

			if ( den == 0 ){
				(*V)[ p ][ j ] = 0;
			} else {
				(*V)[ p ][ j ] = orig * ( num / den );
			}
		}
	}
}

inline void dklMinimizer::updateA( ){
	(*A) = (*U) * (*V);
}

void dklMinimizer::update( ) {

	if ( iCurrent % 2 == 0){
		updateU();
		updateA();
		updateV();
		updateA();
	} else {
		updateV();
		updateA();
		updateU();
		updateA();
	}

	iCurrent ++;

}

void dklMinimizer::run( uint nIterations ) {

	for ( uint i = 0; i < nIterations; i++ ){
		jdbUtils::progressBar( i, nIterations, 60 );
		update();
	}

}

TMatrixD dklMinimizer::species( uint iSpecies ){
	TMatrixD cU = U->GetSub( 0, nRows-1, iSpecies, iSpecies );
	TMatrixD cV = V->GetSub( iSpecies, iSpecies, 0, nCols - 1 );
	
	TMatrixD  cA = cU * cV;
	return cA;
}
TH2D * dklMinimizer::viewSpecies( uint iSpecies ) {
	TMatrixD cA = species( iSpecies );
	return (TH2D*) matrixToHistogram( &cA, "species", minX, maxX, minY, maxY );
}


TH2 * dklMinimizer::speciesProbabilityMap( uint iS ){

	TMatrixD pMap = probabilityMap( species( iS ), (*A) );
	return (TH2*) matrixToHistogram( &pMap, "speciesProbMap", minX, maxX, minY, maxY );

}

TMatrixD * dklMinimizer::rotateMatrix( TMatrixD m, double by ){

	uint rows = m.GetNrows();
	uint cols = m.GetNcols();

	double minX = 1000;
	double minY = 1000;
	double maxX = 0;
	double maxY = 0;

	for ( uint r = 0; r < rows; r++ ){
		for ( uint c = 0; c < cols; c++ ){

			double y = r;
			double x = c;

			double xp = x * TMath::Cos( by ) + y * TMath::Sin( by );
			double yp = y * TMath::Cos( by ) - x * TMath::Sin( by );

			if ( xp < minX )
				minX = xp;
			if ( yp < minY )
				minY = yp;
			if ( yp > maxY )
				maxY = yp;
			if ( xp > maxX )
				maxX = xp;

		}
	}

	int nR = maxY - minY;
	int nC = maxX - minX;
	
	TMatrixD * nM = new TMatrixD( nR, nC );

	for ( uint r = 0; r < rows; r++ ){
		for ( uint c = 0; c < cols; c++ ){

			double y = r;
			double x = c;

			double xp = x * TMath::Cos( by ) + y * TMath::Sin( by );
			double yp = y * TMath::Cos( by ) - x * TMath::Sin( by );
			int newR = (int)(yp - minY);
			int newC = (int)(xp - minX);
			
			if ( newR < nR && newR >= 0 && newC < nC && newC >= 0)
			(*nM)[ newR ][ newC ] = m[ r ][ c ];

		}
	}

	return nM;

}
