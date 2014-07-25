

#include "dklMinimizer.h"
#include <time.h> 

dklMinimizer::dklMinimizer( TH2D* d, uint nSpecies ) {
	cout << "[dklMinimizer." << __FUNCTION__ << "] " << endl;

	this->nSpecies = nSpecies;
	data = d;
	iCurrent = 0;

	nRows = data->GetNbinsY();
	nCols = data->GetNbinsX();

	A = new TMatrixD( nRows, nCols );
	U = new TMatrixD( nRows, nSpecies );
	V = new TMatrixD( nSpecies, nCols );
	
	T = histogramToMatrix( data );

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

TMatrixD * dklMinimizer::histogramToMatrix( TH2* h ) {
	cout << "[dklMinimizer." << __FUNCTION__ << "] " << endl;

	uint rows = h->GetNbinsY();
	uint cols = h->GetNbinsX();

	TMatrixD * m = new TMatrixD( rows, cols );

	for ( int x = 0; x < cols; x ++ ){
		for ( int y = 0; y < rows; y ++ ){
			(*m)[ y ][ x ] = h->GetBinContent( x, y );
		}
	}

	return m;

}

TH2 * dklMinimizer::matrixToHistogram( TMatrixD * m, string name ) {
	cout << "[dklMinimizer." << __FUNCTION__ << "] " << endl;

	uint rows = m->GetNrows();
	uint cols = m->GetNcols();

	TH2D * h = new TH2D( name.c_str(), name.c_str(), cols, 0, cols, rows, 0, rows );

	for ( int x = 0; x < cols; x ++ ){
		for ( int y = 0; y < rows; y ++ ){
			double v = (*m)[ y ][ x ];
			h->SetBinContent( x, y, v );
		}
	}
	return h;
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
				double res = (vT / vA );
				
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
		update();
		if ( iCurrent % 1000 == 0 )
			cout << "[dkl.run] 1000 iterations done : " << iCurrent << endl;
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
	return (TH2D*) matrixToHistogram( &cA, "species" );
}

