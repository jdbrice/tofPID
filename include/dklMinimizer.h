#ifndef DKL_MINIMIZER_H
#define DKL_MINIMIZER_H

#include "allroot.h"
#include "TMatrixD.h"
#include "string.h"
#include <limits>

class dklMinimizer
{
public:
	dklMinimizer( TH2D * data, uint nIterations, double rotate = 0);
	~dklMinimizer();

	void run( uint nIterations );

	static TMatrixD *  histogramToMatrix( TH2* data, bool trim = false  );
	static void trimHistogram( TH2* data, int* xBin1, int* xBin2, int* yBin1, int* yBin2, double threshold = 0 );
	static TH2 *  matrixToHistogram( TMatrixD* m, string name = "hMatrix", double x1 = 0, double x2 = 0, double y1 = 0, double y2 = 0 );
	static TH2 * viewSpeciesClosestTo( double x, double y, TH2 * s1, TH2 * s2, TH2 * s3 = NULL, TH2 * s4 = NULL ); 
	static uint speciesClosestTo( double x, double y, TH2 * s1, TH2 * s2, TH2 * s3 = NULL, TH2 * s4 = NULL ); 
	static inline double distanceFrom( double x, double y, TH2 * h ){
		if ( NULL != h ){
			double mX = h->GetMean( 1 );
			double mY = h->GetMean( 2 );
			return ( TMath::Sqrt( mX*mX + mY*mY ) );
		}
		return -1;
	} 

	static TMatrixD probabilityMap( TMatrixD species, TMatrixD total );
	static TMatrixD * rotateMatrix( TMatrixD m, double by );


	TH2 * speciesProbabilityMap( uint iS );

	uint speciesClosestTo( double x, double y ); 

	void printResult() {
		cout << "inputData : " << endl;
		T->Print();

		cout << "Approximation : " << endl;
		A->Print();
	}

	void printAll() {
		cout << "U" << endl;
		U->Print();
		cout << "V" << endl;
		V->Print();
		printResult();
	}

	TH2D * viewInput( ){
		return (TH2D*)matrixToHistogram( T, (string)"inputData", minX, maxX, minY, maxY );
	}
	TH2D * viewApproximation(  ){
		return (TH2D*)matrixToHistogram( A, (string)"approximation", minX, maxX, minY, maxY ); 
	}
	TH2D * viewSpecies( uint iSpecies );
	TMatrixD species( uint iSpecies );
	double speciesYield( uint iSpecies ){
		return species( iSpecies ).Sum();
	}
	double approximationYield() {
		return A->Sum();
	}
	double inputYield() {
		return T->Sum();
	}
protected:

	uint nSpecies;
	uint iCurrent;
	TMatrixD *A, *T, *U, *V;
	TH2D* data;

	double rot;

	uint nRows, nCols;
	uint nBinsX, nBinsY;
	double minX, minY, maxX, maxY;

	//void createMatrices();
	void updateU( );
	void updateV( );
	void updateA( );
	void update( );



};





#endif