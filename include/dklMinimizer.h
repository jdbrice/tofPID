#ifndef DKL_MINIMIZER_H
#define DKL_MINIMIZER_H

#include "allroot.h"
#include "TMatrixD.h"
#include "string.h"

class dklMinimizer
{
public:
	dklMinimizer( TH2D * data, uint nIterations );
	~dklMinimizer();

	void run( uint nIterations );

	static TMatrixD *  histogramToMatrix( TH2* data );
	static TH2 *  matrixToHistogram( TMatrixD* m, string name = "hMatrix", double x1 = 0, double x2 = 0, double y1 = 0, double y2 = 0 );
	
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

	TH2D * viewInput( double x1 = 0, double x2 = 0, double y1 = 0, double y2 = 0 ){
		return (TH2D*)matrixToHistogram( T, (string)"inputData", x1, x2, y1, y2 );
	}
	TH2D * viewApproximation( double x1 = 0, double x2 = 0, double y1 = 0, double y2 = 0 ){
		return (TH2D*)matrixToHistogram( A, (string)"approximation", x1, x2, y1, y2  );
	}
	TH2D * viewSpecies( uint iSpecies, double x1 = 0, double x2 = 0, double y1 = 0, double y2 = 0 );
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

	uint nRows, nCols;

	//void createMatrices();
	void updateU( );
	void updateV( );
	void updateA( );
	void update( );



};





#endif