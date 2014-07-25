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
	static TH2 *  matrixToHistogram( TMatrixD* m, string name = "hMatrix" );
	
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

	TH2D * viewInput(){
		return (TH2D*)matrixToHistogram( T, "inputData" );
	}
	TH2D * viewApproximation(){
		return (TH2D*)matrixToHistogram( A, "approximation" );
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

	uint nRows, nCols;

	//void createMatrices();
	void updateU( );
	void updateV( );
	void updateA( );
	void update( );



};





#endif