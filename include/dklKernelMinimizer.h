#ifndef DKL_KERNEL_MINIMIZER_H
#define DKL_KERNEL_MINIMIZER_H

#include "allroot.h"
#include "TMatrixD.h"
#include "string.h"
#include "dklMinimizer.h"

class dklKernelMinimizer
{
public:
	dklKernelMinimizer( TH2D * data, uint nSpecies );
	dklKernelMinimizer( TMatrixD data, uint nSpecies );
	~dklKernelMinimizer();

	void run( uint nIterations );
	
	void printResult() {
		cout << "inputData : " << endl;
		T->Print();

		cout << "Approximation : " << endl;
		A->Print();
	}

	void printAll() {
		printResult();
	}
	void printSpecies( uint iS ){
		S[ iS ]->Print();
	}
	void printApproximation() { 
		A->Print();
	}

	TH2D * viewInput(){
		return (TH2D*)dklMinimizer::matrixToHistogram( T, "inputData" );
	}
	TH2D * viewApproximation(){
		return (TH2D*)dklMinimizer::matrixToHistogram( A, "approximation" );
	}
	TH2D * viewSpecies( uint iSpecies );
	//TMatrixD species( uint iSpecies );
	//double speciesYield( uint iSpecies ){
	//	return species( iSpecies ).Sum();
	//}
	double approximationYield() {
		return A->Sum();
	}
	double inputYield() {
		return T->Sum();
	}

	void addKernel( TMatrixD k, uint anchorR, uint anchorC ){
		if ( anchorR >= k.GetNrows() || anchorC >= k.GetNcols() ){
			cout << "invalid kernel " << endl;
			return;
		}

		kernel.push_back( k );
		kernelAnchor a( anchorR, anchorC );
		anchor.push_back( a );

	}

	static TMatrixD makeKernel( uint rows, uint cols, double val = 1.0 );
protected:

	uint nSpecies;
	uint iCurrent;
	TMatrixD *A, *T;
	vector<TMatrixD *> S;
	TH2D* data;

	vector< TMatrixD > kernel;
	
	class kernelAnchor {
	public:
		kernelAnchor( uint r, uint c ){
			row = r;
			col = c;
		}
		uint row;
		uint col;
	};

	vector< kernelAnchor > anchor;

	uint nRows, nCols;

	void createSpecies(  );
	//void updateU( );
	//void updateV( );
	void updateA( );
	inline void updateA( uint r, uint c );
	void updateS( );
	void update( );

	


};





#endif