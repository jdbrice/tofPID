
#include "dklKernelMinimizer.h"
#include "utils.h"


dklKernelMinimizer::dklKernelMinimizer( TH2D* d, uint nSpecies ){

	cout << "[dklMinimizer." << __FUNCTION__ << "] " << endl;

	this->nSpecies = nSpecies;
	data = d;
	iCurrent = 0;

	nRows = data->GetNbinsY();
	nCols = data->GetNbinsX();

	A = new TMatrixD( nRows, nCols );
	T = dklMinimizer::histogramToMatrix( data );	
	createSpecies();

	updateA();
}

dklKernelMinimizer::dklKernelMinimizer( TMatrixD d, uint nSpecies ){

	cout << "[dklMinimizer." << __FUNCTION__ << "] " << endl;

	this->nSpecies = nSpecies;
	iCurrent = 0;

	nRows = d.GetNrows();
	nCols = d.GetNcols();

	A = new TMatrixD( nRows, nCols );
	T = new TMatrixD( d );	
	createSpecies();

	updateA();
}

dklKernelMinimizer::~dklKernelMinimizer() {
	cout << "[dklMinimizer." << __FUNCTION__ << "] " << endl;

	delete T;
	delete A;
	for ( int i = 0; i < S.size(); i++  )
		delete S[ i ];

}

void dklKernelMinimizer::createSpecies(  ){
	cout << "[dklMinimizer." << __FUNCTION__ << "] " << endl;

	
	for ( int iS = 0; iS < nSpecies; iS++ ){
		
		TMatrixD * sp = new TMatrixD( nRows, nCols );

		// apply random initial condition 
		
		Double_t seed = ((double) clock()) + iS; 	
		sp->Randomize( 1.0, 0.0, seed );
		S.push_back( sp );

	}

}

void dklKernelMinimizer::updateA() {
	

	if ( nSpecies <=  0 )
		return;

	(*A) = *(S[ 0 ]);
	for ( uint iS = 1; iS < nSpecies; iS ++ ){
		
		(*A) += *(S[ iS ]);
	}
}
void dklKernelMinimizer::updateA( uint r, uint c ) {
	
	if ( nSpecies <=  0 )
		return;

	(*A)[ r] [c ] = (*(S[ 0 ]))[ r ][ c ];
	for ( uint iS = 1; iS < nSpecies; iS ++ ){
		
		(*A)[r][c] += (*(S[ iS ]))[r][c];
	}
}


void dklKernelMinimizer::updateS(){
	


	// apply each kernel
	uint lenK = kernel.size();
	

		// loop over the elements in the matrix
		for ( uint iR = 0; iR < nRows; iR ++ ){
			for ( uint iC = 0; iC < nCols; iC ++ ){
				for ( uint iS = 0; iS < nSpecies; iS ++ ){

				for ( uint iK = 0; iK < lenK; iK++ ){

		uint kRows = kernel[ iK ].GetNrows();
		uint kCols = kernel[ iK ].GetNcols();
		uint aR = anchor[ iK ].row;
		uint aC = anchor[ iK ].col;
		TMatrixD * pK = &kernel[ iK ];


				int startRow = (int)iR - (int)aR;
				int startCol = (int)iC - (int)aC;
				int endRow = iR + ( kRows - ( aR + 1 ) );
				int endCol = iC + ( kCols - ( aC + 1 ) );

				if ( startRow < 0 )
					startRow = 0;
				if ( startCol < 0 )
					startCol = 0;
				if ( endRow >= nRows )
					endRow = nRows - 1;
				if ( endCol >= nCols )
					endCol = nCols - 1;

				//cout << " AT [ " << iR << "] [ " << iC << " ] " << endl;
				//cout << "\t" << startRow << " -> " << endRow << endl;
				//cout << "\t" << startCol << " -> " << endCol << endl;
				
				double weightedSum = 0;
				double sum = 0;
				// loop on kernel bounds
				for ( uint kR = startRow; kR <= endRow; kR++ ){
					for ( uint kC = startCol; kC <= endCol; kC++ ){

						// get the weight of the kernel at this point
						double kWeight = (*pK)[ kR - startRow ] [ kC - startCol ];
						// fill the unweighted average
						sum += (*(S[ iS ]))[ kR ][ kC ];
						double num = (*(S[ iS ]))[ kR ][ kC ] * ( (*T)[ kR ][ kC ] );
						double den = ( (*A)[ kR ][ kC ] );
						if ( 0 == den )
							weightedSum += 0;
						else 
							weightedSum += kWeight * ( num / den );
					}
				}

				//for ( uint kR = startRow; kR <= endRow; kR++ ){
				//	for ( uint kC = startCol; kC <= endCol; kC++ ){
						if ( 0 == sum || 0 == weightedSum ){
							//cout << "zero " << endl;
							(*(S[ iS ]))[ iR ][ iC ] = 0;
						} else {
							// update the value here
							double v = (*(S[ iS ]))[ iR ][ iC ];
							(*(S[ iS ]))[ iR ][ iC ] = v * ( weightedSum / sum );
						}
						updateA( iR, iC );
				//	}
				//}
				
				}
				}//loop over species
			} // loop over matrix elements
		}


	

}

void dklKernelMinimizer::update(){

	
		updateS();
		updateA();

	iCurrent ++;

}



void dklKernelMinimizer::run( uint nIterations ) {

	for ( uint i = 0; i < nIterations; i++ ){
		update();
		jdbUtils::progressBar( i, nIterations, 60 );
		//if ( iCurrent % 1000 == 0 )
		//	cout << "[dkl.run] 1000 iterations done : " << iCurrent << endl;
	}

}

TH2D* dklKernelMinimizer::viewSpecies( uint iS ){
	return (TH2D*) dklMinimizer::matrixToHistogram( S[iS] );
}

TMatrixD dklKernelMinimizer::makeKernel( uint r, uint c, double val ){
	TMatrixD m( r, c );
	for ( uint i = 0; i < r; i++  ){
		for ( uint j = 0; j < c; j++  ){
			m[ i ][ j ] = val;
		}
	}
	return m;
}




