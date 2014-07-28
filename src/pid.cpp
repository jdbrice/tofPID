
#include <iostream>
#include "allroot.h"
#include "constants.h"
#include "histoBook.h"
#include "xmlConfig.h"
#include "pidHistogramMaker.h"
#include "chainLoader.h"
#include "pidFitter.h"



int main( int argc, char* argv[] ) {

    if ( argc >= 2 ){
        xmlConfig config( argv[ 1 ] );
        config.report();

        string jt = config.getString( "jobType" );
        
        if ( "makeQA" == jt ){
          
            cout << " Making QA plots " << endl;
            TChain * chain = new TChain( "tof" );
            chainLoader::load( chain, config.getString( "input.dataDir" ).c_str(), config.getInt( "input.dataDir:maxFiles" ) );

            pidHistogramMaker* pid = new pidHistogramMaker( chain, &config  );
            pid->makeQA();
            delete pid;

        } else if ( "histogram" == jt ){
          
            TChain * chain = new TChain( "tof" );
            chainLoader::load( chain, config.getString( "input.dataDir" ).c_str(), config.getInt( "input.dataDir:maxFiles" ) );

            pidHistogramMaker* pid = new pidHistogramMaker( chain, &config  );
            pid->make();
            delete pid;

        } else if ( "fit" == jt ) {

            pidFitter * pid = new pidFitter( &config );
            pid->runFit();
            delete pid;

        }
        
    }





  //reporter * rp = new reporter( "kernel.pdf" );

  //TH2D* h = new TH2D( "dan", "dan", 10, 0, 10, 10, 0, 10 );
  
    /*for (Int_t i = 0; i < 500000; i++) {
      double px = gRandom->Gaus( 25, 3 );
      double py = gRandom->Gaus( 25, 3 );
      h->Fill( px, py );
    }
    for (Int_t i = 0; i < 500000; i++) {
      double px = gRandom->Gaus( 15, 3 );
      double py = gRandom->Gaus( 15, 3 );
      h->Fill( px, py );
    }*/
    /*for (Int_t i = 0; i < 500000; i++) {
      double px = gRandom->Gaus( 0, 3 );
      double py = gRandom->Gaus( 5, 3 );
      h->Fill( px, py );
    }
    for (Int_t i = 0; i < 500000; i++) {
      double px = gRandom->Gaus( 9, 3 );
      double py = gRandom->Gaus( 5, 3 );
      h->Fill( px, py );
    }*/
/*
      gStyle->SetOptStat( 0 );
  TMatrixD d( 10, 10 );

  d[ 0 ][ 0 ] = 10;
  d[ 0 ][ 1 ] = 10;
  d[ 1 ][ 0 ] = 10;
  d[ 1 ][ 1 ] = 10;

  d[ 9 ][ 9 ] = 10;
  d[ 9 ][ 8 ] = 10;
  d[ 8 ][ 9 ] = 10;
  d[ 8 ][ 8 ] = 10;
  dklKernelMinimizer dkl( d, 2 );

  TMatrixD k = dklKernelMinimizer::makeKernel( 1, 50 );
  //k[ 1 ][ 1 ] = 2;
  dkl.addKernel( k, 0, 25 );
  k = dklKernelMinimizer::makeKernel( 50, 1 );
  dkl.addKernel( k, 25, 0 );


  rp->newPage( 1, 2 );
  dkl.viewInput()->Draw( "colz" );
  rp->cd( 1, 2 );
  dkl.viewApproximation()->Draw( "colz" );
  rp->savePage();

  rp->newPage( 1, 3 );
  dkl.viewSpecies( 0 )->Draw( "colz" );
  rp->cd( 1, 2 );
  dkl.viewSpecies( 1 )->Draw( "colz" );
  //rp->cd( 1, 3 );
  //dkl.viewSpecies( 2 )->Draw( "colz" );
  rp->savePage();


  dkl.run( 500 );

  rp->newPage( 1, 2 );
  dkl.viewInput()->Draw( "colz" );
  rp->cd( 1, 2 );
  dkl.viewApproximation()->Draw( "colz" );
  rp->savePage();

  rp->newPage( 1, 3 );
  dkl.viewSpecies( 0 )->Draw( "colz" );
  rp->cd( 1, 2 );
  dkl.viewSpecies( 1 )->Draw( "colz" );
  //rp->cd( 1, 3 );
  //dkl.viewSpecies( 2 )->Draw( "colz" );
  rp->savePage();

  delete rp;
*/
	return 0;
}
