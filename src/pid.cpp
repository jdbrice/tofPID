
#include <iostream>
#include "allroot.h"
#include "constants.h"
#include "histoBook.h"
#include "xmlConfig.h"
#include "pidHistogramMaker.h"
#include "chainLoader.h"



int main( int argc, char* argv[] ) {

  if ( argc >= 2 ){
    xmlConfig config( argv[ 1 ] );
    config.report();


    TChain * chain = new TChain( "tof" );
    chainLoader::load( chain, config.getString( "input.dataDir" ).c_str(), config.getInt( "input.dataDir:maxFiles" ) );

    pidHistogramMaker* pid = new pidHistogramMaker( chain, &config  );

    pid->make();

    delete pid;

  }
	return 0;
}
