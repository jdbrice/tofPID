
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
            pid->makePidHistograms();
            delete pid;

        } else if ( "fit" == jt ) {

            pidFitter * pid = new pidFitter( &config );
            pid->runFit();
            delete pid;

        }
        
    }

	return 0;
}
