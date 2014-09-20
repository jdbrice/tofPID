
#include <iostream>
#include "allroot.h"
#include "constants.h"
//#include "histoBook.h"
#include "XmlConfig.h"
#include "pidHistogramMaker.h"
#include "ChainLoader.h"

//#include "pidFitRunner.h"

using namespace jdb;

int main( int argc, char* argv[] ) {

    if ( argc >= 2 ){
        XmlConfig config( argv[ 1 ] );
        config.report();

        string jt = config.getString( "jobType" );
        
        if ( "makeQA" == jt ){
          
            cout << " Making QA plots " << endl;
            TChain * chain = new TChain( "tof" );
            ChainLoader::load( chain, config.getString( "input.dataDir" ).c_str(), config.getInt( "input.dataDir:maxFiles" ) );

            pidHistogramMaker* pid = new pidHistogramMaker( chain, &config  );
            pid->makeQA();
            delete pid;

        } else if ( "histogram" == jt ){
          
            TChain * chain = new TChain( "tof" );
            ChainLoader::load( chain, config.getString( "input.dataDir" ).c_str(), config.getInt( "input.dataDir:maxFiles" ) );

            pidHistogramMaker* pid = new pidHistogramMaker( chain, &config  );
            pid->momentumDistributions();
            //pid->makePidHistograms();
            delete pid;

        } /*else if ( "fit" == jt ) {

            pidFitRunner * pid = new pidFitRunner( &config );
            pid->runFit();
            delete pid;

        }*/
        
    }

	return 0;
}
