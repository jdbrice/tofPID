
#include <iostream>
#include "allroot.h"
#include "constants.h"
//#include "histoBook.h"
#include "XmlConfig.h"
#include "pidHistogramMaker.h"
#include "ChainLoader.h"
#include "SimultaneousFit.h"

#include "Bichsel.h"

//#include "pidFitRunner.h"

using namespace jdb;

int main( int argc, char* argv[] ) {

    if ( argc >= 2 ){
        XmlConfig config( argv[ 1 ] );
        config.report();

        string jt = config.getString( "jobType" );

        Bichsel * bsel = new Bichsel( config.getString( "input.bichselTabels", "dedxBichsel.root" ), 0 );
        cout << bsel->mean( .5, .139 );
        
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
            
            //pid->momentumDistributions();
            pid->makeDedxTofHistograms();
            delete pid;

        } else if ( "fit" == jt ) {

            //pidFitRunner * pid = new pidFitRunner( &config );
            //pid->runFit();
            //delete pid;
        } else if ( "simultaneous" == jt ){

            TFile * in = new TFile( "histogram/nl1k.root", "READ" );
            TH1D * tAll = (TH1D*)in->Get( "tof/tof_a_K_9_0" );

            SimultaneousFit *sf = new SimultaneousFit( 
                                    tAll, NULL, NULL, NULL, 
                                    NULL, NULL, NULL, NULL, 
                                    0.675, &config );

            sf->fitTofAll();
            delete sf;

        }
        
    }

	return 0;
}
