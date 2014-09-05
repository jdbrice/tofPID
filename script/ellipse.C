

#include "../include/pidEllipse.h"

#include "sstream"
#include "string"

#include "TStyle.h"
#include "TCanvas.h"

using namespace std;

void ellipse( string inFile = "allSim.root", int cSpecies = 1 ){

	TCanvas * c = new TCanvas( "c", "c", 800, 800 );
	c->Print( "bPID.pdf[");

	TFile * f = new TFile( inFile.c_str(), "READ" );


	c->Print( "bPID.pdf]");


}