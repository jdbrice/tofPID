
#include "vector"
#include "TFile.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TCanvas.h"

using namespace std;

void epCompare( string f1, string f2, string f3 = "", string f4 = "" ){

	TCanvas * c = new TCanvas( "c", "c", 800, 800 );
	string outName = "epCompare.pdf";
	c->Print( (outName+"[").c_str() );

	vector<string> fName;
	fName.push_back( f1 );
	fName.push_back( f2 );
	if ( "" != f3 )
		fName.push_back( f3 );
	if ( "" != f4 )
		fName.push_back( f4 );

	vector<TFile *> files;
	for ( int i = 0; i < fName.size(); i++ ){
		files.push_back( new TFile( fName[i].c_str(), "READ" ) );
	}
	
	TLegend * leg = new TLegend( 0.7, 0.7, 0.9, 0.9 );
	for ( int i = 0; i < files.size(); i++ ){

		TH1D* eff = (TH1D*)(files[ i ])->Get( "hEffVsP");
		string draw = "";

		eff->SetLineColor( i+1 );

		if ( i >= 1 )
			draw = "same";

		eff->Draw( draw.c_str() );

		leg->AddEntry( eff, fName[i].c_str(), "lpf" );
	}
	leg->Draw();
	c->Print( (outName).c_str() );

	delete leg;
	leg = new TLegend( 0.7, 0.7, 0.9, 0.9 );
	for ( int i = 0; i < files.size(); i++ ){

		TH1D* pure = (TH1D*)(files[ i ])->Get( "hPureVsP");
		string draw = "";

		pure->SetLineColor( i+1 );

		if ( i >= 1 )
			draw = "same";

		pure->Draw( draw.c_str() );

		leg->AddEntry( pure, fName[i].c_str(), "lpf" );
	}
	leg->Draw();
	c->Print( (outName).c_str() );


	c->Print( (outName+"]").c_str() );



	return;
}