

#include "../include/pidSquare.h"

#include "sstream"
#include "string"
#include "TStyle.h"
#include "TCanvas.h"

using namespace std;

void square( string inFile = "allSim.root", int cSpecies = 1 ){

	TCanvas * c = new TCanvas( "c", "c", 800, 800 );
	c->Print( "bPID.pdf[");

	TFile * f = new TFile( inFile.c_str(), "READ" );


	c->Divide( 2, 2 );
	for ( int i = 0; i < 50; i++ ){

		stringstream sstr;
		sstr << "h_dedx_tof_p3_b" << i;
		TH2* sum = (TH2*)f->Get( sstr.str().c_str() );
		
		sstr.str("");
		sstr << "h_dedx_tof_p0_b" << i;

		TH2* p0 = (TH2*)f->Get( sstr.str().c_str() );
		sstr.str("");
		sstr << "h_dedx_tof_p1_b" << i;
		TH2* p1 = (TH2*)f->Get( sstr.str().c_str() );
		sstr.str("");
		sstr << "h_dedx_tof_p2_b" << i;
		TH2* p2 = (TH2*)f->Get( sstr.str().c_str() );

		pidSquare * pid = new pidSquare( sum, cSpecies, 1.0, p0, p1, p2 );

		c->cd( 3 );
		gPad->SetLogz(1);
		sum->Draw("colz");
		

		c->cd( 4 );
		gPad->SetLogx(1);
		TH1* pX = sum->ProjectionX();
		TH1* pY = sum->ProjectionY();
		pY->SetFillColor( kBlue );
		pY->Draw("hbar");


		c->cd(1);
		gPad->SetLogy(1);
		pX->SetFillColor( kBlue );
		pX->Draw();
		c->cd(2);
		

		TH1D* eff = pid->efficiency( 0.0, 5.0, 0.1 );
		//TH1D* pure = pid->purity( 0.0, 5.0, 0.1 );

		gStyle->SetOptStat( 0 );
		eff->SetTitle( "Efficiecy (Blue), Purity (Red)" );
		eff->GetYaxis()->SetRangeUser(0, 1.05);
		eff->SetLineWidth( 2 );
		eff->Draw();
		//pure->SetLineColor( kRed );
		//pure->SetLineWidth( 2 );
		//pure->Draw("same");

		c->Print( "bPID.pdf");
	}


	c->Print( "bPID.pdf]");


}