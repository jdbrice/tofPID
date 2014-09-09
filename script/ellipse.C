

#include "../include/pidEllipse.h"

#include "sstream"
#include "string"
#include "TStyle.h"
#include "TCanvas.h"


using namespace std;

void ellipse( string inFile = "allSim.root", int cSpecies = 1 ){

	TCanvas * c = new TCanvas( "c", "c", 800, 800 );
	string outName = "rpEllipsePid.pdf";
	
	c->Print( (outName+"[").c_str() );

	TFile * f = new TFile( inFile.c_str(), "READ" );
	string rOutName = "rootEllipsePid.root";
	TFile * fOut = new TFile( rOutName.c_str(), "RECREATE" );

	vector<double>effVsP;
	vector<double>pureVsP;


	c->Divide( 2, 2 );
	for ( int i = 0; i < 70; i++ ){

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

		pidEllipse * pid = new pidEllipse( sum, cSpecies, 1.0, p0, p1, p2 );

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

		pid->yield( (TH2D*)sum, 1, 1);
		
		
		c->cd(2);
		sstr.str("");
		sstr << "eff_" << i;
		TH1D* eff = pid->efficiency( sstr.str(), 0.0, 5.0, 0.1 );

		sstr.str("");
		sstr << "pure_" << i;
		TH1D* pure = pid->purity( sstr.str(), 0.0, 5.0, 0.1 );

		gStyle->SetOptStat( 0 );
		eff->SetTitle( "Efficiecy (Blue), Purity (Red)" );
		eff->GetYaxis()->SetRangeUser(0, 1.05);
		eff->SetLineWidth( 2 );
		eff->Draw();
		pure->SetLineColor( kRed );
		pure->SetLineWidth( 2 );
		pure->Draw("same");

		effVsP.push_back( pid->efficiency() );
		pureVsP.push_back( pid->purity() );

		c->Print( outName.c_str());
	}

	int nBins = (3.7 - 0.2) / 0.05;
	TH1D * hEffVsP = new TH1D( "hEffVsP", "Efficiency Vs. P; P [GeV]", nBins, 0.2, 3.7 );
	for ( int i = 0; i < effVsP.size(); i++ ){
		hEffVsP->SetBinContent( i, effVsP[ i ] );
	}
	TH1D * hPureVsP = new TH1D( "hPureVsP", "Purity Vs. P; P [GeV]", nBins, 0.2, 3.7 );
	for ( int i = 0; i < pureVsP.size(); i++ ){
		hPureVsP->SetBinContent( i, pureVsP[ i ] );
	}

	c->Divide( 1 );
	c->cd( 0 );
	hEffVsP->GetYaxis()->SetRangeUser( 0.0, 1.05);
	hEffVsP->SetLineWidth( 2 );
	hEffVsP->Draw( "");
	hPureVsP->SetLineColor( kRed );
	hPureVsP->SetLineWidth( 2 );
	hPureVsP->Draw( "same" );

	c->Print( outName.c_str());


	c->Print( (outName+"]").c_str() );

	fOut->Write();


}