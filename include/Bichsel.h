
#ifndef DEDX_BICHSEL_H
#define DEDX_BICHSEL_H

#include <string>
#include <math.h>
#include "TH1D.h"
#include "TFile.h"

using namespace std;

class Bichsel
{
protected:
	string tableFile;
	int method;
	TFile* table;

	TH1D *hP, *hPi, *hK;

	static constexpr double piMass = 0.1395702;
	static constexpr double kaonMass = 0.493667;
	static constexpr double protonMass = 0.9382721;

public:
	Bichsel( string tableFile = "dedxBichsel.root", int method = 0 ){

		this->tableFile = tableFile;
		this->method = method;
		table = new TFile( tableFile.c_str(), "READ" );

		getTables();

	}
	~Bichsel(){

	}

	double mean( double p, double m, int method = -1 ){

		if ( method >= 0 ){
			getTables();
		}

		// use the mass to determine the particle species to use
		// This is to maintain the same usage 
		TH1D * h = tableFor( m );

		if ( h ){
			
			int bin = h->GetXaxis()->FindBin( p );
			return ( h->GetBinContent( bin ) );

		} else 
			cout << " error no table found " << endl;

		return -999.999;
	}

	double getFromTable( string plc, double p );
	
	void getTables(){

		if ( 1 == method  ){
			hP = (TH1D*)table->Get( "t70P" );
			hK = (TH1D*)table->Get( "t70K" );
			hPi = (TH1D*)table->Get( "t70Pi" );
		} else if ( 2 == method ){
			hP = (TH1D*)table->Get( "polP" );
			hK = (TH1D*)table->Get( "polK" );
			hPi = (TH1D*)table->Get( "polPi" );
		} else {
			hP = (TH1D*)table->Get( "mpmP" );
			hK = (TH1D*)table->Get( "mpmK" );
			hPi = (TH1D*)table->Get( "mpmPi" );
		}
	}

	static constexpr double epsilon = 0.01;
	TH1D* tableFor( double mass ){
		//cout << " Mass : " << mass << endl;
		if ( abs(mass - piMass) < epsilon  ){
			//cout << "\tid: Pi" << endl; 
			return hPi;
		}
		else if ( abs(mass - kaonMass) < epsilon  ){
			//cout << "\tid: K" << endl;
			return hK;
		}
		else if ( abs(mass - protonMass) < epsilon  ){
			//cout << "\tid: P" << endl;
			return hP;
		}
		return hP;
	}

	static const int MPM = 0;
	static const int T70 = 1;
	static const int POL = 2;
	
};

#endif