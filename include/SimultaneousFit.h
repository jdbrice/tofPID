#ifndef SIMULTANEOUSFIT
#define SIMULTANEOUSFIT

#include "XmlConfig.h"
#include "tofGenerator.h"
#include "Logger.h"
#include "LoggerConfig.h"

using namespace jdb;

#include "TH1D.h"


class SimultaneousFit
{
protected:

	// describe the state
	TH1D *tofAll, *tofS0, *tofS1, *tofS2;
	TH1D *dedxAll, *dedxS0, *dedxS1, *dedxS2;

	tofGenerator *tGen;

	double avgP;

	XmlConfig * config;
	string path;

	Logger * lg;

public:
	SimultaneousFit( 	TH1D* tAll, TH1D* tS0, TH1D* tS1, TH1D* tS2,
						TH1D* dAll, TH1D* dS0, TH1D* dS1, TH1D* dS2,
						double p, XmlConfig * config = NULL, string nodepath ="" 
		);
	~SimultaneousFit();

	void fitTofAll( );
	
};



#endif