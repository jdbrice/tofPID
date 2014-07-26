#ifndef PID_FITTER_H
#define PID_FITTER_H

#include "allroot.h"
#include "histoBook.h"
#include "constants.h"
#include "xmlConfig.h"
#include "utils.h"
#include "reporter.h"
#include <vector>
// clock_t, clock, CLOCKS_PER_SEC 
#include <time.h>  

class pidFitter
{
public:
	pidFitter( xmlConfig * config  );
	~pidFitter();

protected:
	// the canvas used to draw report hidtos
	reporter* report;

	// the histobook that stores all of our pidHistogramMaker histograms
	histoBook *book;

	// config file
	xmlConfig* config;

	clock_t startTime;

	map< string, reporter * > pReport;

	/*
	*	Utility functions that should be moved soon
	*/ 
	void startTimer( ) { startTime = clock(); }
	double elapsed( ) { return ( (clock() - startTime) / (double)CLOCKS_PER_SEC ); }
	

	string sName( string pType, int charge );

	void runFit();
	void runDkl( TH2D*, reporter * rp, uint nS, uint nIt );


};



#endif