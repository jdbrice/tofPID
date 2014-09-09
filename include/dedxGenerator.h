

/**
 * 	Original Author: Evan Sangaline
 *  Author: Daniel Brandenburg
 *  Date: September 3rd, 2014
 *  Description:
 *  Generates random dedx values as would be observed by the STAR TPC detector. 
 *  Attempts to add realistic detector effects.
 */

#ifndef DEDX_GENERATOR_H
#define DEDX_GENERATOR_H

#include "TROOT.h"
#include "TRandom3.h"
#include "TMath.h"

/*
	For landau_qauntile
 */
#include "Math/Math.h"
#include "Math/QuantFuncMathCore.h"

class dedxGenerator
{
protected:

	/**
	 * Parameters for the 8th order fit for dedx
	 */
	static const double pars[];
	/**
	 * Random number generator with periodicity 10**600
	 */
	TRandom3 * rGen;

public:
	dedxGenerator(){

		// setup the random number generator
		// Use the unique seed provided by ROOT
		rGen = new TRandom3( 0 );

	}
	~dedxGenerator(){

	}

	/**
	 * Calculates the expected mean dedx value ( in [KeV/cm] )
	 * using an 8th order polynomial parametrization
	 * @param  p    Momentum [KeV/c]
	 * @param  mass Mass
	 * @return      Expected dedx value for given p, m
	 */
	double mean( double p, double mass ){
		
    
    	const double bg = TMath::Log10( p / mass );
    	double rv = pars[0];
    	double xpow = 1;
	    for(int i = 1; i < 9; i++) {
	        xpow *= bg;
	        rv += xpow*pars[i];
	    }
	    return TMath::Power( 10, rv );
	}

	/**
	 * Generates a single randomized dedx value
	 * @param  p    Momentum [GeV/c]
	 * @param  mass Mass [GeV/c^2]
	 * @return      Random dedx value for the given p, m
	 */
	double random( double p, double mass ){
		double noise = 0;
	    double mu = mean( p, mass );
	    double sigma = mu * 0.15;
	    mu *= 0.9;

	    static double dedx_vals[50];
	    static int order[50];

	    const int dedx_points = rGen->Integer(30) + 15;

	    // Generate each point
	    for(int i = 0; i < dedx_points; i++) {
	        do {
	            
	            dedx_vals[i] = ROOT::Math::landau_quantile(rGen->Uniform(0, 1))*sigma + mu;
	            dedx_vals[i] += rGen->Gaus(0, noise*dedx_vals[i]);

	        } while ( dedx_vals[i] < 0 || dedx_vals[i] > 1e10 );
	    }
	    
	    // Sort the points into order
	    TMath::Sort(dedx_points, dedx_vals, order, false);

	    const int max_point = int(double(dedx_points)*0.7);

	    double sum = 0;
	    for(int i = 0; i < max_point; i++) {
	        sum += dedx_vals[order[i]];
	    }

	    sum /= double(max_point);

	    return sum;
	}
	
};

/**
 * An 8th order polynomial fit to b70
 */
const double dedxGenerator::pars[] = {-2.38469, -0.968702, 1.08892, -0.0461616, -0.449233, 0.15573, 0.0512726, -0.0348713, 0.00484003};


#endif