


void fitLUT(){



	TFile * f = new TFile( "lut.root", "READ" );

	TH1D* xMeanK = f->Get( "xMeanPi" );
	TH1D* yMeanK = f->Get( "yMeanK" );
	TH1D* xSigmaK = f->Get( "xSigmaK" );
	TH1D* ySigmaK = f->Get( "ySigmaK" );


	TF1 * meanFit = new TF1( "fitMean", fitMean, 0.2, 1.1, 3);
	xMeanK->Fit( meanFit, "R" );

	xMeanK->Draw();


}


double fitMean( Double_t *x, Double_t *par ){
	float xx = x[0];
	return (par[0] + par[1] * TMath::Power(xx, par[2] ));
}