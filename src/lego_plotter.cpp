#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TH2F.h"
#include "TBox.h"
#include "TObject.h"
#include "TColor.h"
#include "TStyle.h"


#include <iostream>
#include <stdio.h>
#include <list>



void lego_plotter()
{
	TH2F *sipm_3d = new TH2F("Sipm_Graphical", "Graphical representation of the Sipm", 4,0,4,4,0,4); //for lego plot of the sipm
	TFile* fout=new TFile("results.root","RECREATE");
	TH2F *sipm_text = new TH2F("2d_sipm","2d_sipm",4,0,4,4,0,4);
	double avg_ph;
	TCanvas *c1=new TCanvas("sipm","sipm",600,600); //canvas for 2d drawing of sipm
	TCanvas *c3=new TCanvas("sipm3d","3dsipm",3000,3000); //canvas for 3d drawing of sipm (simply storing the histogram doesn't work properly)
	sipm_3d->GetXaxis()->SetNdivisions(4);
	sipm_3d->GetYaxis()->SetNdivisions(4);
	sipm_3d->SetStats(0);
	c3->cd();
	int k=0;
	for(int y=0;y<4;y++)
	{
		for(int x=0;x<4;x++)
		{
			cout<<"Enter value for channel "<<k;
			cin>>avg_ph;
			sipm_3d->Fill(x,y,avg_ph);
			sipm_text->Fill(x,y,k+1);
			k++;	
		}
	}

	sipm_3d->SetStats(0);
	c3->SetBottomMargin(0.08);
	sipm_3d->Draw("LEGO2 PLC");
	simp_3d->Write();
	
}
