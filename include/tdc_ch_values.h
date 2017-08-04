#ifndef __TDC_CH_VALUES_H__
#define __TDC_CH_VALUES_H__

#include<iostream>
#include<fstream>
#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
#include<TSpectrum.h>
#include<TH1F.h>
#include<TCanvas.h>
#include<TPad.h>
#include<TFile.h>
#include<time.h>
#include<math.h>


const double PLL_REFCLK_FREQ = 622.08e6;

const unsigned int MAPPING_N_BINS = 128;
//const unsigned int MAPPING_N_BINS = 64;
//const unsigned int MAPPING_N_BINS = 32;

struct dithermap{
	double prob[50];
	int assign[50];
};


class tdc_ch_values{

   public:

	//FUNCTIONS TO RETURN THE CREATED HISTOGRAMS

	TH1F* get_cdt() {return this->cdt;};
	TH1F* get_DNL() {return this->DNL;};
	TH1F* get_INL() {return this->INL;};
	TH1F* get_realtimedist() {return this->realtimedist;};
	TH1F* get_mps() {return this->mps;};

	//FUNCTIONS TO FILL AND CREATE THE HISTOGRAMS

	tdc_ch_values(int min,int max, const char* hname); //CONSTRUCTOR: ALLOCATE ALL THE NEEDED HISTOGRAMS

	void fill(unsigned int tdc_value);	
	double get_ps_value(unsigned int tdc_value);
	unsigned int get_dith_value(unsigned int tdc_value);
	int update_histos();

   private:
	void create_map(TH1F* realtime); //CREATES THE MAP FOR PSEUDORANDOM BIN DITHERING
	unsigned int bin_dice(int bin); //RETURNS THE ASSIGNED VALUE OF THE GIVEN BIN
	
	//START AND STOP OF THE BINS CONTAINING HITS
	int start_bin;
	int stop_bin;

	struct dithermap *interpol;
	TH1F *cdt;
	TH1F *DNL;
	TH1F *INL;

	TH1F* realtimedist;
	TH1F* mps;

};
#endif
