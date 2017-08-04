#include "tdc_ch_values.h"

//#define DEBUG_DNL
//#define DEBUG_DITH

extern unsigned const int MAPPING_N_BINS; // set the MAPPING to N Bins in file: include/tdc_ch_values.h
extern const double PLL_REFCLK_FREQ; // set the PLL_REFCLK_FREQ in the file include/tdc_ch_values.d



tdc_ch_values::tdc_ch_values(int min,int max,const char* hname) //Constructor: allocate all the needed histograms
{

	char name[31];

	sprintf(name,"cdt_%s",hname);
	cdt = new TH1F(name,"TDC Channel CDT",(max-min)+1,min,max+1);

	sprintf(name,"dnl_%s",hname);
	DNL = new TH1F(name,"DNL Evaluation without remapping",(max-min)+1,min,max+1);

	sprintf(name,"inl_%s",hname);
	INL = new TH1F(name,"INL Evaluation",(max-min)+1,min,max+1);

	sprintf(name,"real_%s",hname);
	realtimedist= new TH1F(name,"real bin delays calculated from CDT",(max-min)+1,min,max+1);

	sprintf(name,"mps_%s",hname);
	mps = new TH1F(name,"probability distribution of binsizes",400,0,200e-12);

	//allocate the map for bin dithering
	interpol = new struct dithermap[32];

}

void tdc_ch_values::create_map(TH1F* realtime)
{
#ifdef DEBUG_DNL
	FILE* fout;
	fout =fopen("create_map.log","w");
#endif

	memset(interpol,0,sizeof(interpol));
	int i;

	for(i=0; i<32;i++)
	{
		memset(interpol[i].prob,0,sizeof(interpol[i].prob));
		memset(interpol[i].assign,0,sizeof(interpol[i].assign));
	}

	int nbins,j;
	double tlow,thigh,tsize;

	int binlow, binhigh;
	double blr, bhr;

	double binsize;

	//binsize=(1.0/(PLL_REFCLK_FREQ)/128.0);	//Binsize of the bins we map to, 128 bins
	//binsize=(1.0/(PLL_REFCLK_FREQ)/64.0);	//Binsize of the bins we map to, 64 bins
	//binsize=(1.0/(PLL_REFCLK_FREQ)/32.0);	//Binsize of the bins we map to, 32 bins
	binsize=(1.0/(PLL_REFCLK_FREQ)/MAPPING_N_BINS);	//Binsize of the bins we map to
#ifdef DEBUG_DITH
	std::cout << "New binsize we map to: " << binsize << std::endl;
#endif
//	binsize=1.0;
	tlow=0;	//lower and higher timeborder of the bin
	thigh=0;

	for (i=0;i<realtime->GetNbinsX();i++)
	{
		tlow=thigh;
		thigh += realtime->GetBinContent(i+1);	//Get the Content of bin+1 because STUPID ROOT USES BIN 0 AS UNDERFLOW!!!! 
		tsize = realtime->GetBinContent(i+1); 
		binlow = (int) (tlow/binsize);
		binhigh = (int) (thigh/binsize);


		blr = binsize - (tlow-binlow*binsize); //much better calculation of the remainders
		bhr = thigh - binhigh*binsize;


#ifdef DEBUG_DNL
		fprintf(fout,"creating dithermap: Real bin %u : size: %E tlow: %E thigh: %E binlow: %u binhigh: %u \n", i,tsize,tlow,thigh, binlow, binhigh);
#endif


//		if(i==237)
//		{
//			printf("\n Additional information: Calculated binlowremain: %f, binhighremain: %f\n\n",blr,bhr);
//		}


		if(binlow == binhigh)	//THE BIN IS COMPLETLY WITHIN A BIN IT GETS MAPPED TO SO NO NEED TO CALCULATE THE PROBABILITIES
		{
			interpol[i].prob[0] = 1.0;
			interpol[i].assign[0] = binlow;
		}
		else
		{
			nbins = binhigh-binlow-1;	//how many bins except the start and end bin are occupied?
			interpol[i].prob[0] = blr/tsize;	//assign the probabilities of the low and high bin	
			interpol[i].assign[0] = binlow;
//			if(i==231) printf("\n\n%f %f %u %u ",interpol[i].prob[0],interpol[i].prob[1],interpol[i].assign[0],interpol[i].assign[1]);
//			printf(" nbins: %u \n",nbins);
			if(nbins != 0)
			{
				for(j=0;j<nbins;j++)		//assign the probabilites for the rest of the occupied bins
				{
					interpol[i].prob[1+j] = blr/tsize+(j+1)*(tsize-bhr-blr)/(tsize*nbins); //Take care with the probabilities
					interpol[i].assign[1+j] = binlow+1+j;
//					if(i==231) printf("%f %u \n",interpol[i].prob[1+j],interpol[i].assign[1+j]);
				}
			}
#ifdef DEBUG_DNL
			else {
			  printf("dnl create_map:nbins==0; this bin is mapping to two new bins.\n");
			}
#endif

			interpol[i].prob[nbins+1] = interpol[i].prob[nbins]+bhr/tsize;
			interpol[i].assign[nbins+1] = binhigh;


		}
#ifdef DEBUG_DNL
		for (int ii = 0; ii <= (binhigh - binlow); ii++) {
			if (ii == 0) {
				fprintf(fout,"Mapping current bin to bin: %u probability: %e accumulated probablily: %e \n", binlow+ii, interpol[i].prob[ii], interpol[i].prob[ii]);
			}
			else {
				fprintf(fout,"Mapping current bin to bin: %u probability: %e accumulated probablily: %e \n", binlow+ii, interpol[i].prob[ii] - interpol[i].prob[ii-1], interpol[i].prob[ii]);
			}
		}
#endif

	}
#ifdef DEBUG_DNL
	fclose(fout);
#endif
}

double tdc_ch_values::get_ps_value(unsigned int tdc_value)
{
	double ps1=0;
	int i;
	for(i=start_bin+1;i<=(0x1F&tdc_value)+1;i++)
	{
		ps1 += realtimedist->GetBinContent(i);
	}
#ifdef DEBUG_DNL
//	std::cout << "TDC Value: " << tdc_value << " CC_Value: " << (tdc_value >> 5)/622.02e6 << " picosecond addition: " << ps1 << std::endl;
#endif
	ps1 = (tdc_value >> 5)/PLL_REFCLK_FREQ + ps1;
	return ps1;
}

//USE RANDOM NUMBER GENERATOR TO GET THE DITHERED BIN NUMBER
unsigned int tdc_ch_values::bin_dice(int bin)
{
	double rnumber;

	rnumber = (double)(rand() % 10000)/10000.0;
	FILE *ranlog;
	ranlog = fopen("random.log","a");
	fprintf(ranlog,"%f\n",rnumber);
	fclose(ranlog);
	float last_prob;
	last_prob = 0.0;
	for(int i = 0; i<50;i++)
	{
		if(last_prob <= rnumber && interpol[bin].prob[i] > rnumber) 
		{
			return interpol[bin].assign[i];
		}
		last_prob = interpol[bin].prob[i];
	}
	printf("bin_dice:loop over 50 is not enough!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
}


unsigned int tdc_ch_values::get_dith_value(unsigned int tdc_value)
{
	unsigned int temp;
	temp = 0x1F & tdc_value;
	temp = bin_dice(temp);
	//temp = ((tdc_value - (0x1F&tdc_value))<<2) + temp; //128 bins mapped
	//temp = ((tdc_value - (0x1F&tdc_value))<<1) + temp; //64 bins mapped
	//temp = ((tdc_value - (0x1F&tdc_value))<<0) + temp; //32 bins mapped
	//temp = ((tdc_value - (0x1F&tdc_value))<<(unsigned int)log2(MAPPING_N_BINS/32)) + temp; // MAPPING_N_BINS bins mapped
	temp = ((tdc_value - (0x1F&tdc_value))<<(unsigned int)(log2(MAPPING_N_BINS/32))) + temp; // MAPPING_N_BINS bins mapped
	//std::cout<<"bin shift: "<<(unsigned int)(log2(MAPPING_N_BINS/32))<<" !!!!!!!!!!!!!!!!!!!!!\n";
	return temp;

}

int tdc_ch_values::update_histos()
{

	int i,bin_number;
	double sum;
	start_bin = 64;
	stop_bin = 0;

	for(i=1; i <= cdt->GetNbinsX();i++)
	{
		bin_number = i - 1;		//~₢!"☠! ROOT -- Maybe I should stick to the bin number of ROOT but...
		if(cdt->GetBinContent(i)>2 && bin_number < start_bin) start_bin = bin_number;
		if(cdt->GetBinContent(i)>2 && bin_number > stop_bin) stop_bin = bin_number;
	}

#ifdef DEBUG_DITH
	std::cout << "Starting bin: " << start_bin << " Stop Bin: " << stop_bin << std::endl;
#endif


	float average_size = 0.0;
	unsigned int event_number;
	event_number = cdt->GetEntries();
	if(event_number < 100) return -1;
	average_size = (double)event_number/(double)(stop_bin-start_bin+1);

	//CALCULATE THE DNL AND THE REALTIME DISTRIBUTION OF THE BINS
	for(i=1; i <= cdt->GetNbinsX();i++)
	{
		DNL->SetBinContent(i, ((cdt->GetBinContent(i) - average_size)/average_size));
		realtimedist->SetBinContent(i,(cdt->GetBinContent(i)*((1.0/PLL_REFCLK_FREQ)/(double)event_number)));
	}

	//CALCULATE THE BIN WIDTH DISTRIBUTION
	for (i = 1; i <= realtimedist->GetNbinsX();i++)
	{
		if(realtimedist->GetBinContent(i) != 0) mps->Fill(realtimedist->GetBinContent(i));
	}

	//CREATE THE BIN DITHERING MAP
	create_map(realtimedist);


	//CALCULATE THE INL -- ONLY NEEDED FOR VECTOR CORRECTION 	
	sum=0;
	for(i=start_bin+1;i<=stop_bin+1;i++)	//First calculate the INL without remapping. Remapping should be the latest stage!!
	{
		sum = DNL->GetBinContent(i) + sum;
		INL->SetBinContent(i, sum);
	}


}


void tdc_ch_values::fill(unsigned int tdc_value)
{
	unsigned int time;
	time = 0x1F&tdc_value;
	cdt->Fill(time);
}
