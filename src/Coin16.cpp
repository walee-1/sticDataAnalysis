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
//#include "EventDict.h"
#include "EventType.h"


#include "tdc_ch_values.h" // for DNL correction

using namespace std;

double physical_channel[16]={0}; //array for 16 physical channel histogram mapping
char* file;
char treename[16] = "dump";
unsigned int num_ch = 2;	//default number of channels (= required minimum)
char* outputName;
int CTW = 32; 			// Coincidence timing window given in fine counter bins (1 bin = 50.2 ps)
int start_channel;
int coin_counter;

TH1F *histo_FE = new TH1F("frame_entries","frame_entries",20,0,20); //histo for number of events per frame (FE= frame events)
			// name, name, num_bins, start, end)
TH1F *histo_coin = new TH1F("num_coin", "Number of coincident channels", 16, 0, 16);
TH1F *histo_coincident_channels = new TH1F("coin_ch", "channels involved in coincident events", 64, 0, 64);
TH1F *histo_cpf = new TH1F("cpf", "coincidences per frame", 10,0,10);
TH1F *histo_frames_skipped = new TH1F("frames_skipped","number of skipped frames between two events (i.e. frame doesn't exist)",50,0,50);


TH1F *histo_energy_coin = new TH1F("Energy_coincident_channels","Energy of coincident channels",50000,0,50000); //historgram for coincident channel's energies for coincident events only (so to get an energy spectra)
TH1F *histo_energy_coinlist1 = new TH1F("Energy_coincident_1st","Energy of first coincident channel",5000,0,5000);
TH1F *histo_energy_coinlist2 = new TH1F("Energy_coincident_2nd","Energy of second coincident channel",15000,0,10000);
//TH1F *histo_energy_ch2 = new TH1F("Energy_ch2","Energy of ch2", 2000,0,2000); //histogram for debugging purposes, shows the full energy of a specific channel (chosen in the code)

TH1F* histo_timediff_1_2 = new TH1F("Histogram", "Time difference between first and second coincident channel", 32, 0, 32*50.2);
TH1F* histo_timediff_2_3 = new TH1F("timeDiff_2_3", "Time difference between second and third coincident channel", 32, 0, 32*50.2);
TH1F* histo_timediff_3_4 = new TH1F("timeDiff_3_4", "Time difference between third and fourth coincident channel", 32, 0, 32*50.2);
TH1F* histo_timediff_4_5 = new TH1F("timeDiff_4_5", "Time difference between fourth and fifth coincident channel", 32, 0, 32*50.2);
TH1F* histo_timediff_1_5 = new TH1F("timeDiff_1_5", "Time difference between fourth and fifth coincident channel", 32, 0, 32*50.2); //last and first histo time diff

TH2F *sipm_3d = new TH2F("Sipm_Graphical", "Graphical representation of the Sipm", 4,0,4,4,0,4); //for lego plot of the sipm

TH2F *sipm_text = new TH2F("2d_sipm","2d_sipm",4,0,4,4,0,4);

TCanvas *c1=new TCanvas("sipm","sipm",600,600); //canvas for 2d drawing of sipm
TCanvas *c3=new TCanvas("sipm3d","3dsipm",3000,3000); //canvas for 3d drawing of sipm (simply storing the histogram doesn't work properly)
TCanvas *c2=new TCanvas("Coincident_Energies"); //canvas for the coincident energies.


void sortList(list<stic3_data_t> *eventList);
void getT_Diff_First_Last(list<stic3_data_t> *eventList);
void coinSearch(list<stic3_data_t> *eventList);
void boxdraw(double x1,double y1,double x2,double y2, int index, double max);//--test--might change it with a 2d histogram contour or colored 2d histogram
void sipmdraw(double max);
void sipmdraw3d(double max); //function for 3d sipm drawing
void channel_mapping();
double find_max_array(int length);

//***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***
int main(int argc, char *argv[]){

	if (argc<5){
		printf("use: %s file outputName starting_channel_number_from_STiC number_of_coincidences\n",argv[0]);
		printf("(file extention for outputName will be added automatically) \n\n");
		return -1; //simple display or passed parameters
	}
	coin_counter=atoi(argv[4])-1;
	file = argv[1]; //input file passed from argument
	outputName = argv[2]; //output file passed from argument
	start_channel=atoi(argv[3]); //starting channel number from argument

	// get the tree from the root file
	TFile f(file);
	TTree *tree=(TTree*)f.Get(treename);
	if (tree==NULL) {
		printf("No tree is found!\n");
		return -1;
	}


	//*************************************************************************
	// set branch address to address of declared instance of class stic3_data_t
	//*************************************************************************
	stic3_data_t* event=NULL;
	tree->SetBranchAddress("br",&event);

// 	Create List of stic3_data_t types. i.e. create list of events to be filled later in this programme

	//**********************************************************
	// Get entries of tree 
	// make analysis with all entries that have same frame number
	// invoke corresponting sorting and coincidence search functions
	//  fill histogramms there
	//**********************************************************
	int num_entries = tree->GetEntries();
	cout << "Entries: " << num_entries<< endl;
	int i=0;
	cout << i << endl;
	cout<< "Percentage of events processed ... " << endl;


	tree->GetEntry(i);
        int fr_num = event->frame_number;
	i=0;
	cout <<"fr_num at start: " << fr_num << endl;

	std::list<stic3_data_t> temp_list; //TEMPORARY LIST FOR COINCIDENCE SEARCH
        temp_list.clear();

	while( i < (num_entries)) {

                tree->GetEntry(i);
		//if(event->channel == 41 && event->energy>10)
//{
		//	histo_energy_ch2->Fill(event->energy);
//}
		if(event -> frame_number == fr_num){			// all events from same frame are put in the list 
			temp_list.push_back(*event);			// regardsless of its channel number
			i++;
		}else{
			histo_FE->Fill(temp_list.size());
			histo_frames_skipped->Fill(((event->frame_number)-fr_num)-1);
			sortList(&temp_list);				// sort the elements in the list by their time stamp (CALL BY REFERNCE!)
//			getT_Diff_First_Last(&temp_list);
			coinSearch(&temp_list);

			fr_num=event->frame_number;		// set new frame number
		        temp_list.clear();
		}
	        if (i % 100000 == 0) cout<< "\t" << (double)i/num_entries*100 << "\t" << "%" << endl;
	}
	cout << "\t" << 100 << "\t" << "%" << endl;

	//***********************************************************
	// write the histograms to file
	//***********************************************************

	cout << "Writing histos...." << endl;

	char out[128];
	sprintf(out, "%s_results.root", outputName);
	TFile* fout=new TFile(out,"RECREATE");


	histo_FE->Write();
	histo_coin->Write();
	histo_coincident_channels->Write();
	histo_cpf->Write();
	histo_timediff_1_2->Write();
	histo_timediff_2_3->Write();
	histo_timediff_3_4->Write();
	histo_timediff_4_5->Write();
	//histo_energy_coin->Write();

	//histo_energy_ch2->Write();
	c2->cd();
	gStyle->SetOptStat(111111);//--test--not doing it for the histogram only for the canvas, so might change this altogether
	
	histo_energy_coin->Draw();
	
	char out1[128]; //--test--if you want to save an output pdf of coincident energy graph or need a canvas for it.
	sprintf(out1, "%s_Coincident_E.pdf", outputName);
	c2->SaveAs(out1);
	c2->Write();

	histo_frames_skipped->Write();
	

//************** for drawing SiPM***********************

	cout<<"Assigning channels"<<endl;
	channel_mapping();
	double max=find_max_array(16); //just to find the maximum value for normalization which is used in coloring the graphical drawing
	cout << "max:\t" << max << endl;
	cout<<"Drawing SIPM"<<endl;
	sipmdraw3d(max);
	c1->cd();
	sipm_3d->GetXaxis()->SetTitle("Bottom of SiPM");
	sipm_3d->GetXaxis()->SetNdivisions(4);
	sipm_3d->GetYaxis()->SetNdivisions(4);
	sipm_3d->SetStats(0);
	c1->SetRightMargin(0.18);
	sipm_3d->Draw("COLZ");
	sipm_text->Draw("SAME TEXT");
	c1->Write();
	c3->cd();
	sipm_3d->SetStats(0);
	c3->SetBottomMargin(0.08);
	sipm_3d->Draw("LEGO2 PLC");
	
	c3->Write();
//******************************************************
	fout->Close();
cout << "...finished! " << endl;
	return 0;
}

//************************************************************************************************************************************************************
//**************END OF MAIN****************END OF MAIN*********************END OF MAIN************************************END OF MAIN*************************
//************************************************************************************************************************************************************

void sortList(list<stic3_data_t> *eventList){

//	THIS JUST PRINTS THE LIST BEFORE IT WAS SORTED
/*
	list<stic3_data_t>::iterator test_it;
	test_it=eventList->begin();
	do{
		cout << test_it->time << "\t" << test_it->T_CC << "\t" << test_it->T_fine <<  endl;
		++test_it;
	}while (test_it != eventList->end());
	cout << " -------------------------------------- " << endl;
*/

// SORT THE LIST WITH BUBBLE SORT
	list<stic3_data_t>::iterator it;
	for(it = eventList->end(); it != eventList->begin(); --it){
//		cout << "for1" << endl;
		list<stic3_data_t>::iterator it2= eventList->begin();
		while (it2 != it) {
//			cout << "it2 invoked" << endl;

			list<stic3_data_t>::iterator next_it = it2;
			++next_it;
			if(next_it !=it && it2->time > next_it->time){
				std::swap(*it2, *next_it);
				++it2;
			}else{
			++it2;
			}

		}
	}

/*		EXAMPLE WITH INTEGERS:
		for(int n = 4; n>1; n--){
			for(int i = 0; i<n-1; i++){
				if(times[i]>times[i+1]){
					int a = times[i];
					times[i]=times[i+1];
					times[i+1]=a;
					int b = cn[i];
					cn[i]=cn[i+1];
					cn[i+1]=b;
				}
			}
		}
*/
}

void getT_Diff_First_Last(list<stic3_data_t> *eventList){			// prints out the time difference between first and last event in one frame

	list<stic3_data_t>::iterator it = eventList->begin();
	unsigned int time_1 = it->time;
	unsigned int time_last =0;

	for(it = ++eventList->begin(); it != eventList->end(); ++it){

		time_last = it->time;

	}
	cout << "time 1:\t" << time_1 << "\ttime last:\t" << time_last << "\ttime diff:\t" << (time_last-time_1) << endl;
}


void coinSearch(list<stic3_data_t> *eventList){

	list<stic3_data_t>::iterator it;
	int coinChannels[64];
	unsigned int ct_coincidences =0;
	unsigned int timeStamps[16];	// stores the time stamps of each single event involved in a coincidence 
	int cpf = 0; 	// coincidences per frame, i.e. how many "bunches" of events were there 
			// that belonged to ONE coincident LED event (or Cherenkov Light Cone respectively)

	for(it = eventList->begin(); it != eventList->end(); ++it){
		unsigned int ct_coincidences =0;	//counter for coincidences
		double dummy_energy=it->energy;

		coinChannels[ct_coincidences] = it->channel;
		timeStamps[ct_coincidences] = it->time;	// saves the time stamp of each event in this coincidence
		ct_coincidences++;
		list<stic3_data_t>::iterator it2=it;
		++it2;
		do {
			int time_diff = (it->time)-(it2->time); 
			if ( time_diff < CTW && time_diff > -CTW && it2->energy>70 ) {		// CTW = coincidence timing window
				coinChannels[ct_coincidences] = it2->channel;
				timeStamps[ct_coincidences] = it2->time;
				dummy_energy=dummy_energy+it2->energy;
				/*histo_energy_coinlist1->Fill(it->energy);
				histo_energy_coinlist1->SetLineColor(2);//red
				histo_energy_coinlist2->Fill(it2->energy);
				histo_energy_coinlist2->SetLineColor(3);//green*/
				ct_coincidences++;
				if(it2->T_badhit !=33){				// Asking if it's a bad hit...
					it2 = eventList->erase(it2);
				}else{
					cout << "BAD HIT = 33! Corrupted event... Returning to main function..." << endl;
					cout << "Channel where bad Hit Encountered: "<<it2->channel <<endl;
					cout << "Energy where bad Hit Encountered: "<<it2->energy <<endl;
					cout << "time of bad hit "<<it2->time << endl;
					return;
				}
			}else{
				it2 = --it2;		// go one position back
				it= it2;		// take this position for iterator 1
				break;			// this position will be incremented by the for loop
			}				// ending up at the first position not to be a concidence
		}while ( it2 != eventList->end() );

		if ( ct_coincidences > coin_counter) {		// Define a coincident event to have 3 or more events!
			histo_coin->Fill(ct_coincidences);
			cpf++;				// Increment counter for coincidences per frame
			histo_timediff_1_2->Fill(50.2*(timeStamps[1]-timeStamps[0]));	// Fill histo with first time difference
			histo_timediff_2_3->Fill(50.2*(timeStamps[2]-timeStamps[1]));
			if (ct_coincidences > 3) histo_timediff_3_4->Fill(50.2*(timeStamps[3]-timeStamps[2]));
			if (ct_coincidences > 4) {
				histo_timediff_4_5->Fill(50.2*(timeStamps[4]-timeStamps[3]));
				histo_timediff_1_5->Fill(50.2*(timeStamps[4]-timeStamps[0]));	//fill histo with first and last time stamp
			}	
			for (int i = 0; i < ct_coincidences; i++){
				histo_coincident_channels->Fill(coinChannels[i]);
			}
		}
		if (ct_coincidences > coin_counter)
		{
			histo_energy_coin->Fill(dummy_energy);
		}

	}// END OF FOR LOOP

	histo_cpf->Fill(cpf);	//coincidences per frame


} // END OF coinSearch(..){}


/*


	//**********************************************************
	// Here is a demo on how to do a DNL correction.
	//**********************************************************
	cout << "starting dnl (not included in results yet!!) " << endl;

	TH1F *hist_dnl_tfine = new TH1F("T_fine after DNL","T_fine after DNL",128,0,128); // DNL correction mapping 32 bins to 128 bins
	char hname[64];
	sprintf(hname, "channel_%d",ch);
	tdc_ch_values *dnl_correction = new tdc_ch_values(0,31,hname);

	// fill all the fine counter values for code density test
	for(int i=0; i<tree->GetEntries();i++){
		tree->GetEntry(i);
		if(event->channel == ch) {
			dnl_correction->fill(event->T_fine);
			dnl_correction->fill(event->E_fine);
		}
	}

	// calculate the mapping for bin dithering
	dnl_correction->update_histos();


	//***********************************************************
	// Use the corrected time information below for serious timing analysis.
	// Implement your own analysis in the loop below
	//***********************************************************
	unsigned int t_fine_dnl;
	double time_dnl;
	for(int i=0; i<tree->GetEntries();i++){
		tree->GetEntry(i);
		if(event->channel == ch) {
			t_fine_dnl = dnl_correction->get_dith_value(event->T_fine); // Get T_fine after DNL correction
			hist_dnl_tfine->Fill(t_fine_dnl);

			time_dnl = (event->T_CC *128 + t_fine_dnl)*(1.0/622.08e6/128);  // calculate the real time value with new T_fine value, The PLL clock frequency is 622.08e6 Hz
		//	printf("old time: %e s; new time: %e s \n", event->time*50.2e-12, time_dnl);
		//	fflush(stdout);

		}
	}
	hist_dnl_tfine->Write();
*/

void channel_mapping() //function for mapping stic channels to physical sipm channels
{
	double out=0;
	double temp=0;
	int sticch[64]={0};
	int sipmch[64]={0};
	int x,y;
	int count=0;
	ifstream file2("../src/stic2channel.txt");
	string line; 
 	while(getline(file2, line)) //reads data from file and stores it
   	{
      		if (!line.length())
         		continue;
      		sscanf(line.c_str(), "%d  %d", &x, &y);
		sticch[count]=x;
		sipmch[count]=y;
		count++;
   	}
	for (int i=start_channel; i<start_channel+16;i++)
	{	
		temp=histo_coincident_channels->GetXaxis()->FindBin(sticch[i]); //stores the stic channel number to a temp variable
		out=histo_coincident_channels->GetBinContent(temp); //the channel number's coincident hit values are then stored in another variable
		physical_channel[sipmch[i]]=out; //the out variable, containing the bin content of the stic channel number is finally stored in the corresponding physical sipm channel number
	}
	
}

double find_max_array(int length) //just finds the max value out of an array as well as prints the results while searching
{
	double max=0;
	for(int i=0; i<length;i++)
	{
		if(max<physical_channel[i])
		{
			max = physical_channel[i];
		}
		cout <<i << "\t" << physical_channel[i] << "\t" << max << endl;
	}
	return(max);

}
void sipmdraw3d(double max)
{
	int k=0;
	for(int y=0;y<4;y++)
	{
		for(int x=0;x<4;x++)
		{
			sipm_3d->Fill(x,y,physical_channel[k]);
			sipm_text->Fill(x,y,k+1);
			k++;
			
		}
	}
}
