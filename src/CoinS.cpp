#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"

#include <iostream>
#include <stdio.h>
#include <list>

//#include "EventDict.h"
#include "EventType.h"


#include "tdc_ch_values.h" // for DNL correction

using namespace std;

char* file;
char treename[16] = "dump";
unsigned int ch = 0;
unsigned int ch2 = 0;
unsigned int ch3 = 0;
unsigned int ch4 = 0;
unsigned int num_ch = 2;	//default number of channels (= required minimum)
char* outputName;

int ch1_ct = 0;
int ch2_Ct = 0;
int ch3_ct = 0;
int ch4_ct = 0;

std::list<stic3_data_t> temp_list; //TEMPORARY LIST FOR COINCIDENCE SEARCH

TH1F *histo_time_diff = new TH1F("T_diff","T_diff",200,-100,100);
TH1F *histo_FE = new TH1F("frame_entries","frame_entries",20,0,20); //histo for number of events per frame (FE= frame events)
			// name, name, num_bins, start, end)
TH1F *histo_coin = new TH1F("coin", "Number of coincident channels", 8, 0, 8);

TH1F *t_diff_histo_1 = new TH1F("t_diff_1", "time difference between 1st and 2nd coincident channel", 50,0,50);
TH1F *t_diff_histo_2 = new TH1F("t_diff_2", "time difference between 2nd and 3rd coincident channel", 50,0,50);
TH1F *t_diff_histo_3 = new TH1F("t_diff_3", "time difference between 3rd and 4th coincident channel", 50,0,50);
TH1F *t_diff_histo_4 = new TH1F("t_diff_4", "time difference between first and last coincident channel", 50,0,50);
TH1F *first_coincident_channel = new TH1F("first_coin", "first of four coincident channels", 20, 15, 35);

void time_diff(list<stic3_data_t> *eventList, unsigned int one, unsigned int two);
void coinSearch(list<stic3_data_t> *eventList);

//***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***MAIN***
int main(int argc, char *argv[]){

	if (argc<5){
		printf("use: %s file outputName ch_1 ch_2 (ch_3) (ch_4) \n",argv[0]);
		printf("(file extention for outputName will be added automatically) \n\n");
		printf("Please note: \nThe ordering of the channels is important: Histogram for time difference will only be created for the first two channels! \n");
		return -1;
	}
	file = argv[1];
	outputName = argv[2];
	ch = atoi(argv[3]);
	ch2 = atoi(argv[4]);
	if (argc == 6) {
		ch3 = atoi(argv[5]);
		num_ch = 3;
	}
	if (argc == 7) {
		ch4 = atoi(argv[6]);
		num_ch = 4;
	}


	
/*
	int num_ch;
	cout << "Enter number of channels: ..." << endl;
	cin >> num_ch;
	cout << "You've entered the following number of channels: " << num_ch << endl;
*/
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
	// Here is a demo of how to loop over the tree and to fill a histogram
	// Serious timing analysis should not be done here
	// Do timing analysis in the loop with DNL correction below
	//**********************************************************
	int num_entries = tree->GetEntries();
	cout << "Entries: " << num_entries<< endl;
	int i=0;
	cout<< "Percentage of events processed ... " << endl;

        tree->GetEntry(i);
        int fr_num = event->frame_number;
        temp_list.clear();

	while( i < num_entries) {						// change here: Number of entries to be processed! 
                tree->GetEntry(i);	
		if(event -> frame_number == fr_num){			// all events from same frame are put in the list 
			temp_list.push_back(*event);			// regardsless of its channel number
			i++;	
		}else{	
			histo_FE->Fill(temp_list.size());
			time_diff(&temp_list, ch, ch2);
			coinSearch(&temp_list);
//			cout << "size: " << eventList->size()<< endl;
			fr_num=event->frame_number;		// set new frame number
		        temp_list.clear();
		}
	        if (i % 10000 == 0) cout<< "\t" << (double)i/num_entries*100 << "\t" << "%" << endl;
	}
	cout << "\t" << 100 << "\t" << "%" << endl;
	
	TH1F *hist = new TH1F("T_CC","T_CC",33000,0,33000);
//	TH1F *hist_name = new TH1F("name","head_line",num_bins,start,end);
	for(int i=0; i<tree->GetEntries();i++){
		tree->GetEntry(i);
		if(event->channel == ch) {
			hist->Fill(event->T_CC);
//			hist->Fill(event->T_fine);
		}
	}

	//***********************************************************
	// Create Energy histogramms for all channels:
	//***********************************************************

cout << "Creating Histos.... "<< endl;
	
        TH1F *hist_E[16];
	char temp_char[265];
	unsigned int num_channels_for_histo=0;
	for(int i = 3; i<argc; i++){	
		sprintf(temp_char, "E_%s",argv[i]);
		TH1F *temp_hist = new TH1F(temp_char, temp_char,1000,0,1000);
		for(int j=0; j<tree->GetEntries();j++){
            		           tree->GetEntry(j);
	     		   if( event->channel == atoi(argv[i]) ) {
		           temp_hist->Fill(event->energy);
		       	   }
       		}
		hist_E[num_channels_for_histo]=temp_hist;
		num_channels_for_histo++;
	}

	//***********************************************************
	// fit gauss to histogram and print results
	//***********************************************************
	
	histo_time_diff->Fit("gaus");
	TF1 *fit= histo_time_diff->GetFunction("gaus");

	double mean = fit->GetParameter(1);
	double FWHM = (fit->GetParameter(2))*2.355;
	cout << "\n" << "-------------------- RESULTS --------------------" << endl;
	cout << "Entries:" << "\t" << num_entries << endl;
	cout << "MEAN"<< "\t\t" << mean*50.2 <<"\t" << "ps" << endl;	
	cout << "FWHM" << "\t\t" << FWHM*50.2 << "\t" << "ps" << endl;

	//***********************************************************
	// write the histograms to file
	//***********************************************************

	cout << "Writing histos...." << endl;

	char out[128];
	sprintf(out, "%s_results.root", outputName);
	TFile* fout=new TFile(out,"RECREATE");
	hist->Write();
	t_diff_histo_1->Write();
	t_diff_histo_2->Write();
	t_diff_histo_3->Write();
	t_diff_histo_4->Write();
	first_coincident_channel->Write();

for(int n=0; n<num_ch;n++){
		hist_E[n]->Write();
}
	histo_time_diff->Write();
	histo_time_diff->GetXaxis()->SetTitle("Time in units of 50.2 ps");
	histo_FE->Write();
	histo_coin->Write();
	fout->Close();
cout << "...finished! " << endl;
	return 0;
}

//************************************************************************************************************************************************************
//**************END OF MAIN****************END OF MAIN*********************END OF MAIN************************************END OF MAIN*************************
//************************************************************************************************************************************************************

void time_diff(list<stic3_data_t> *eventList, unsigned int one, unsigned int two){
// cout << "Called function time_diff..." << endl;
//        cout << i << "\t" <<"frame " << event->frame_number << " ch: " << event->channel << " T_CC: " << event->T_CC << endl;
	list<stic3_data_t>::iterator it;

	for(it = eventList->begin(); it != eventList->end(); ++it){
		list<stic3_data_t>::iterator inner_it = it; 	// Define second iterator to compare the first one with
		++inner_it;					// Set iterator to next position

		while(inner_it != eventList->end() ){
                	if(((inner_it->channel == one) && (it->channel == two)) || ((inner_it->channel == two) && (it->channel == one))){
				int time_diff = (it->time)-(inner_it->time);        
				if( (time_diff < 32 ) && (time_diff > -32) ){	
//				if(it->T_CC == inner_it->T_CC){			// Coincidence window: One bin in the corse counter (= 1.6064 ns)
			//		cout << "Time-diff: " << time_diff << endl;
					
					histo_time_diff->Fill(time_diff);
				}
			}
		++inner_it;
		}
	}
}

void coinSearch(list<stic3_data_t> *eventList){
	list<stic3_data_t>::iterator it;
	for(it = eventList->begin(); it != eventList->end(); ++it){
		unsigned int ct_coincidences = 1; 	// define counter for the number of channels that show coincident signals
//		cout << "Fr " << it->frame_number << "\t" << "Channel " << it->channel << endl;
		list<stic3_data_t>::iterator it2 = it; 
		++it2;
		unsigned int ch_A;
		unsigned int ch_B;
		unsigned int ch_C;
		while(it2 != eventList->end()){			// go through other entries with second iterator
			if(it->channel != it2->channel){	// if channel is different to first channel
				int time_diff = (it->time)-(it2->time);        
				if( (time_diff < 32 ) && (time_diff > -32) ){
//				if (it->T_CC == it2->T_CC){	// check coincidence
//					cout << "First two channels: "<< it->time << "\t" << it2->time << endl;
					ct_coincidences++;	// increase counter for this loop
//					cout << "Found " << ct_coincidences << " coincidences." << endl;
					ch_A = it->channel;
					ch_B = it2->channel; 
					list<stic3_data_t>::iterator it3 = it2;	// define third iterator. 
					++it3;
					while(num_ch>2 && it3 != eventList->end()){
						if (it3->channel != ch_A && it3->channel != ch_B){
//							if(it3->T_CC == it->T_CC){
							int time_diff_2 = (it->time)-(it3->time);        
							if( (time_diff_2 < 32 ) && (time_diff_2 > -32) ){
								ct_coincidences++;							// increase counter for this loop
//								cout << "Found " << ct_coincidences << " coincidences." << endl;
								ch_C = it3->channel;
								list<stic3_data_t>::iterator it4 = it3;	// define fourth iterator. 
								++it4;
								while(num_ch>3 && it4 != eventList->end()){
									if (it4->channel != ch_A && it4->channel != ch_B && it4->channel != ch_C){
//										if(it4->T_CC == it->T_CC){
										int time_diff_3 = (it->time)-(it4->time);        
										if( (time_diff_3 < 32 ) && (time_diff_3 > -32) ){
											ct_coincidences++;
											/*for (int n = 4; n>1; n--){
												for(int i = 0; i<n-1; i++){
													
												}
											}*/
											int cn[4] = {it->channel,it2->channel,it3->channel,it4->channel};
											int times[4] = {it->time, it2->time, it3->time, it4->time};	// Calculate time differences
											for(int n = 4; n>1; n--){					// bubble sort algorithm	
												for(int i = 0; i<n-1; i++){				// sort with increasing time
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

											t_diff_histo_1->Fill(times[1]-times[0]);
											t_diff_histo_2->Fill(times[2]-times[1]);
											t_diff_histo_3->Fill(times[3]-times[2]);
											t_diff_histo_4->Fill(times[3]-times[0]);
											first_coincident_channel->Fill(cn[0]);
//					cout << "Channel ordering: " << cn[0] << "\t" << cn[1] << "\t" << cn[2] << "\t" << cn[3] << endl;
//					cout << "Times: " << times[0] << "\t" << times[1] << "\t" << times[2] <<  "\t" << times[3] << endl;
//					cout << "Coincidences: " <<times[1]-times[0]<< "\t" << times[2]-times[1]<< "\t" <<times[3]-times[2] <<  "\t" <<times[0]-times[3] << endl;
										}
									}
									if (ct_coincidences > 3){
										it4 = eventList->erase(it4);
										break;
									}else{
										++it4;
									}
								}		
							}
						}
						if(ct_coincidences > 2){
							it3 = eventList->erase(it3);
							break;
						}else{
							++it3;
						}
					}
				}
			}
			if (ct_coincidences > 1) {
				histo_coin->Fill(ct_coincidences);		
								// at this point all coincidences for this event have been found (if present)
				it2 = eventList->erase(it2);	// delete this entry to prevent double counting
				break;				// break out of while loop
			}else{
				++it2;
			}
		}
		if ( ct_coincidences == 1 ) {
			histo_coin->Fill(1);
//			cout << "Found " << ct_coincidences << " coincidence." << endl;
		}
	}

}



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
