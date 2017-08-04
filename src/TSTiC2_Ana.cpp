#include "TSTiC2_Ana.h"
#include <list>
#include <map>
#include "TSpectrum.h"
#include <fstream>

//#define DEBUG
//#define DEBUG_COIN 
//#define DNL_DEBUG
#define ECUT

extern unsigned const int MAPPING_N_BINS;
extern const double PLL_REFCLK_FREQ;

using namespace std;

const float Ecut = 1.5;
const float binSize = 1.0e12/PLL_REFCLK_FREQ/32.0; // in ps


bool flagFirst[64]={true};

//Constructor
TSTiC2_Ana::TSTiC2_Ana(TTree* tree)
{
	event_tree=tree;
	DNL_correction_present = false;
}

TSTiC2_Ana::~TSTiC2_Ana()
{

	map<unsigned short, TH1F*>::iterator it;
	for(it=map_energy_hist.begin();it!=map_energy_hist.end();++it)
	{
		delete it->second; //Delete the allocated histograms
	}

	map<unsigned short, TSpectrum*>::iterator it_p;
	for(it_p=map_peakfinder.begin();it_p!=map_peakfinder.end();++it_p)
	{
		delete it_p->second; //Delete the allocated TSpectrum
	}
}



int TSTiC2_Ana::add_coincidence_channel(unsigned short channel, unsigned short side )
{
  printf("in add coincidence_channel...\n");
	list<unsigned short>::iterator find_chan;

	if(side != 0 && side != 1){
		std::cout << "add_coincidence_channel(): Invalid side specified\n";
		return -1;
	}

	if(side == 0){
		if( find(channels_A.begin(),channels_A.end(),channel) != channels_A.end()){
			cout << "add_coincidence_channel(): Channel already added to side A\n";
			return -1;
		}
		if( find(channels_B.begin(),channels_B.end(),channel) != channels_B.end()){
			cout << "add_coincidence_channel(): Channel already added to opposite side B\n";
			return -1;
		}
		printf("Add %d to side %d\n",channel,side);
		channels_A.push_back(channel);
	}

	if(side == 1){
		if( find(channels_B.begin(),channels_B.end(),channel) != channels_B.end()){
			cout << "add_coincidence_channel(): Channel already added to side B\n";
			return -1;
		}
		if( find(channels_A.begin(),channels_A.end(),channel) != channels_A.end()){
			cout << "add_coincidence_channel(): Channel already added to opposite side A\n";
			return -1;
		}
		printf("Add %d to side %d\n",channel,side);
		channels_B.push_back(channel);
	}
	return 0;
}

void* TSTiC2_Ana::eval_coincidence( list< stic3_data_t > *ev_list)
{

	if (ev_list->size() < 2) return NULL; //ONLY SINGLE EVENTS IN THE FRAME
	//values of the current energy cut
	unsigned int maxE_ch1, maxE_ch2;
	unsigned int minE_ch1, minE_ch2;
	unsigned short ch1,ch2;

	unsigned int time0,energy0,point;

	list<unsigned short>::iterator it_ch1;
	list<unsigned short>::iterator it_ch2;

	list<stic3_data_t>::iterator findit;
	double ctr;
	//DOUBLE LOOP, SHORTENED TO HAVE LESS INDENTATION
	for(it_ch1=channels_A.begin();it_ch1!=channels_A.end();++it_ch1){
		ch1=*it_ch1;
		if (find(ev_list->begin(),ev_list->end(),ch1) == ev_list->end()) continue; //Next if ch1 is not present

		for(it_ch2=channels_B.begin();it_ch2 != channels_B.end();++it_ch2){

			ch2=*it_ch2;

			if( map_energy_hist.count(ch1) == 0 || map_energy_hist.count(ch2)==0 )
			{
#ifdef DEBUG
				std::cout << "Channel " << ch1 << " or " \
					<< ch2 << " is missing, aborting coincidence search\n";
#endif
				continue;

			}
			else if( map_einfo.count(ch1) == 0 || map_einfo.count(ch2)==0 )
			{
#ifdef DEBUG
				std::cout << "No energy information present!";
				std::cout <<" Did you run the find_peak function?" \
					<< ch1 << " " << ch2 << std::endl;
#endif
				continue;
			}


			//GET THE ENERGY CUTS FOR THE CHANNELS
			maxE_ch1 = map_einfo[ch1].mean+Ecut*map_einfo[ch1].sigma;
			minE_ch1 = map_einfo[ch1].mean-Ecut*map_einfo[ch1].sigma;
			
			maxE_ch2 = map_einfo[ch2].mean+Ecut*map_einfo[ch2].sigma;
			minE_ch2 = map_einfo[ch2].mean-Ecut*map_einfo[ch2].sigma;

			if (find(ev_list->begin(),ev_list->end(),ch2) != ev_list->end())
			{

				findit=find(ev_list->begin(),ev_list->end(),ch1);
#ifdef ECUT
				//APPLY E-CUT TO CHANNEL 1
				if(findit->energy > maxE_ch1 || findit->energy< minE_ch1) continue;
#endif

 				time0=findit->time;

				energy0=findit->energy;
				findit=find(ev_list->begin(),ev_list->end(),ch2);
#ifdef ECUT
				//APPLY E-CUT TO CHANNEL 2
				if(findit->energy > maxE_ch2 || findit->energy < minE_ch2) continue;
#endif
				
				//Allocate a new histogram if required
				if(map_cspect.count(ch1)==0) new_hist_cspect(ch1,ch2);
				if(map_cspect[ch1].count(ch2)==0) new_hist_cspect(ch1,ch2);
				
				if(DNL_correction_present==true){
				  map_cspect[ch1][ch2]->Fill((int)map_dnl_correct[ch1]->get_dith_value(time0) - (int)map_dnl_correct[ch2]->get_dith_value(findit->time) );
				}else{
				  map_cspect[ch1][ch2]->Fill( (int)time0 - (int)(findit->time) );
				}
			}
			
		}
	}
}

void TSTiC2_Ana::search_coincidence()
{
	TSTiC2_Ana* cp = this;
	void* (TSTiC2_Ana::*eval)(list< stic3_data_t >*) = &TSTiC2_Ana::eval_coincidence;
	tree_walk(eval);
	
}

void TSTiC2_Ana::fit_cspects()
{
  //std::cout << "Not yet implemented\n";

  std::map< unsigned short, map < unsigned short, TH1F* > >::iterator it_ch1;
  std::map < unsigned short, TH1F* >::iterator it_ch2;

  TH1F* current_hist;
  //TF1 *f;
  unsigned int max_bin;
  int low_edge;
  double para[4];
  float NsigmaT=3.0;
  
  for(it_ch1 = map_cspect.begin();it_ch1 != map_cspect.end(); ++it_ch1)
    {
      for( it_ch2 = it_ch1->second.begin();it_ch2!=it_ch1->second.end(); ++it_ch2)
	{
	  current_hist=it_ch2->second;
	  
	  if ( current_hist->GetEntries() > 30 ) // if enough entries 
	    {
	     
	      TF1 *f = new TF1("f","gaus");
	      //GET THE PEAK POSITION
	      max_bin=current_hist->GetMaximumBin();
	      low_edge=current_hist->GetBinLowEdge(max_bin);
	      //if(current_hist->GetBinContent(max_bin)<30) continue;
	      
	      f->SetRange(low_edge-20,low_edge+20);
	      
	      current_hist->Fit("f","QMR");
	      f->GetParameters(&para[0]);
	      f->SetRange(para[1]-5*para[2], para[1]+5*para[2]);
	      current_hist->GetXaxis()->SetRangeUser(para[1]-200,para[1]+200);
	      
	      //current_hist->Draw();
	      current_hist->Fit("f","QMR");
	      
	      if(f!=NULL )
		{
		  //current_hist->SetAxisRange(f->GetParameter(1)-400, f->GetParameter(1)+400);
		  
		  if(DNL_correction_present == true) //print out right information for dnl correction case
		    {
		      std::cout << "mapping n bin:"<< MAPPING_N_BINS <<endl;
		      std::cout << "Attention: New bin size: "<< binSize/(MAPPING_N_BINS/32) << " ps!!!!!"<<endl;
		      std::cout << "Coincidence between Ch" << it_ch1->first << " and Ch" << it_ch2->first;
		      std::cout << " with event number: " << current_hist->GetEntries();
		      std::cout << std::endl;
		      
		      std::cout << "Mean difference: " << f->GetParameter(1)*binSize/(MAPPING_N_BINS/32) << " ps"<< std::endl;
		      std::cout << "FWHM: " <<  f->GetParameter(2)*2.35*binSize/(MAPPING_N_BINS/32) << " ps"<< std::endl;
		      std::cout << "Sigma: " <<  f->GetParameter(2)*binSize/(MAPPING_N_BINS/32) << " ps"<< std::endl;
		      std::cout << std::endl;
		      
		      std::cout << "FWHM_error: " <<  f->GetParError(2)*2.35*binSize/(MAPPING_N_BINS/32)<< " ps"<< std::endl;
		      std::cout << "Chi square reduced: " <<  f->GetChisquare()/f->GetNDF()<< " prob("<<f->GetProb()<<")"<<std::endl;
		      std::cout << std::endl;
		      
		      std::ofstream logfile;
		      logfile.open("cspect_fit.log", std::ofstream::out | std::ofstream::app );
		      logfile <<  it_ch1->first << " " << it_ch2->first << " " << f->GetParameter(1)*binSize/(MAPPING_N_BINS/32) << " " << f->GetParameter(2)*binSize/(MAPPING_N_BINS/32)*2.35 ;
		      logfile << " " << f->GetParError(2)*binSize/(MAPPING_N_BINS/32)*2.35;
		      logfile << " " << f->GetChisquare()/f->GetNDF();
		      logfile << " " << current_hist->GetEntries();
		      		      
		      logfile << std::endl;
		      logfile.close();
		    }
		  else
		    {
		      std::cout << "Coincidence between Ch" << it_ch1->first << " and Ch" << it_ch2->first;
		      std::cout << " with event number: " << current_hist->GetEntries();
		      std::cout << std::endl;
		      
		      std::cout << "Mean difference: " << f->GetParameter(1)*binSize << " ps"<< std::endl;
		      std::cout << "FWHM: " <<  f->GetParameter(2)*2.35*binSize << " ps"<< std::endl;
		      std::cout << "Sigma: " <<  f->GetParameter(2)*binSize << " ps"<< std::endl;
		      std::cout << std::endl;
		      
		      std::cout << "FWHM_error: " <<  f->GetParError(2)*2.35*binSize<< " ps"<< std::endl;
		      std::cout << "Chi square reduced: " <<  f->GetChisquare()/f->GetNDF()<< " prob("<<f->GetProb()<<")"<<std::endl;
		      std::cout << std::endl;
		      
		      std::ofstream logfile;
		      logfile.open("cspect_fit.log", std::ofstream::out | std::ofstream::app );
		      logfile <<  it_ch1->first << " " << it_ch2->first << " " << f->GetParameter(1)*binSize << " " << f->GetParameter(2)*binSize*2.35 ;
		      logfile << " " << f->GetParError(2)*binSize*2.35;
		      logfile << " " << f->GetChisquare()/f->GetNDF();
		      logfile << " " << current_hist->GetEntries();
		      		      
		      logfile << std::endl;
		      logfile.close();

		    }

		} // end of if(f!=NULL)

	    } // end of  if ( current_hist->GetEntries() > 30 )
	}
    }
}

int TSTiC2_Ana::find_peak(unsigned int channel)
{
	int n_peaks,i;
	float_t *x_positions;
	float_t fXpos[20],fYpos[20];
	TH1F* current_hist;

	if ( map_energy_hist.count(channel) == 0)
	{
		std::cout << "find_peak(" << channel << "): Channel not found\n";
		return -1;
	}
	else current_hist = map_energy_hist[channel];

	std::cout << "Searching for photonpeak in energy spectrum of channel " << channel << std::endl;
	int bins=current_hist->GetNbinsX()-1;
	float source[bins];
	float dest[bins];
	
	for(i=0;i<bins;i++) source[i]=current_hist->GetBinContent(i+1);
	
	//Search with high resolution, use Markov average window to adjust number of found peaks in the algorithm
	map_peakfinder[channel]->SearchHighRes(source,dest,bins,10,5,kTRUE,3,kTRUE,3);

	n_peaks=map_peakfinder[channel]->GetNPeaks();
//	std::cout << "Found " << n_peaks << " peaks in the energy spectrum of channel " << channel << std::endl;

	x_positions=map_peakfinder[channel]->GetPositionX();

	//Remove old polymarkers when calling the function multiple times
	TPolyMarker* pm=(TPolyMarker*)(current_hist->GetListOfFunctions()->FindObject("TPolyMarker"));
	if(pm){
		current_hist->GetListOfFunctions()->Remove(pm);
		delete pm;
	}

	for(i=0;i<n_peaks;i++)
	{
	//	std::cout << "Peak " << ": " << x_positions[i] << std::endl << std::endl;
		fXpos[i]=current_hist->GetBinCenter(Int_t(x_positions[i]+0.5)+1);
		fYpos[i]=current_hist->GetBinContent(Int_t(x_positions[i]+0.5)+1);
	}

	TF1 *f;
	
	for(i=0;i<n_peaks;i++)
	  //for(i=0;i<1;i++)
	  {
	    printf("channel %d:par:%f %f\n",channel,fYpos[i],fXpos[i]);
	    int width=20;
	    //if(channel==53) width =50 ;
	    f= new TF1("gauss","gaus(0)",fXpos[i]-width,fXpos[i]+width);
	    f->SetParameter(0,fYpos[i]);
	    f->SetParameter(1,fXpos[i]);
	    f->SetParameter(2,0.2);
	    printf("npeaks(%d):par:%f %f\n",n_peaks,fYpos[i],fXpos[i]);
	    current_hist->Fit("gauss","M","",fXpos[i]-width,fXpos[i]+width);
	    current_hist->Fit("gauss","M","",fXpos[i]-Ecut*f->GetParameter(2),fXpos[i]+Ecut*f->GetParameter(2));
	    //int ret;
	    //cin >>ret;
	    if(f->GetParameter(2) < 100 && f->GetParameter(2) > 2) break;
	  }
	struct energy_info einfo;
	
	if( f!=NULL )
	{
	  einfo.mean=f->GetParameter(1);
	  einfo.sigma=f->GetParameter(2);
	  map_einfo[channel]=einfo;
	}
	else
	{
	  std::cout << "Warning: Channel " << channel << " has bad energy information\n";
	}

	//Add Polymarkers to the histogram
	pm = new TPolyMarker(n_peaks,fXpos,fYpos);
	current_hist->GetListOfFunctions()->Add(pm);
	pm->SetMarkerStyle(23);
	pm->SetMarkerColor(kRed);
	pm->SetMarkerSize(1.2);

	//f->Delete();
	//pm->Delete();
	
	return 0;
}


/*
 * Print the energy information table
 */
void TSTiC2_Ana::print_einfo()
{
	std::map<unsigned short,struct energy_info>::iterator it;
	for(it=map_einfo.begin();it!=map_einfo.end();++it)
	{
		std::cout << "Channel " << it->first;
		std::cout << " E-Mean: " << it->second.mean;
		std::cout << " E-Sigma: " << it->second.sigma;
		std::cout << endl;
	}
}


/*
 * Transfer the TTree into the event container
 */
void TSTiC2_Ana::tree_walk(void* (TSTiC2_Ana::*eval_function)(list<stic3_data_t> *ev_list))
{


	TSTiC2_Ana *cp = this;

	pbar = new ProgressBar();

	stic3_data_t* event=NULL;
	std::list<stic3_data_t> tmp_list; //TEMPORARY LIST FOR INSERTING IN THE MAP
	tmp_list.clear();

	TTree* tree=event_tree; //Pointer renaming...

	tree->SetBranchAddress("br",&event);

	unsigned int entries;
	unsigned int last_frame, frame_seq;
	frame_seq=0;

	tree->GetEntry(0);
	last_frame=event->frame_number; //Start with the first frame number


	entries=tree->GetEntries();

	pbar->start();
	for (int i=0; i<tree->GetEntries(); i++)
	{

		tree->GetEntry(i);
		if (last_frame != event->frame_number)
		{
			frame_seq++;


			(cp->*eval_function)(&tmp_list); //Call the eval function for the current list
			tmp_list.clear(); //Clear the map for the next frame

			last_frame = event->frame_number;
		}

		tmp_list.push_back(*event);	//COPY THE EVENT TO THE LIST
		if(i%2000==0) pbar->status(i,entries);
	}

	pbar->stop();
	delete pbar;

}


/*
 * Allocate a new cspect histogram if required
 */
TH1F* TSTiC2_Ana::new_hist_cspect(unsigned short ch1, unsigned short ch2)
{

#ifdef DEBUG
	std::cout << "Allocating new coincidence spectrum for channel " << ch1 << " to channel " << ch2 << std::endl;
#endif
	ostringstream hist_name;
	ostringstream hist_title;

	hist_name.clear();
	hist_name.str("");
	hist_name << "Coincidence" << ch1 << "_" << ch2 << "_";
	//hist_name << "Coincidence" << ch1 << "_" << ch2;
	hist_title.clear();
	hist_title.str("");
	hist_title << "T-Spectrum between " << ch1 << " and " << ch2;
	TH1F* cspect_hist;

	if (DNL_correction_present==false)
		cspect_hist=new TH1F(hist_name.str().c_str(),hist_title.str().c_str(),10000000,-5000000,5000000); //allocate the new histogram
	else
		cspect_hist=new TH1F(hist_name.str().c_str(),hist_title.str().c_str(),10000000,-5000000,5000000); //allocate the new histogram
#ifdef DEBUG
	std::cout<< "bin width of coincidence spectrum: " << cspect_hist->GetXaxis()->GetBinWidth(1) <<""<<std::endl;
#endif

	map_cspect[ch1][ch2]=cspect_hist;
	return cspect_hist;
}

/*
 * Allocate a new histogram if required
 */
TH1F* TSTiC2_Ana::new_hist_energy(unsigned int channel)
{

	ostringstream hist_name;
	ostringstream hist_title;

	TH1F* energy_spect;
	//Create histograms
	hist_name.clear();
	hist_name.str("");
	hist_name << "ToT" << channel;

	hist_title.clear();
	hist_title.str("");
	hist_title << "ToT Spectrum Ch" << channel;

	energy_spect = new TH1F(hist_name.str().c_str(),hist_title.str().c_str(),2000,0,2000);
	energy_spect->GetXaxis()->SetTitle("ToT [TDC Bins]");
	energy_spect->GetYaxis()->SetTitle("Entries");

	return energy_spect;
}

void* TSTiC2_Ana::eval_energy(list<stic3_data_t> *ev_list)
{

	list<stic3_data_t>::iterator it;
	for(it = ev_list->begin(); it!=ev_list->end(); ++it)
	{
		//Create a new histogram for the channel if required
		if (map_energy_hist.count(it->channel) == 0)
		{
			map_energy_hist[it->channel]=new_hist_energy(it->channel);
			map_peakfinder[it->channel]= new TSpectrum(5,1.0); //Allocate the corresponding peak finder
		}
		map_energy_hist[ it->channel ]->Fill( it->energy );
	}

}

/*
 * Calculate the DNL and INL plots, call tree_walk with eval_dnl
 */
void TSTiC2_Ana::gen_dnl()
{
  //#ifdef DEBUG
	std::cout << "Generating DNL correction tables" << std::endl;
	//#endif
	TSTiC2_Ana* cp = this;
	void* (TSTiC2_Ana::*eval)(list< stic3_data_t >*) = &TSTiC2_Ana::eval_dnl;
	tree_walk(eval);

	map<unsigned short, tdc_ch_values*>::iterator it_dnl;
	for(it_dnl = map_dnl_correct.begin(); it_dnl!=map_dnl_correct.end(); ++it_dnl)
	{
		it_dnl->second->update_histos();
	}

#ifdef DEBUG
	std::cout << "Done generating DNL correction tables" << std::endl;
#endif
	DNL_correction_present = true;
}


/*
 * Evaluation function for the TDC DNL and INL plots
 */
void* TSTiC2_Ana::eval_dnl(list<stic3_data_t> *ev_list)
{


	ostringstream hist_number;

	list<stic3_data_t>::iterator it;
	for(it = ev_list->begin(); it!=ev_list->end(); ++it)
	{
		//Create a new histogram for the channel if required
		if (map_dnl_correct.count(it->channel) == 0)
		{
			//Create histograms
			hist_number.clear();
			hist_number.str("");
			hist_number << it->channel;
			map_dnl_correct[it->channel]=new tdc_ch_values(0,31,hist_number.str().c_str());
		}
		map_dnl_correct[it->channel]->fill(it->T_fine);
		map_dnl_correct[it->channel]->fill(it->E_fine);
		
	}

}


/*
 * Generate the energy spectra from the map, call tree_walk with eval_energy
 */
void TSTiC2_Ana::gen_hist_energy()
{
  printf("generating hist_energy...\n");
	TSTiC2_Ana* cp = this;
	void* (TSTiC2_Ana::*eval)(list< stic3_data_t >*) = &TSTiC2_Ana::eval_energy;
	tree_walk(eval);
}


/*
 * Return the energy histogram if available else NULL
 */
TH1F* TSTiC2_Ana::get_hist_energy(unsigned int channel)
{
	if(map_energy_hist.count(channel) > 0) return map_energy_hist[channel];
	return NULL;
}

TH1F* TSTiC2_Ana::get_hist_cspect(unsigned short ch1, unsigned short ch2)
{
	if(map_cspect.count(ch1)==0) return NULL;
	if(map_cspect[ch1].count(ch2)==0) return NULL;
	return map_cspect[ch1][ch2];
}


/*
 * Comparison function between stic events which channel is smaller
 */
bool TSTiC2_Ana::compare_channel(stic3_data_t& first, stic3_data_t& second)
{
	return first.channel < second.channel;
}


//void* TSTiC2_Ana::eval_user_func(void* (*user_func)(list<stic3_data_t> *ev_list))
//{
//	//tree_walk(user_func);
//	return NULL;
//}
