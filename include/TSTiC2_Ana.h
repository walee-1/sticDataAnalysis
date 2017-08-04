/*
 * Data analysis class for STiC2 readout data
 *
 */

#ifndef __TSTIC2_ANA_H_
#define __TSTIC2_ANA_H_


#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "stdio.h"
#include "TSystem.h"

#ifndef __CINT__
#include "EventDict.h"
#endif

#include "EventType.h"
#include "TFile.h"
#include "TKDE.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TPolyMarker.h"
#include "Math/RootFinder.h"
#include "Math/WrappedTF1.h"
#include <list>
#include <map>
#include <sstream>
#include <iostream>
#include <TF2.h>
//progress bar for the impatient
#include "progressbar.h"
#include "tdc_ch_values.h"

#include "TStyle.h"
#include "TROOT.h"



//Compare a stic3_data_t to a channel
inline bool operator==(const stic3_data_t& lhs, const short& rhs) {return lhs.channel==rhs;}
inline bool operator==(const short& lhs, const stic3_data_t& rhs) {return lhs==rhs.channel;}


#ifdef __CINT__
	typedef struct energy_info
#else
//STRUCT CONTAINING MEAN AND SIGMA OF THE COINCIDENCE ENERGY PEAK
struct energy_info{
	unsigned int mean;
	float sigma;
};
#endif

class TSTiC2_Ana{

	public:
		//Functions

		TSTiC2_Ana(TTree *tree);
		~TSTiC2_Ana(); //!< Destructor, removing the allocated root histograms

		//! Generate the energy histograms
		void gen_hist_energy();

		//! Generate the energy histograms
		void gen_dnl();
		
		//! Returns the energy histogram for a channel
		/*!
		 * \param channel The STiC channel number
		 * \return The TH1F energy histogram
		 */
		TH1F* get_hist_energy(unsigned int channel);//!< Return the energy histogram for the channel

		//! Returns the coincidence histogram for a channel combination
		/*!
		 * \param ch1 The first STiC channel number
		 * \param ch2 The second STiC channel number
		 * \return The TH1F energy histogram
		 */
		TH1F* get_hist_cspect(unsigned short ch1, unsigned short ch2);

		
		/*! Return the realtime bin distribution of specific channel
		 * \param channel Channel for which to return the distribution
		 */
		TH1F* get_realtime_distribution(unsigned int channel){return map_dnl_correct[channel]->get_realtimedist();};

		/*! Return the DNL bin distribution of specific channel
		 * \param channel Channel for which to return the distribution
		 */
		TH1F* get_DNL_distribution(unsigned int channel){return map_dnl_correct[channel]->get_DNL();};

		/*! Return the DNL bin distribution of specific channel
		 * \param channel Channel for which to return the distribution
		 */
		TH1F* get_INL_distribution(unsigned int channel) {return map_dnl_correct[channel]->get_INL();};

		/*! Find the 511 keV Peaks in the spectrum of the channel
		 * \return 0 if successful, -1 on error
		 */
		int find_peak(unsigned int channel);


		/*! Build the coincidences between the two channels
		 * \param ch1 First channel
		 * \param ch2 Second channel
		 * \return TH1F* histogram if successful, NULL on error
		 */
		void search_coincidence();

		//! Print energy information table
		void print_einfo();

		//! Fit the spectra an print the found resolutions
		void fit_cspects();

		//!Use the supplied user function to evaluate the events
		//void* eval_user_func( void* (*user_func)(list<stic3_data_t>*ev_list) );


		//!Use the supplied user function to evaluate the events
		/*!
		 *\param channel Channel number
		 *\param side which side the channel is on (0: side A, 1: opposite of side A)
		 */
		int add_coincidence_channel(unsigned short channel, unsigned short side );

		//! Function to measure the event period of a channel
		/*!
		 * This function calculates the period for all channels defined on side A
		 */
	private:
		struct energy_info{
			unsigned int mean;
			float sigma;
		};


		char outfile[255];
		unsigned int time0[64];
		unsigned int frame0[64];

		bool DNL_correction_present;
		
		ProgressBar *pbar;
		//Functions

		//!Transfering a tree into the event_seq map
		//!Processing the events in the tree, grouping by frames and calling a function
		//!for evaluation of the events in the current frame
		/*!
		 * \param tree TTree with the dump data from STiC2
		 * \param eval_function Function to be called for the list of events
		 */
		void tree_walk(void* (TSTiC2_Ana::*eval_function)(list<stic3_data_t> *ev_list));

		//! Allocate a new energy histogram for the channel
		TH1F* new_hist_energy(unsigned int channel);

		//! Allocate a new timing spectrum histogram for the channels
		TH1F* new_hist_cspect(unsigned short ch1, unsigned short ch2);


		//! Function to compare the channel values of two events, used to sort lists of events
		/*!
		 *\return True if channel of first is smaller than channel of second
		 */
		bool compare_channel(stic3_data_t& first, stic3_data_t& second);

		//! Function to process the energy information in the tree and generate the energy spectra
		void* eval_energy( list<stic3_data_t> *ev_list);

		//! Function to search for coincidences between all opposite placed channels
		void* eval_coincidence( list< stic3_data_t > *ev_list);

		//! Function to process the energy information in the tree and generate the energy spectra
		void* eval_dnl( list<stic3_data_t> *ev_list);

		//Variables

		std::map<unsigned short, TH1F*> map_energy_hist; //!< Container mapping a channel to its energy histogram
		std::map<unsigned short, TSpectrum*> map_peakfinder; //!<Container mapping a channel to the corresponding TSpectrum
		std::map<unsigned short, struct energy_info> map_einfo; //!< Container for information about the 511keV peaks

		std::map<unsigned short, std::map<unsigned short, TH1F*> > map_cspect;	//!< Container with the generated coincidence spectra

		std::list<unsigned short> channels_A;
		std::list<unsigned short> channels_B;

		std::map<unsigned short, tdc_ch_values*> map_dnl_correct; //!< Container for TDC DNL values and correction


		TTree* event_tree; //!< Original tree containing the recorded STiC2 event data


};

typedef void* (TSTiC2_Ana::* func_ptr)( list< stic3_data_t > *ev_list);

#endif
