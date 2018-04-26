#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdint.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include <cstdio>
#include <signal.h>
#include <dirent.h>
#include <sys/stat.h>
#include <stdio.h>

using namespace std;

int limit = 1500; 	// sets a lower limit for the number of coincidences a parameter configuration must produce in order to be taken into account. 
			// Suggestion: Measurement_Time * Frequency * 0.01 -> at least 1 percent of events need to be detected by both channels

int main(int argc, char *argv[]){

	if(argc<4){cout << "wrong use of parameters" << endl; return -1;}

	char* fileLocation = argv[1];
	int channel_1 = atoi(argv[2]);
	int channel_2 = atoi(argv[3]);

	char input[128];
	sprintf(input, "%s%d_%d_Calib_Results.txt", fileLocation, channel_1, channel_2);

	ifstream in(input);
	if(!in){cout << input << ": file does not exist." << endl; return -1;}

	char line[1000];
	in.getline(line, 1000);

	double best_FWHM=1000000.0;		// just pre-initialize
	double best_FWHM_error=0;
	unsigned int highest_num_coins = 0;	// just pre-initialize
	int bestParameters[2];			// first entry: Best TThresh, second entry: Best EThresh

	int TThresh;
	int EThresh;
	double mean;
	double mean_error;
	double FWHM;
	double FWHM_error;
	unsigned int num_coins;

	while(!in.eof()){
		in >> TThresh >> EThresh >> mean >> mean_error >> FWHM >> FWHM_error >> num_coins;
		if (FWHM < best_FWHM && num_coins > limit) {
			best_FWHM = FWHM;
			best_FWHM_error = FWHM_error;
			bestParameters[0]=TThresh;
			bestParameters[1]=EThresh;
		}
		if (num_coins > highest_num_coins) highest_num_coins = num_coins;
		}
	cout << "Best parameter configuration:\t" << bestParameters[0] <<"\t" << bestParameters[1] << endl;
	cout << "CTR:\t" << best_FWHM << "\t+/-" << best_FWHM_error << endl;
	in.close();
	return 0;
}
