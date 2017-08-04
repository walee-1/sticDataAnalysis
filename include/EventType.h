#ifndef EVENT_TYPE_H__
#define EVENT_TYPE_H__

#include "TObject.h"

class stic3_data_t: public TObject
{
public:
	unsigned char	handleID;	//ChipID
	unsigned short	packetID;	//The Packet ID the event was transmitted in
	unsigned short	frame_number;	//The frame number of the event (for debug)
	unsigned short	channel;	//The channel number
	unsigned int	T_CC;		//The Time Coarse Counter value
	bool		T_badhit;	//The bad hit flag of the time measurement
	unsigned short	T_fine;		//The fine counter value of the time measurement
	unsigned int	E_CC;		//The coarse counter value of the energy measurement
	unsigned short	E_fine;		//The fine counter value of the energy measuremet
	bool		E_badhit;	//The bad hit flag of the energy measurement (not needed here)

	/* higher level entries */
	unsigned int	time;		//time including energy
	unsigned int	energy;		//energy (coarse counter only)			both select the correct CM/CS values
ClassDef(stic3_data_t,1);
};
//ClassImp(stic3_data_t);
#endif

