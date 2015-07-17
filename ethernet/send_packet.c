/*=================================================================
 *      Sends ethernet packets for the synthesiser FPGA
 *
 *      The calling syntax is:
 *
 *        error_code = send_packet(destination_address, sender_address, packet)
 *
 *      Error codes are
 *          0 - success
 *          1 - dest_addr or send_addr wrong length
 *          2 - cannot find devices
 *          3 - unable to send packet
 *
 *      This is a MEX-file for MATLAB.
 *=================================================================*/

#include <math.h>
#include "pcap.h"
#include "mex.h"

unsigned int sendPacket(unsigned char *dest_addr, unsigned char *send_addr, unsigned char *data){
	int length = sizeof(data);
	unsigned char *packet = (unsigned char*) malloc(length + 14);
    
    
	pcap_if_t *alldevs;
	pcap_if_t d;
	pcap_t *fp;
	char errbuf[PCAP_ERRBUF_SIZE];
	int i;

	if (pcap_findalldevs(&alldevs, errbuf) == -1) {
		return 2; // cannot find devices
	}
	d = alldevs[0];

	/* Open the adapter */
	if ((fp = pcap_open_live(d.name,		// name of the device
		65536,			// portion of the packet to capture. It doesn't matter in this case 
		1,				// promiscuous mode (nonzero means promiscuous)
		1000,			// read timeout
		errbuf			// error buffer
		)) == NULL)
	{
		return 2; // unable to open the adapter
	}

	for (i = 0; i < 6; i++) {
		packet[i] = dest_addr[i];
	}

	for (i = 0; i < 6; i++) {
		packet[i + 6] = send_addr[i];
	}
    char lengthHex[3];
	if (length <= 0xFF)
	{
		sprintf(&lengthHex[0], "%02x", length);
	}
	packet[12] = lengthHex[0];
	packet[13] = lengthHex[1];

	for (i = 0; i < length; i++) {
		packet[i + 14] = data[i];
	}

	if (pcap_sendpacket(fp,	// Adapter
		packet,				// buffer with the packet
		100					// size
		) != 0)
	{
		return 3; // unable to send packet
	}
	pcap_close(fp);
	return 0;
}

void errorCheck(int nlhs, int nrhs)
{
    if (nrhs != 3) {
        mexErrMsgTxt("I receive 3 args.");
    } else if (nlhs != 1) {
        mexErrMsgTxt("I return an error code, which you need to assign to a variable.");
    }
    return;
}

//main entry point
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    unsigned char *dest_addr, *send_addr, *data;
    double *error_code;
    //errorCheck(nlhs, nrhs);

    mxGetString(prhs[0], dest_addr, 6);
    mxGetString(prhs[1], send_addr, 6);
    mxGetString(prhs[2], data, 1486);
    
    error_code = mxGetPr(plhs[0]);
    
    error_code[0] = sendPacket(dest_addr, send_addr, data);   
    return;
}


