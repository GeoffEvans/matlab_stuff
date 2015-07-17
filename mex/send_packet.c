/*=================================================================
 *      Sends ethernet packets for the synthesiser FPGA
 *
 *      The calling syntax is:
 *
 *        error_code = send_packet(destination_address, sender_address, packet)
 *
 *      Error codes are
 *          0 - success
 *          1 - cannot find devices
 *          2 - cannot open adaptor
 *          3 - unable to send packet
 *
 *      This is a MEX-file for MATLAB.
 *=================================================================*/

#include <math.h>
#include "pcap.h"
#include "mex.h"

unsigned int sendPacket(unsigned char *dest_addr, unsigned char *send_addr, unsigned char *data, int packet_length){
	unsigned char *packet = (unsigned char*)malloc(packet_length);

	pcap_if_t *alldevs;
	pcap_t *fp;
	char errbuf[PCAP_ERRBUF_SIZE];
	int i;

	if (pcap_findalldevs(&alldevs, errbuf) == -1) {
		return 1; // cannot find devices
	}

	if ((fp = pcap_open_live(alldevs[0].name,
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
	for (i = 6; i < 12; i++) {
		packet[i] = send_addr[i - 6];
	}     
	for (i = 12; i < packet_length; i++) {
		packet[i] = data[i - 12];
	}
	if (pcap_sendpacket(fp,	// Adapter
		packet,				// buffer with the packet
		packet_length		// size
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
	double error_code;
	unsigned char *data, *dest_addr, *send_addr;
	int packet_length;
	
	errorCheck(nlhs, nrhs);
    
    dest_addr = mxGetPr(prhs[0]);
    send_addr = mxGetPr(prhs[1]);
    data = mxGetPr(prhs[2]);
    packet_length = mxGetN(prhs[2]) + 12;
    
    error_code = sendPacket(dest_addr, send_addr, data, packet_length);   
    plhs[0] = mxCreateDoubleScalar(error_code);
    return;
}