/***************************************************
This is the matlab interface code to the grabber code.
It just wraps the grabber functions and does some error
conversion.

Written by Micah Richert.
07/14/2005
**************************************************/

#include "DDGrab.h"

// unfortunately to detect memory leaks, we have to put code in each .cpp file
#if defined(_DEBUG) && defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__,__LINE__)
#endif

#include "mex.h"

TCHAR str[200];

TCHAR* message(HRESULT hr)
{
	if (hr == S_OK)
	{
		return "";
	} else {
		if (AMGetErrorText(hr,str,200) != 0) return str;
		return "Unknown error";
	}
}

DDGrabber DDG;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || !mxIsChar(prhs[0])) mexErrMsgTxt("First parameter must be the command (a string)");

	char cmd[100];
	mxGetString(prhs[0],cmd,100);

	if (!strcmp("buildGraph",cmd))
	{
		if (nrhs < 2 || !mxIsChar(prhs[1])) mexErrMsgTxt("buildGraph: second parameter must be the filename (as a string)");
		if (nlhs > 0) mexErrMsgTxt("buildGraph: there are no outputs");
		int filenamelen = mxGetN(prhs[1])+1;
		char* filename = new char[filenamelen];
		if (!filename) mexErrMsgTxt("buildGraph: out of memory");
		mxGetString(prhs[1],filename,filenamelen);

		char* errmsg =  message(DDG.buildGraph(filename));
		delete[] filename;

		if (strcmp("",errmsg)) mexErrMsgTxt(errmsg);
	} else if (!strcmp("doCapture",cmd)) {
		if (nlhs > 0) mexErrMsgTxt("doCapture: there are no outputs");
		char* errmsg =  message(DDG.doCapture());
		if (strcmp("",errmsg)) mexErrMsgTxt(errmsg);
	} else if (!strcmp("getVideoInfo",cmd)) {
		if (nrhs < 2 || !mxIsNumeric(prhs[1])) mexErrMsgTxt("getVideoInfo: second parameter must be the video stream id (as a number)");
		if (nlhs > 6) mexErrMsgTxt("getVideoInfo: there are only 5 output values: width, height, rate, nrFramesCaptured, nrFramesTotal");

		unsigned int id = mxGetScalar(prhs[1]);
		int width,height,nrFramesCaptured,nrFramesTotal;
		double rate, totalDuration;
		char* errmsg =  message(DDG.getVideoInfo(id, &width, &height,&rate, &nrFramesCaptured, &nrFramesTotal, &totalDuration));

		if (strcmp("",errmsg)) mexErrMsgTxt(errmsg);

		if (nlhs >= 1) {plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[0])[0] = width; }
		if (nlhs >= 2) {plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[1])[0] = height; }
		if (nlhs >= 3) {plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[2])[0] = rate; }
		if (nlhs >= 4) {plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[3])[0] = nrFramesCaptured; }
		if (nlhs >= 5) {plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[4])[0] = nrFramesTotal; }
		if (nlhs >= 6) {plhs[5] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[5])[0] = totalDuration; }
	} else if (!strcmp("getAudioInfo",cmd)) {
		if (nrhs < 2 || !mxIsNumeric(prhs[1])) mexErrMsgTxt("getAudioInfo: second parameter must be the audio stream id (as a number)");
		if (nlhs > 7) mexErrMsgTxt("getAudioInfo: there are only 6 output values: nrChannels, rate, bits, nrFramesCaptured, nrFramesTotal, subtype");

		unsigned int id = mxGetScalar(prhs[1]);
		int nrChannels,bits,nrFramesCaptured,nrFramesTotal;
		GUID subtype;
		double rate, totalDuration;
		char* errmsg =  message(DDG.getAudioInfo(id, &nrChannels, &rate, &bits, &nrFramesCaptured, &nrFramesTotal, &subtype, &totalDuration));

		if (strcmp("",errmsg)) mexErrMsgTxt(errmsg);

		if (nlhs >= 1) {plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[0])[0] = nrChannels; }
		if (nlhs >= 2) {plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[1])[0] = rate; }
		if (nlhs >= 3) {plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[2])[0] = bits; }
		if (nlhs >= 4) {plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[3])[0] = nrFramesCaptured; }
		if (nlhs >= 5) {plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[4])[0] = nrFramesTotal; }
		if (nlhs >= 6) {plhs[5] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[5])[0] = subtype==MEDIASUBTYPE_PCM?0:subtype==MEDIASUBTYPE_IEEE_FLOAT?1:-1; }
		if (nlhs >= 7) {plhs[6] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[6])[0] = totalDuration; }
	} else if (!strcmp("getCaptureInfo",cmd)) {
		if (nlhs > 2) mexErrMsgTxt("getCaptureInfo: there are only 2 output values: nrVideo, nrAudio");

		int nrVideo, nrAudio;
		DDG.getCaptureInfo(&nrVideo, &nrAudio);

		if (nlhs >= 1) {plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[0])[0] = nrVideo; }
		if (nlhs >= 2) {plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[1])[0] = nrAudio; }
	} else if (!strcmp("getVideoFrame",cmd)) {
		if (nrhs < 3 || !mxIsNumeric(prhs[1]) || !mxIsNumeric(prhs[2])) mexErrMsgTxt("getVideoFrame: second parameter must be the audio stream id (as a number) and third parameter must be the frame number");
		if (nlhs > 2) mexErrMsgTxt("getVideoFrame: there are only 2 output value: data");

		unsigned int id = mxGetScalar(prhs[1]);
		int frameNr = mxGetScalar(prhs[2]);
		BYTE* data;
		int nrBytes;
		double time;
		int dims[] = {1,1};
		char* errmsg =  message(DDG.getVideoFrame(id, frameNr, &data, &nrBytes, &time));

		if (strcmp("",errmsg)) mexErrMsgTxt(errmsg);

		dims[0] = nrBytes;
		plhs[0] = mxCreateNumericArray(2, dims, mxUINT8_CLASS, mxREAL); // empty 2d matrix
		memcpy(mxGetPr(plhs[0]),data,nrBytes);
		free(data);
		if (nlhs >= 2) {plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[1])[0] = time; }
	} else if (!strcmp("getAudioFrame",cmd)) {
		if (nrhs < 3 || !mxIsNumeric(prhs[1]) || !mxIsNumeric(prhs[2])) mexErrMsgTxt("getAudioFrame: second parameter must be the audio stream id (as a number) and third parameter must be the frame number");
		if (nlhs > 2) mexErrMsgTxt("getAudioFrame: there are only 2 output value: data");

		unsigned int id = mxGetScalar(prhs[1]);
		int frameNr = mxGetScalar(prhs[2]);
		BYTE* data;
		int nrBytes;
		double time;
		int dims[] = {1,1};
		mxClassID mxClass;
		char* errmsg =  message(DDG.getAudioFrame(id, frameNr, &data, &nrBytes, &time));

		if (strcmp("",errmsg)) mexErrMsgTxt(errmsg);

		int nrChannels,bits,nrFramesCaptured,nrFramesTotal;
		GUID subtype;
		double rate, totalDuration;
		DDG.getAudioInfo(id, &nrChannels, &rate, &bits, &nrFramesCaptured, &nrFramesTotal, &subtype, &totalDuration);

		switch (bits)
		{
			case 4:
			{
				BYTE* tmpdata = (BYTE*)malloc(nrBytes*2);
				int i;

				for (i=0;i<nrBytes;i++)
				{
					tmpdata[i*2] = ((0x8&data[i])?-1:1)*(0x7&data[i]);
					tmpdata[i*2+1] = ((0x80&data[i])?-1:1)*((0x70&data[i])>>4);
				}

				free(data);
				data = tmpdata;
				mxClass = mxUINT8_CLASS;
				dims[0] = nrBytes*2;
                nrBytes = nrBytes*2;
				break;
			}
			case 8:
			{
		 		dims[0] = nrBytes;
				mxClass = mxUINT8_CLASS;
				break;
			}
			case 16:
			{
				mxClass = mxINT16_CLASS;
				dims[0] = nrBytes/2;
				break;
			}
			case 24:
			{
				int* tmpdata = (int*)malloc(nrBytes/3*4);
				int i;

				//I don't know how 24bit float data is organized...
				for (i=0;i<nrBytes/3;i++)
				{
					tmpdata[i] = (((0x80&data[i*3+2])?-1:0)&0xFF000000) | ((data[i*3+2]<<16)+(data[i*3+1]<<8)+data[i*3]);
				}

				free(data);
				data = (BYTE*)tmpdata;

				mxClass = mxINT32_CLASS;
				dims[0] = nrBytes/3;
                nrBytes = nrBytes/3*4;
				break;
			}
			case 32:
			{
				mxClass = subtype==MEDIASUBTYPE_PCM?mxINT32_CLASS:subtype==MEDIASUBTYPE_IEEE_FLOAT?mxSINGLE_CLASS:mxUINT32_CLASS;
				dims[0] = nrBytes/4;
				break;
			}
			case 64:
			{
				mxClass = mxDOUBLE_CLASS;
				dims[0] = nrBytes/8;
				break;
			}
			default:
			{
		 		dims[0] = nrBytes;
				mxClass = mxUINT8_CLASS;
				break;
			}
		}

		plhs[0] = mxCreateNumericArray(2, dims, mxClass, mxREAL); // empty 2d matrix
		memcpy(mxGetPr(plhs[0]),data,nrBytes);
		free(data);
		if (nlhs >= 2) {plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[1])[0] = time; }
	} else if (!strcmp("setFrames",cmd)) {
		if (nrhs < 2 || !mxIsDouble(prhs[1])) mexErrMsgTxt("setFrames: second parameter must be the frame numbers (as doubles)");
		if (nlhs > 0) mexErrMsgTxt("setFrames: has no outputs");
		int nrFrames = mxGetN(prhs[1]) * mxGetM(prhs[1]);
		int* frameNrs = new int[nrFrames];
		if (!frameNrs) mexErrMsgTxt("setFrames: out of memory");
		double* data = mxGetPr(prhs[1]);
		for (int i=0; i<nrFrames; i++) frameNrs[i] = data[i];

		DDG.setFrames(frameNrs, nrFrames);

		delete[] frameNrs;
	} else if (!strcmp("setTime",cmd)) {
		if (nrhs < 3 || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2])) mexErrMsgTxt("setTime: start and stop time are required (as doubles)");
		if (nlhs > 0) mexErrMsgTxt("setTime: has no outputs");

		DDG.setTime(mxGetScalar(prhs[1]), mxGetScalar(prhs[2]));
	} else if (!strcmp("setTrySeeking",cmd)) {
		if (nrhs < 2 || !mxIsDouble(prhs[1])) mexErrMsgTxt("setTrySeeking: must specify true/false or 1/0");
		if (nlhs > 0) mexErrMsgTxt("setTrySeeking: has no outputs");

		DDG.setTrySeeking(mxGetScalar(prhs[1]));
	} else if (!strcmp("setMatlabCommand",cmd)) {
		if (nrhs < 2 || !mxIsChar(prhs[1])) mexErrMsgTxt("setMatlabCommand: the command must be passed as a string");
		if (nlhs > 0) mexErrMsgTxt("setMatlabCommand: has no outputs");

		char * matlabCommand = (char*)calloc(100,1);
		mxGetString(prhs[1],matlabCommand,100);

		if (strlen(matlabCommand)==0) 
		{
			DDG.setMatlabCommand(NULL);
			free(matlabCommand);
		} else DDG.setMatlabCommand(matlabCommand);

	} else if (!strcmp("disableVideo",cmd)) {
		if (nlhs > 0) mexErrMsgTxt("disableVideo: there are no outputs");
		DDG.disableVideo();
	} else if (!strcmp("disableAudio",cmd)) {
		if (nlhs > 0) mexErrMsgTxt("disableAudio: there are no outputs");
		DDG.disableAudio();
	} else if (!strcmp("cleanUp",cmd)) {
		if (nlhs > 0) mexErrMsgTxt("cleanUp: there are no outputs");
		DDG.cleanUp();
	}
}