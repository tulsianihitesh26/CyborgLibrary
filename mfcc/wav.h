/*
 * wav.h
 *
 *  Created on: Mar 7, 2014
 *      Author: sharath
 */

#ifndef WAV_H_
#define WAV_H_

#ifdef WAVGLOBAL
#define WAVGLOBAL
#else
#define WAVGLOBAL extern
#endif

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include "../cyborg/commonUtils.h"
using namespace std;

/*! Header format of WAV file
 header size-44 bytes: char-1 byte, short-2 bytes, long-4 bytes */
typedef struct {
	char ChunkID[4];
	int ChunkSize;
	char Format[4];
	char Subchunk1ID[4];
	int Subchunk1Size;
	short AudioFormat;
	short NumChannels;
	int SampleRate;		//Sampling Rate
	int ByteRate;
	short BlockAlign;
	short BitsPerSample;
	char Subchunk2ID[4];
	int Subchunk2Size;		//Number of bytes in data
} WAV_HEADER;

/*! Structure for data samples read from wav file */
typedef struct {
	int DataSize; /*!< number of data samples in audio file */
	int SampleRate; /*!< sampling rate of audio samples */
	NumberArray Data; /*!< */
} WAV_SIGNAL;

/* Structure WAV: contains
 1) the file pointer to the wav file.
 2) The header data structure
 */
typedef struct {
	FILE *fp;
	WAV_HEADER hdr;
} WAV;
/*one shot processing*/WAVGLOBAL WAV_SIGNAL *WavLoadFile(
		const char * InputFileName,
		WAV_SIGNAL *WDat);
WAVGLOBAL WAV_SIGNAL *FreeWAV_SIGNAL(
		WAV_SIGNAL *foo);
void getWavChunk(
		char *InputFileName,
		WAV_SIGNAL *WDat,
		Number startTime_s,
		Number endTime_s);
/*frame by frame processing*/WAVGLOBAL WAV WavFileOpen(
		char *fname,
		int mode);
WAVGLOBAL void WavFileClose(
		WAV w);
WAVGLOBAL void WavFrameRead(
		WAV w,
		Number *arr,
		int fno,
		int fsize,
		int osize);
WAVGLOBAL void WavFrameReadnew(
		WAV w,
		Number *arr,
		int fsize,
		int osize);
WAVGLOBAL int WavGetNumberOfFrames(
		WAV w,
		int fsize,
		int osize);
WAVGLOBAL int WavToSamples(
		WAV w,
		Number time);
//WAVGLOBAL Number WavToTime (WAV w, int samples);
WAVGLOBAL Number FrameNumToTime(
		int frame,
		Number timestep,
		Number framesize);
WAVGLOBAL void WavPrintHeader(
		WAV w);
WAVGLOBAL void WavHeaderWrite(
		WAV w);
WAVGLOBAL void WavHeaderInit(
		WAV *w,
		int n,
		int ch,
		int sr,
		int bps);
WAVGLOBAL void WavFrameWrite(
		WAV w,
		Number *arr,
		int fsize);
//WAVGLOBAL Number WavGetMaxima (WAV w);

#endif /* WAV_H_ */
