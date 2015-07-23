/*
 * wav.c
 *
 *  Created on: Mar 7, 2014
 *      Author: sharath
 */

#define WAVGLOBAL
#include "wav.h"

/* ############################################################# */
/* ------------------------------------------------------------- */
/*!LoadWavFile: reads data samples from wav file		 */
/* ------------------------------------------------------------- */
/* ############################################################# */

//read wavfile from given startSample and endSample
//assumes the wav file starts from sample 0
// code block added by sharath on 13/4/2013: adaptation of vishu's code block in WaveFunctions.cpp
void getWavChunk(
		char *InputFileName,
		WAV_SIGNAL *WDat,
		Number startTime_s,
		Number endTime_s) {

	FILE *fpInput;
	WAV_HEADER header;
	short *Signal;
	int i;

	/*!Read wav file header */
	fpInput = fopen(InputFileName, "rb");
	if (fpInput == NULL) {
		printf("%s File open failed.", InputFileName);
		exit(0);
	}
	fread(&header, sizeof(WAV_HEADER), 1, fpInput);
	WDat->DataSize = (int) header.Subchunk2Size / 2; /* for 16-bit (monophonic) data */
	WDat->SampleRate = (int) header.SampleRate; /* sampling rate (required for pitch calculation etc.) */

	if ((endTime_s + 1) < 0.05) // tweak to load the full song, for entire song endTime_s==-1
			{
		endTime_s = (double) WDat->DataSize / (double) WDat->SampleRate;
	}

	int dataSizeToRead = (endTime_s - startTime_s) * WDat->SampleRate;
	int startSample = (int) (startTime_s * WDat->SampleRate) * header.BlockAlign;

	if (!(startTime_s < endTime_s)) {
		printf("wrong indices to read start time= %g and end time=%g\n", startTime_s, endTime_s);
		printf("song length is %g second", (double) WDat->DataSize / WDat->SampleRate);
		exit(1);
	}

	WDat->DataSize = dataSizeToRead; //rewrite the Datasize, to the amount of data read

	WDat->Data.resize(dataSizeToRead, 0.0);
	Signal = (short *) calloc(dataSizeToRead, sizeof(short));

	fseek(fpInput, startSample, SEEK_CUR);

	/*!Read data samples from the wav file */
	fread(Signal, sizeof(short), dataSizeToRead, fpInput);
	for (i = 0; i < dataSizeToRead; i++) {
		WDat->Data[i] = (Number) Signal[i] / 32768;
	}

	fclose(fpInput);
	free(Signal);
	Signal = (short *) NULL;
}

WAV_SIGNAL *WavLoadFile(
		const char *InputFileName,
		WAV_SIGNAL *WDat) {
	FILE *fpInput;
	WAV_HEADER header;
	short *Signal;
	int i;

	/*!Read wav file header */
	fpInput = fopen(InputFileName, "rb");
	if (fpInput == NULL) {
		printf("%s File open failed.", InputFileName);
		exit(0);
	}
//	printf(" File open Successful.");

	fread(&header, sizeof(WAV_HEADER), 1, fpInput);
	WDat->DataSize = (int) header.Subchunk2Size / 2; /* for 16-bit (monophonic) data */
	WDat->SampleRate = (int) header.SampleRate; /* sampling rate (required for pitch calculation etc.) */

	WDat->Data.resize(WDat->DataSize, 0.0);
	Signal = (short *) malloc(sizeof(short) * WDat->DataSize);

	/*!Read data samples from the wav file */
	fread(Signal, sizeof(short), WDat->DataSize, fpInput);
	for (i = 0; i < WDat->DataSize; i++) {
		WDat->Data[i] = (Number) Signal[i];
	}

	fclose(fpInput);
	free(Signal);
	Signal = (short *) NULL;
	return WDat;
}

WAV_SIGNAL *FreeWAV_SIGNAL(
		WAV_SIGNAL *foo) {
	free(foo);
	return (WAV_SIGNAL *) NULL;
}
/********************************************************************/
/* WavFileOpen reads headers of the specified wav file
 mode: 1: read ---- 2: write
 */
/********************************************************************/
WAV WavFileOpen(
		char *fname,
		int mode) {
	WAV w;
	//FILE *fp;
	if (mode == 2) {
		w.fp = fopen(fname, "wb");
		if (w.fp == NULL) {
			printf("\nUnable to open wav file in write mode");
			exit(1);
		}
		fseek(w.fp, 0, 0);
		return w;
	}

	w.fp = fopen(fname, "rb");
	if (w.fp == NULL) {
		printf("\nError opening wav file.");
		exit(1);
	}
	fread(&(w.hdr), sizeof(WAV_HEADER), 1, w.fp);

	/*test*/
	//fp = fopen("wavinfo.txt", "wt");
	//fprintf(fp, "\nWavFileOpen Debug Dump:\n");
	//fprintf(fp, "\nChunkID       : %c%c%c%c",(w.hdr).ChunkID[0], (w.hdr).ChunkID[1], (w.hdr).ChunkID[2], (w.hdr).ChunkID[3]);
	//fprintf(fp, "\nFileSize(-8b) : %ld", w.hdr.ChunkSize);
	//fprintf(fp, "\nFormat        : %c%c%c%c", w.hdr.Format[0], w.hdr.Format[1], w.hdr.Format[2], w.hdr.Format[3]);
	//fprintf(fp, "\nFMT chunkID   : %c%c%c%c", w.hdr.Subchunk1ID[0], w.hdr.Subchunk1ID[1], w.hdr.Subchunk1ID[2], w.hdr.Subchunk1ID[3]);
	//fprintf(fp, "\nSubchunk1Size : %ld", w.hdr.Subchunk1Size);
	//fprintf(fp, "\nAudioFormat   : %d", w.hdr.AudioFormat);
	//fprintf(fp, "\nNumChannels   : %d", w.hdr.NumChannels);
	//fprintf(fp, "\nSampleRate    : %ld", w.hdr.SampleRate);
	//fprintf(fp, "\nByteRate      : %ld", w.hdr.ByteRate);
	//fprintf(fp, "\nBlockAlign    : %d", w.hdr.BlockAlign);
	//fprintf(fp, "\nBitsPerSample : %d", w.hdr.BitsPerSample);
	//fprintf(fp, "\nDATA chunkID  : %c%c%c%c", w.hdr.Subchunk2ID[0], w.hdr.Subchunk2ID[1], w.hdr.Subchunk2ID[2], w.hdr.Subchunk2ID[3]);
	//fprintf(fp, "\nData(in bytes): %ld", w.hdr.Subchunk2Size);
	//fclose(fp);
	return w;
}

/********************************************************************/
/* WavFileClose closes the wavfile
 */
/********************************************************************/
void WavFileClose(
		WAV w) {
	fclose(w.fp);
}
/********************************************************************/
/* WavFrameRead reads Framenumber fno into array arr.
 Frame size in samples: fsize
 overlap size in samples : osize
 */
/********************************************************************/
void WavFrameRead(
		WAV w,
		double *arr,
		int fno,
		int fsize,
		int osize) {
	int startpt, endpt, i, j, center;
	short int *signal;

	center = 0 + (fno) * (fsize - osize); /*the 0th frame is centered around 0 sec*/
	if (fsize % 2 != 0) /*fsize is odd*/
	{
		startpt = center - fsize / 2;
		endpt = center + fsize / 2;
	} else /*fsize is even*/
	{
		startpt = center - fsize / 2;
		endpt = center + fsize / 2 - 1;
	}
	//printf("\nfno = %d, startpt = %d, endpt = %d, tstep = %d", fno, startpt, endpt, fsize-osize);
	signal = (short int *) malloc(fsize * sizeof(short int));
	if (startpt < 0) {
		for (i = startpt, j = 0; i < 0; i++, j++)
			signal[j] = 0; /*prepad zeroes*/
		fseek(w.fp, (44 + (0 * (w.hdr.BitsPerSample / 8))), SEEK_SET);
		fread(&signal[j], sizeof(short), endpt, w.fp); /*read from file till the endpt*/
	} else {
		fseek(w.fp, (44 + (startpt * (w.hdr.BitsPerSample / 8))), SEEK_SET);
		fread(signal, sizeof(short), fsize, w.fp);
	}

	for (i = 0; i < fsize; i++)
		arr[i] = (double) ((Number) signal[i] / (Number) 32768);
	free(signal);

	return;
}
/********************************************************************/
/*A better implementation of wavframeread:
 Read in only the hopsize amount of data and append to the previously
 read data
 0-10*/
/********************************************************************/
void WavFrameReadnew(
		WAV w,
		double *arr,
		int fsize,
		int osize) {
	int i, j, endpt;
	static int fno = 0;
	short int *signal;

	if (fsize % 2 != 0) /*fsize is odd*/
		endpt = fsize / 2;
	else
		/*fsize is even*/
		endpt = fsize / 2 - 1;
	//printf("\nfno = %d, endpt = %d", fno, endpt);
	signal = (short int *) malloc(fsize * sizeof(short int));
	if (fno == 0) {
		for (i = (-fsize / 2), j = 0; i < 0; i++, j++)
			signal[j] = 0; /*prepad zeroes*/
		fseek(w.fp, (44 + (0 * (w.hdr.BitsPerSample / 8))), SEEK_SET);
		fread(&signal[j], sizeof(short), endpt, w.fp); /*read from file till the endpt*/
	} else {
		for (i = fsize - osize, j = 0; i < fsize; i++, j++)
			arr[j] = arr[i];
		fread(&signal[osize], sizeof(short), (fsize - osize), w.fp); /*read only new data from file*/
	}
	for (i = osize; i < fsize; i++)
		arr[i] = (double) ((Number) signal[i] / (Number) 32768);
	free(signal);
	fno++;
	return;
}
/********************************************************************/
/*WavFrameWrite:
 w: wav structure with file pointer initialed to write point in file,
 arr: array containing one frame of data.
 fsize: frame size in samples*/
/********************************************************************/
void WavFrameWrite(
		WAV w,
		double *arr,
		int fsize) {
	short *temp;
	int i;
	temp = (short *) malloc(fsize * sizeof(short));
	for (i = 0; i < fsize; i++) {
		if (arr[i] > 1.0)
			printf("\nwarning: data clipped while writing to wav file");
		temp[i] = (short) (arr[i] * 32768);
	}

	//fwrite args: 1) buffer 2)bytes per entity 3)no. of entities
	fwrite(temp, 2, fsize, w.fp);
	free(temp);
}
/********************************************************************/
/*n: number of samples written
 ch: 1: mono, 2: stereo
 sr: sampling rate
 bps: bits per sample
 */
/********************************************************************/
void WavHeaderInit(
		WAV *w,
		int n,
		int ch,
		int sr,
		int bps) {
	w->hdr.ChunkID[0] = 'R';
	w->hdr.ChunkID[1] = 'I';
	w->hdr.ChunkID[2] = 'F';
	w->hdr.ChunkID[3] = 'F';

	w->hdr.ChunkSize = 36 + n * 2;

	w->hdr.Format[0] = 'W';
	w->hdr.Format[1] = 'A';
	w->hdr.Format[2] = 'V';
	w->hdr.Format[3] = 'E';

	w->hdr.Subchunk1ID[0] = 'f';
	w->hdr.Subchunk1ID[1] = 'm';
	w->hdr.Subchunk1ID[2] = 't';
	w->hdr.Subchunk1ID[3] = ' ';

	w->hdr.Subchunk1Size = 16;
	w->hdr.AudioFormat = 1;
	w->hdr.NumChannels = ch;
	w->hdr.SampleRate = sr;		//Sampling Rate

	w->hdr.ByteRate = ch * sr * bps / 8;
	w->hdr.BlockAlign = ch * bps / 8;
	w->hdr.BitsPerSample = bps;

	w->hdr.Subchunk2ID[0] = 'd';
	w->hdr.Subchunk2ID[1] = 'a';
	w->hdr.Subchunk2ID[2] = 't';
	w->hdr.Subchunk2ID[3] = 'a';

	w->hdr.Subchunk2Size = n * ch * (bps / 8);		//Number of bytes in data
}
/********************************************************************/
/*WavHeaderWrite
 WavHeaderInit is expected to be called before writing*/
/********************************************************************/
void WavHeaderWrite(
		WAV w) {

	fseek(w.fp, 0, SEEK_SET);

	fwrite(w.hdr.ChunkID, 1, 4, w.fp);

	fwrite(&w.hdr.ChunkSize, 1, 4, w.fp);
	fwrite(w.hdr.Format, 1, 4, w.fp);
	fwrite(w.hdr.Subchunk1ID, 1, 4, w.fp);

	fwrite(&w.hdr.Subchunk1Size, 1, 4, w.fp);	//format length

	fwrite(&w.hdr.AudioFormat, 1, 2, w.fp);		//format pcm header
	fwrite(&w.hdr.NumChannels, 1, 2, w.fp);		//mono - stereo

	fwrite(&w.hdr.SampleRate, 1, 4, w.fp);			//Sample rate

	fwrite(&w.hdr.ByteRate, 1, 4, w.fp);			//bytes/sec

	fwrite(&w.hdr.BlockAlign, 1, 2, w.fp);			//Block align

	fwrite(&w.hdr.BitsPerSample, 1, 2, w.fp);		//bits/sample
	fwrite(&w.hdr.Subchunk2ID, 1, 4, w.fp);		//Data header

	fwrite(&w.hdr.Subchunk2Size, 1, 4, w.fp);		//the long int of data size..
}
/********************************************************************/
/* WavGetNumberOfFrames : returns the number of frames
 w: initialized WAV structure
 fsize: frame size
 osize: overlap size (not timestep)
 */
/********************************************************************/
int WavGetNumberOfFrames(
		WAV w,
		int fsize,
		int osize) {
	int noframes;
	FILE *fp;
	fp = fopen("frameinfo.txt", "wt");

	noframes = (int) (floor(((0.5 * w.hdr.Subchunk2Size) - fsize) / (fsize - osize)) + 1);

	fprintf(fp, "\nframesize   (in samples) : %d (in time) %f ms", fsize,
			(double) fsize * 1000 / (double) w.hdr.SampleRate);
	fprintf(fp, "\noverlapsize (in samples) : %d (in time) %f ms", osize,
			(double) osize * 1000 / (double) w.hdr.SampleRate);
	fprintf(fp, "\nNumber of Frames : %d", noframes);
	fclose(fp);
	return noframes;
}
/********************************************************************/
/*number of samples for a given time s*/
int WavToSamples(
		WAV w,
		double time) {
	return (time * w.hdr.SampleRate);
}

/********************************************************************/
/*time in ms corresponding to a given number if samples*/
double WavToTime(
		WAV w,
		int samples) {
	return (double) ((double) samples / (double) w.hdr.SampleRate); //time in ms.
}

/********************************************************************/
/*
 inputs: frame : frame number
 timestep : time step in ms
 framesize : framesize in ms

 returns: time instant corresponding to the frame center.
 */
/********************************************************************/
double FrameNumToTime(
		int frame,
		double timestep,
		double framesize) {
	return ((framesize / 2.0) + (frame - 1) * (double) timestep);
}
/********************************************************************/
void WavPrintHeader(
		WAV w) {
	printf("\nHeader info:\n");
	printf("\nChunkID       : %c%c%c%c", (w.hdr).ChunkID[0], (w.hdr).ChunkID[1], (w.hdr).ChunkID[2],
			(w.hdr).ChunkID[3]);
	printf("\nFileSize(-8b) : %d", w.hdr.ChunkSize);
	printf("\nFormat        : %c%c%c%c", w.hdr.Format[0], w.hdr.Format[1], w.hdr.Format[2], w.hdr.Format[3]);
	printf("\nFMT chunkID   : %c%c%c%c", w.hdr.Subchunk1ID[0], w.hdr.Subchunk1ID[1], w.hdr.Subchunk1ID[2],
			w.hdr.Subchunk1ID[3]);
	printf("\nSubchunk1Size : %d", w.hdr.Subchunk1Size);
	printf("\nAudioFormat   : %d", w.hdr.AudioFormat);
	printf("\nNumChannels   : %d", w.hdr.NumChannels);
	printf("\nSampleRate    : %d", w.hdr.SampleRate);
	printf("\nByteRate      : %d", w.hdr.ByteRate);
	printf("\nBlockAlign    : %d", w.hdr.BlockAlign);
	printf("\nBitsPerSample : %d", w.hdr.BitsPerSample);
	printf("\nDATA chunkID  : %c%c%c%c", w.hdr.Subchunk2ID[0], w.hdr.Subchunk2ID[1], w.hdr.Subchunk2ID[2],
			w.hdr.Subchunk2ID[3]);
	printf("\nData(in bytes): %d", w.hdr.Subchunk2Size);
}
/********************************************************************/
double WavGetMaxima(
		WAV w) {
	int i = 0, bps;
	short value, max = -1;
	bps = w.hdr.BitsPerSample / 8;
	while (!feof(w.fp)) {
		fseek(w.fp, (44 + i), 0);
		fread(&value, sizeof(short), 1, w.fp);
		/*testing*/
		//printf("\n%f",(double)value/(double)32768);
		if (value > max)
			max = value;
		i += bps;
	}
	return (double) max / (double) 32768;
}
/********************************************************************/
