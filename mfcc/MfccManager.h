/*
 * MfccManager.h
 *
 *  Created on: Jul 25, 2014
 *      Author: sharath
 */

#ifndef MFCCMANAGER_H_
#define MFCCMANAGER_H_
#include <iostream>
#include <vector>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include "wav.h"
#include "Dsp.h"
#include "Mfcc.h"
#include "../cyborg/commonUtils.h"
namespace std {

class MfccManager {

private:
	int m_Fs;
	int m_winLen;
	int m_numCep;
	Mfcc mfcc;
	NumberArray hammWindow;
public:
	MfccManager() {
	}
	;
	void MfccManagerInit(
			Number winLen_sec,
			int Fs,
			Number lowerf,
			Number upperf,
			int numFilt,
			int nfft,
			int dctType);
	Number2DArray computeFeatures(
			string filePath,
			string fileName);
	Number2DArray computeFeatures(
			NumberArray audioIn);
	Number2DArray cepstralMeanNormalize(
			Number2DArrayRef feat_s);
	Number2DArray computeDeltaFeatures(
			Number2DArrayRef);
	Number2DArray readFeatFile(
			string featDir,
			string fileId);
	Number2DArray readPostFeatFile(
				string featDir,
				string fileId);
	void printFeatures(
			Number2DArrayRef feat);
	void writeMFCC(
			string fileName,
			Number2DArrayRef feat);
	void getContextMFCCFrames(
			Number2DArrayRef mfccMatrix);
	~MfccManager() {
	}
	;
};

} /* namespace std */
#endif /* MFCCMANAGER_H_ */
