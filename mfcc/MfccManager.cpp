/*
 * MfccManager.cpp
 *
 *  Created on: Jul 25, 2014
 *      Author: sharath
 */

#include "MfccManager.h"

namespace std {

void MfccManager::MfccManagerInit(
		Number winLen_sec,
		int Fs,
		Number lowerf,
		Number upperf,
		int numFilt,
		int nfft,
		int dctType) {

	m_Fs = Fs;
	m_winLen = (int) (winLen_sec * Fs);
	m_numCep = 13;
	mfcc.mfccInit(m_winLen, Fs, lowerf, upperf, numFilt, m_numCep, nfft, dctType);
	hammWindow.resize(m_winLen, 0.0);
	Dsp::hamming(hammWindow);
}


//Compute MFCC's of the audio file, specified by filePath and fileName
Number2DArray MfccManager::computeFeatures(
		string filePath,
		string fileName) {

	Number2DArray mfcc_13dim;
	// --------- Read audio file ---------//
	string refFilename = filePath + '\\' + fileName + ".wav";
	WAV_SIGNAL audioIn;
	WavLoadFile(refFilename.c_str(), &audioIn);

	if (m_Fs != audioIn.SampleRate) {
		cout << "ERROR: " << refFilename << " sampling frequency is " << audioIn.SampleRate << "Hz and not " << m_Fs
				<< "Hz" << endl;
		exit(0);
	}

	NumberArray filtAudio(audioIn.DataSize, 0);
	Dsp::prEmpFilter(audioIn.Data, filtAudio);
	int hopLen = 0.01 * m_Fs;

	// --------- Main Algo starts here ---------//
	int frame, newest;
	for (newest = 0, frame = 0; newest < (audioIn.DataSize - m_winLen); newest = newest + hopLen, frame++) //
			{

		NumberArray input(m_winLen, 0.0);
		for (int i = 0; i < m_winLen; i++) {
			input[i] = filtAudio[newest + i] * hammWindow[i];
		}
		NumberArray mfccValues = mfcc.calculate(input);
		mfcc_13dim.push_back(mfccValues);
	}

	Number2DArray mfcc_13_cmn;
	mfcc_13_cmn = cepstralMeanNormalize(mfcc_13dim);

	// Calculate Delta/Double-Delta features
	return computeDeltaFeatures(mfcc_13_cmn);
}


//Compute MFCC's of the audio vector specified by audioIn
Number2DArray MfccManager::computeFeatures(
		NumberArray audioIn) {

	Number2DArray mfcc_13dim;
	NumberArray filtAudio(audioIn.size(), 0);
	Dsp::prEmpFilter(audioIn, filtAudio);
	int hopLen = 0.01 * m_Fs;

	// --------- Main Algo starts here ---------//
	int frame, newest;
	for (newest = 0, frame = 0; newest < ((int) audioIn.size() - m_winLen); newest = newest + hopLen, frame++) //
			{

		NumberArray input(m_winLen, 0.0);
		for (int i = 0; i < m_winLen; i++) {
			input[i] = filtAudio[newest + i] * hammWindow[i];
		}
		NumberArray mfccValues = mfcc.calculate(input);
		mfcc_13dim.push_back(mfccValues);
	}

	Number2DArray mfcc_13_cmn;
	mfcc_13_cmn = cepstralMeanNormalize(mfcc_13dim);

	// Calculate Delta/Double-Delta features
	return computeDeltaFeatures(mfcc_13_cmn);
}


Number2DArray MfccManager::cepstralMeanNormalize(
		Number2DArrayRef feat_s) {

	int frames = feat_s.size();
	NumberArray mean(13, 0.0);
	// computing the mean for each of the 13 features for the current utterance
	for (int i = 0; i < 13; i++) {
		float avg = 0;
		for (int j = 0; j < frames; j++) {
			avg += feat_s[j][i];
		}
		mean[i] = avg / frames;
	}
	// subtracting the mean from each features
	for (int i = 0; i < 13; i++) {
		for (int j = 0; j < frames; j++) {
			feat_s[j][i] -= mean[i];
		}
	}
	return feat_s;

}

void MfccManager::printFeatures(
		Number2DArrayRef feat) {
	cout << endl;
	for (int i = 0; i < (int) feat.size(); i++) {
		for (int j = 0; j < (int) feat[i].size(); j++) {
			cout << feat[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void MfccManager::writeMFCC(
		string fileName,
		Number2DArrayRef feat) {

	ofstream wdfp(fileName.c_str());
	if (!wdfp.is_open()) {
		cout << "ERROR: " << fileName.c_str()
				<< " is in use by some other program. Please close the program before you run this code" << endl;
		exit(0);
	}
	for (int i = 0; i < (int) feat.size(); i++) {
		for (int j = 0; j < (int) feat[i].size(); j++) {
			wdfp << feat[i][j] << " ";
		}
		wdfp << endl;
	}
	wdfp << endl;
}

Number2DArray MfccManager::computeDeltaFeatures(
		Number2DArrayRef feat) {
	int frames = feat.size();
	Number2DArray feat_s = Number2DArray(frames + 6, NumberArray(13));

	// Appending 3 rows in the beginning and the end of feat.
	for (int i = 0, n = 3; i < frames; i++, n++)
		for (int j = 0; j < 13; j++)
			feat_s[n][j] = feat[i][j];

	// replicating the first 3 rows.
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 13; j++)
			feat_s[i][j] = feat_s[3][j];

	// replicating the last 3 rows.
	for (int i = frames + 3; i < frames + 6; i++)
		for (int j = 0; j < 13; j++)
			feat_s[i][j] = feat_s[3 + frames - 1][j];

	// computing the delta features
	Number2DArray feat_d = Number2DArray(frames, NumberArray(13));
	for (int i = 0, n = 3; i < frames; i++, n++)
		for (int j = 0; j < 13; j++)
			feat_d[i][j] = feat_s[n + 2][j] - feat_s[n - 2][j];

	// computing the double delta features
	Number2DArray feat_dd = Number2DArray(frames, NumberArray(13));
	for (int i = 0, n = 3; i < frames; i++, n++)
		for (int j = 0; j < 13; j++)
			feat_dd[i][j] = (feat_s[n + 3][j] - feat_s[n - 1][j]) - (feat_s[n + 1][j] - feat_s[n - 3][j]);

	// merging delta and double-delta with single features
	Number2DArray feat_s_d_dd = Number2DArray(frames, NumberArray(39, 0));
	for (int i = 0; i < frames; i++) {
		for (int j = 0; j < 13; j++) {
			feat_s_d_dd[i][j] = feat_s[i + 3][j];
			feat_s_d_dd[i][j + 13] = feat_d[i][j];
			feat_s_d_dd[i][j + 26] = feat_dd[i][j];
		}
	}
	// returning the array with delta and double delta features.
	return feat_s_d_dd;
}

Number2DArray MfccManager::readFeatFile(
		string featDir,
		string fileId) {
	Number2DArray featMat;
	ifstream featFile((featDir + "/" + fileId + ".mfc").c_str());
	string line;

	if (featFile.is_open()) {
		while (getline(featFile, line)) {
			NumberArray featVec;
			Number featVal;
			stringstream ss(line);
			while (ss >> featVal) {
				featVec.push_back(featVal);
			}
			featMat.push_back(featVec);
		}
		featFile.close();
	} else {
		cout << "ERROR: Failed to open " << (featDir + "/" + fileId + ".mfc").c_str() << endl;
		exit(0);
	}

	return featMat;
}

void MfccManager::getContextMFCCFrames(
		Number2DArrayRef mfccMatrix) {
	Number2DArray oldMFCCMatrix(mfccMatrix);
	// clear contents of mfccMatrix
	for (int i = 0; i < (int) mfccMatrix.size(); i++)
		mfccMatrix[i].clear();
	mfccMatrix.clear();

	// perform zero mean unit variance normalization
	Dsp::zeroMeanUnitVarianceNormalization(oldMFCCMatrix);

	// add context frames
	int numContextframes = 3;
	for (int i = numContextframes; i < ((int) oldMFCCMatrix.size() - numContextframes); i++) {
		NumberArray contextFrame;
		for (int j = (i - numContextframes); j < (i + numContextframes + 1); j++) {
			contextFrame.insert(contextFrame.end(), oldMFCCMatrix[j].begin(), oldMFCCMatrix[j].end());
		}
		mfccMatrix.push_back(contextFrame);
	}

	//Append first frame numContextframes times in the beginning
	NumberArray firstFrame = mfccMatrix[0];
	for (int i = 0; i < numContextframes; i++) {
		mfccMatrix.insert(mfccMatrix.begin(), firstFrame);
	}

	//Append last frame numContextframes times in the end
	NumberArray lastFrame = mfccMatrix[mfccMatrix.size() - 1];
	for (int i = 0; i < numContextframes; i++) {
		mfccMatrix.push_back(lastFrame);
	}

}

} /* namespace std */
