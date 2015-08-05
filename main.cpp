/*
 * main.cpp
 *
 *  Created on: Jul 23, 2014
 *      Author: sharath
 */

#include <iostream>
#include "cyborg/commonUtils.h"
#include "cyborg/CyborgManager.h"
#include "cyborg/CyborgResults.h"
using namespace std;
void printLog() {
	cout << " ERROR : Not enough input arguments " << endl;
	cout << "Arguments list definition for Aligner :" << endl;
	cout << "[NAME]                   [DESCRIPTION]             " << endl;
	cout << "-hmm            Input model files directory" << endl;
	cout << "-ctl            Control file listing utterances to be processed" << endl;
	cout << "-insent         Input transcriptions file corresponds to audio" << endl;
	cout << "-dict           Main pronunciation dictionary (lexicon) input file" << endl;
	cout << "-fdict          Filler dictionary input file" << endl;
	cout << "-wav            Input audio wav file(s) directory" << endl;
	cout
			<< "-fstdir         Path where fst files exit, If -fst = 0 then the directory should contain fileName.fst for each of the files in -ctl"
			<< endl;
	cout << "-fst            Flag to tell the decoder whether to calculate fst or not, " << endl;
	cout << "                    0- does not calculate, uses the fst in the fileid.fst format, " << endl;
	cout << "                    1- calculates FST, " << endl;
	cout
			<< "                    2- does not calculate fst, uses the global fst - expects an input for -fsg variable.					"
			<< endl;
	cout
			<< "-fsg            Finite state grammar mode. Name of FSG file in -fstdir path to be provided. When -fsg is provided -fst, -dict, -fdict is ignored."
			<< endl;
	cout << "-btmode         Backtrace mode. " << endl;
	cout << "                    0- Backtrace from last silence phone" << endl;
	cout << "                    1- Backtrace from the best likelihood node- dont bother about complete alignment"
			<< endl;
	cout << "                    2-Force align to the given transcript." << endl;
	cout << "-fa              Force alignment mode" << endl;
	cout << "                    0 - Forward force align" << endl;
	cout << "                    1 - Backward force align" << endl;
	cout << "                    2 - Limited vocabulary decoding mode" << endl;
	cout << "-phseg         Output directory for phone segmentation files" << endl;
	cout << "-nn			   NN-HMM mode" << endl;
	cout << "                    0- runs in GMM-HMM mode. Uses posteriors from GMM for recognition." << endl;
	cout << "                    1- runs in NN-HMM mode. Uses posteriors from NN for recognition." << endl;

	exit(0);
}

int main(
		int argc,
		char *argv[]) {

	string filePath;
	string transcriptPath;
	string fstDir;

	string amPath;
	string dictPath;
	string fillerDictPath;

	string audDir;
	string phsegDir;
	string posteriorDir;

	int fstFlag = 1;
	string globalFstName; // if fstFlag==2, then default fst to be used for all
	bool globalFstFlag = false;
	bool calcFst = true;
	bool calcMfcc = true;
	bool backwardFA = false;
	bool freeDecodeFlag = false;

	bool neuralNetMode = false;
	bool posteriorMode = false;

	int backTraceMode = 2; //		a) backTraceMode =0 : choose leaf only if reached SIL
						   //		b) backTraceMode =1 : choose best likelihood - dont bother what alignment
						   //		c) backTraceMode =2 : choose best likelihood - Force align to given transcript

	Number beamWidth = 1e-38; //Beam selecting active HMMs (relative to  best) in each frame [0(widest)..1(narrowest)]
							  //For widest beam, we cant have 0, because log(0) =Infinity. The lowest float value supported 1e-38, so we use it.

	CyborgManager myDecoder;
	int minArgCnt = 0;
	for (int i = 0; i < argc; i++) {
		if (strcmp(argv[i], "-hmm") == 0) {
			amPath = argv[++i];
			minArgCnt++;
		}
		if (strcmp(argv[i], "-ctl") == 0) {
			filePath = argv[++i];
			minArgCnt++;
		}
		if (strcmp(argv[i], "-wav") == 0) {
			audDir = argv[++i];
			minArgCnt++;
		}
		if (strcmp(argv[i], "-insent") == 0) {
			transcriptPath = argv[++i];
			minArgCnt++;
		}
		if (strcmp(argv[i], "-dict") == 0) {
			dictPath = argv[++i];
			minArgCnt++;
		}
		if (strcmp(argv[i], "-fdict") == 0) {
			fillerDictPath = argv[++i];
			minArgCnt++;
		}
		if (strcmp(argv[i], "-phseg") == 0) {
			phsegDir = argv[++i];
			minArgCnt++;
		}
		if (strcmp(argv[i], "-post") == 0) {
			posteriorDir = argv[++i];
		}
		if (strcmp(argv[i], "-fst") == 0) {
			fstFlag = atoi(argv[++i]);
			if (fstFlag == 1)
				calcFst = true;
			else if (fstFlag == 0)
				calcFst = false;
			else if (fstFlag == 2) {
				calcFst = false;
				globalFstFlag = true;
			} else {
				cout << "ERROR : -fst value can only be 0, 1 or 2. You gave " << argv[++i] << endl;
				exit(0);
			}

		}
		if (strcmp(argv[i], "-fsg") == 0) {
			globalFstName = argv[++i];
			globalFstFlag = true;
		}
		if (strcmp(argv[i], "-fstdir") == 0) {
			fstDir = argv[++i];
			minArgCnt++;
		}
		if (strcmp(argv[i], "-feat") == 0) {
			int featCalc = atoi(argv[++i]);
			if (featCalc == 1)
				calcMfcc = true;
			else if (featCalc == 0)
				calcMfcc = false;
			else {
				cout << "ERROR : -feat value can only be 0 or 1. You gave " << argv[++i] << endl;
				exit(0);
			}
		}
		if (strcmp(argv[i], "-btmode") == 0) {
			backTraceMode = atoi(argv[++i]);
		}
		if (strcmp(argv[i], "-beam") == 0) {
			beamWidth = atof(argv[++i]);
		}
		if (strcmp(argv[i], "-fa") == 0) {
			int faFlag = atoi(argv[++i]);

			if (faFlag == 0) {
				backwardFA = false;
				freeDecodeFlag = false;
			} else if (faFlag == 1) {
				backwardFA = true;
				freeDecodeFlag = false;
			} else {
				backwardFA = false;
				freeDecodeFlag = true;
			}
		}
		if (strcmp(argv[i], "-nn") == 0) {
			int nnFlag = atoi(argv[++i]);
			if (nnFlag == 1) {
				neuralNetMode = true;
				posteriorMode = false;
			} else if (nnFlag == 2){
				neuralNetMode = true;
				posteriorMode = true;
			} else {
				neuralNetMode = false;
				posteriorMode = false;
			}
		}
		if (strcmp(argv[i], "-help") == 0) {
			printLog();
		}
	}

	if ((fstFlag == 2) & (globalFstFlag == false)) {
		cout << "ERROR: -fst set to 2 but -fsg not mentioned." << endl;
		exit(0);
	}

	myDecoder.setNeuralNetworkMode(neuralNetMode);
	myDecoder.setBackTraceMode(backTraceMode);
	myDecoder.setAM(amPath);
	myDecoder.setBeamWidth(beamWidth);
	myDecoder.setBackwardFA(backwardFA);
	myDecoder.setDecoderMode(freeDecodeFlag);
	myDecoder.setPosteriorMode(posteriorMode);

	myDecoder.setWorkingDirectory(phsegDir);
	if (globalFstFlag) {
		calcFst = false;
		myDecoder.setFST(fstDir, globalFstName);
	}

	if (calcFst) {
		myDecoder.loadDictionary(dictPath, fillerDictPath);
	}

	// LOAD FILEID FILE
	ifstream fileIdFile(filePath.c_str());
	if (!fileIdFile.is_open()) {
		cout << "ERROR: Unable to open : " << filePath << endl;
		exit(0);
	}

	// LOAD TRANSCRIPTION FILE
	ifstream transcriptFile(transcriptPath.c_str());
	if (!transcriptFile.is_open()) {
		cout << "ERROR: Unable to open : " << transcriptPath << endl;
		exit(0);
	}

	string fileId;
	string transcriptLine;
	int lineNo = 0;
	while (getline(fileIdFile, fileId)) { // Iterate through filenames in fileid file
		fileId = Dictionary::string2words(fileId)[0];
		lineNo++;
		string transcriptFileLine;
		StringArray transcriptLine;
		if (getline(transcriptFile, transcriptFileLine)) {// Iterate through transcriptions in transcription file
			transcriptLine = Dictionary::getFileNameAndTranscript(transcriptFileLine);

			//Check if fileName specified in transcription matches with the fileName in fileid file.
			if (Dictionary::getFolderFileName(fileId).back().compare(transcriptLine[1]) == 0) {
				cout << endl << endl << transcriptLine[1] << " " << transcriptLine[0] << endl;

				// Calculate FST if necessary
				if (calcFst) {
					myDecoder.createFST(transcriptLine[0], fstDir, fileId);
				}
				if (globalFstFlag == false) { // When each audio file is given a precomputed FST
					myDecoder.setFST(fstDir, fileId);
				}

				// Calcualte MFCC if necessary
				Number2DArray mfccMatrix;
				if (posteriorMode){
					mfccMatrix = myDecoder.readPosteriorFile(posteriorDir, fileId);
				} else {
					if (calcMfcc) {
						mfccMatrix = myDecoder.computeMFCC(audDir, fileId);
					} else {
						mfccMatrix = myDecoder.readMFCCFile(audDir, fileId);
					}
				}

				// For neural network mode, add context frames
				if (neuralNetMode) {
					if (posteriorMode){
						myDecoder.setPosteriorPath(posteriorDir,fileId);
					} else {
						myDecoder.getContextMFCCFrames(mfccMatrix);
					}
				}


				// Perform alignment using Viterbi tree
				myDecoder.doAlignment(mfccMatrix);

				// Collect the state, phone and word segmentations
				CyborgResultsContainer wdseg = myDecoder.getWdseg();
				myDecoder.writeWdseg(transcriptLine[1]);
				myDecoder.writePhseg(transcriptLine[1]);
				myDecoder.writeStseg(transcriptLine[1]);
				myDecoder.printPhseg();
//				myDecoder.writeMFCC(transcriptLine[1], mfccMatrix);
			} else {
				cout << "ERROR: transcription and fileID mismatch at line " << lineNo << "!!!" << endl;
				cout << "FileName :  " << fileId << " transcript FileName  : " << transcriptLine[1] << endl;
				exit(0);
			}
		}
	}
	return 0;
}
