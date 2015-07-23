/*
 * CyborgManager.h
 *
 *  Created on: Feb 3, 2015
 *      Author: sharath
 */

#ifndef CYBORGMANAGER_H_
#define CYBORGMANAGER_H_
#include "commonUtils.h"
#include "Dictionary.h"
#include "../mfcc/MfccManager.h"
#include "ViterbiAligner.h"
#include "CyborgResults.h"
#include "params.h"
namespace std {

class CyborgManager {
private:
	MfccManager mfccManager;
	ViterbiAligner aligner;
	string m_workingDir;

public:
	CyborgManager();
	CyborgManager(
			string amPath,
			string fstPath,
			string fstName,
			string workingDir);

	void setWorkingDirectory(
			string workingDir);
	void setAM(
			string amPath);
	void loadDictionary(
			string dictPath,
			string fillerDictPath);
	void createFST(
			string transcription,
			string fstPath,
			string fstName);
	void setFST(
			string fstPath,
			string fstName);
	void setNeuralNetworkMode(
			bool neuralNetMode);
	void setPosteriorMode(
				bool posteriorMode);
	void setBackTraceMode(
			int backTraceMode);
	void setBeamWidth(
			Number beamWidth);
	void setBackwardFA(
			bool backwardFA);
	void setDecoderMode(
			bool backwardFA);

	Number2DArray computeMFCC(
			string filePath,
			string fileName);
	Number2DArray computeMFCC(
			NumberArray audioIn);
	Number2DArray readMFCCFile(
			string filePath,
			string fileName);
	void getContextMFCCFrames(
			Number2DArrayRef mfccMatrix);
	void printMFCC(
			Number2DArrayRef feat);
	void writeMFCC(
			string fileName,
			Number2DArrayRef feat);

	void doAlignment(
			Number2DArrayRef mfcc);

	CyborgResultsContainer getWdseg() {
		return aligner.getWdseg();
	}
	;
	CyborgResultsContainer getPhseg() {
		return aligner.getPhseg();
	}
	;
	CyborgResultsContainer getStseg() {
		return aligner.getStseg();
	}
	;

	void printWdseg();
	void printPhseg();
	void printStseg();

	void writeWdseg(
			string fileName);
	void writePhseg(
			string fileName);
	void writeStseg(
			string fileName);

	Number getLyricScore(
			CyborgResultsContainerRef wdseg);
	virtual ~CyborgManager();
};

} /* namespace std */
#endif /* CYBORGMANAGER_H_ */
