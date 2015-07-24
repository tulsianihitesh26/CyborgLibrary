/*
 * CyborgManager.cpp
 *
 *  Created on: Feb 3, 2015
 *      Author: sharath
 */

#include "CyborgManager.h"

namespace std {

CyborgManager::CyborgManager() :
				aligner(VAR_FLOOR, MIXWT_FLOOR, TMAT_FLOOR, INCOMPLETE_MODE, LOGBASE, N_BEST, BEAM_WIDTH) {
	mfccManager.MfccManagerInit(WINLEN_SEC, SYS_FS, LOWER_F, UPPER_F, NUM_FILT, NFFT, DCT_TYPE);
}

CyborgManager::CyborgManager(
		string amPath,
		string fstPath,
		string fstName,
		string workingDir) :
				m_workingDir(workingDir) {
	aligner.initAM(amPath);
	aligner.initFST(fstPath, fstName);
	aligner.setWorkDir(workingDir);
}

void CyborgManager::setAM(
		string amPath) {
	aligner.initAM(amPath);
}

void CyborgManager::setPosteriorPath(
		string posteriorPath,
		string fileId) {
	aligner.setPosteriorPath(posteriorPath,fileId);
}

void CyborgManager::loadDictionary(
		string dictPath,
		string fillerDictPath) {
	aligner.loadDictionary(dictPath, fillerDictPath);
}

void CyborgManager::createFST(
		string transcription,
		string fstPath,
		string fstName) {
	aligner.createFST(transcription, fstPath, fstName);
}

void CyborgManager::setFST(
		string fstPath,
		string fstName) {
	aligner.initFST(fstPath, fstName);
}

void CyborgManager::setNeuralNetworkMode(
		bool neuralNetMode) {
	aligner.setNeuralNetworkMode(neuralNetMode);
}

void CyborgManager::setPosteriorMode(
		bool posteriorMode) {
	aligner.setPosteriorMode(posteriorMode);
}

void CyborgManager::setBackTraceMode(
		int backTraceMode) {
	aligner.setBackTraceMode(backTraceMode);
}

void CyborgManager::setBeamWidth(
		Number beamWidth) {
	aligner.setBeamWidth(beamWidth);
}

void CyborgManager::setDecoderMode(
		bool backwardFA) {
	aligner.setDecoderMode(backwardFA);
}

void CyborgManager::setBackwardFA(
		bool backwardFA) {
	aligner.setBackwardFA(backwardFA);
}

void CyborgManager::setWorkingDirectory(
		string workingDir) {
	aligner.setWorkDir(workingDir);
	m_workingDir = workingDir;
}

Number2DArray CyborgManager::computeMFCC(
		string filePath,
		string fileName) {
	return mfccManager.computeFeatures(filePath, fileName);
}

Number2DArray CyborgManager::computeMFCC(
		NumberArray audioIn) {
	return mfccManager.computeFeatures(audioIn);
}

Number2DArray CyborgManager::readMFCCFile(
		string filePath,
		string fileName) {
	return mfccManager.readFeatFile(filePath, fileName);
}

void CyborgManager::getContextMFCCFrames(
		Number2DArrayRef mfccMatrix) {
	mfccManager.getContextMFCCFrames(mfccMatrix);
}

void CyborgManager::printMFCC(
		Number2DArrayRef feat) {
	mfccManager.printFeatures(feat);
}

void CyborgManager::writeMFCC(
		string fileName,
		Number2DArrayRef feat) {

	mfccManager.writeMFCC(m_workingDir + "//" + fileName + ".mfc", feat);
}

void CyborgManager::doAlignment(
		Number2DArrayRef mfcc) {
	aligner.cleanResultVector();
	aligner.doAlignment(mfcc);
}

void CyborgManager::printWdseg() {
	cout << endl << endl;
	cout << "Start\tEnd\tPhone\tPost./Phone" << endl;
	CyborgResultsContainer wdseg = getWdseg();
	for (int rInd = 0; rInd < wdseg.size; rInd++) {
		cout << wdseg.frameStart[rInd] << "\t" << wdseg.frameEnd[rInd] << "\t" << wdseg.phone[rInd] << "\t"
				<< wdseg.posterior[rInd] << endl;
	}

}
void CyborgManager::printPhseg() {
	cout << endl << endl;
	cout << "Start\tEnd\tPhone\tPost./Phone" << endl;
	CyborgResultsContainer phseg = getPhseg();
	for (int rInd = 0; rInd < phseg.size; rInd++) {
		cout << phseg.frameStart[rInd] << "\t" << phseg.frameEnd[rInd] << "\t" << phseg.phone[rInd] << "\t"
				<< phseg.posterior[rInd] << endl;
	}

}
void CyborgManager::printStseg() {
	cout << endl << endl;
	cout
			<< "Frame\tPhone\tphoneState\tState\t Acoustic Scr\tTransition Scr\tNode Cost (AScr+TScr)\twordEndFlag\twordPosition\tword"
			<< endl;
	CyborgResultsContainer stseg = getStseg();
	for (int rInd = 0; rInd < stseg.size; rInd++) {
		cout << stseg.frameStart[rInd] << "\t" << stseg.phone[rInd] << "\t" << stseg.phoneState[rInd] << "\t"
				<< stseg.state[rInd] << "\t" << stseg.acousticScore[rInd] << "\t" << stseg.transitionScore[rInd] << "\t"
				<< stseg.score[rInd] << "\t" << stseg.wordEnd[rInd] << "\t" << stseg.position[rInd] << "\t"
				<< stseg.word[rInd] << endl;
	}

}

void CyborgManager::writeWdseg(
		string fileName) {
	ofstream wdfp((m_workingDir + "//" + fileName + ".wdseg").c_str());
	if (!wdfp.is_open()) {
		cout << "ERROR: " << (m_workingDir + "//" + fileName + ".wdseg").c_str()
				<< " is in use by some other program. Please close the program before you run this code" << endl;
		exit(0);
	}

	wdfp << "Start\tEnd\tPhone\tPost./Phone" << endl;
	CyborgResultsContainer wdseg = getWdseg();
	for (int rInd = 0; rInd < wdseg.size; rInd++) {
		wdfp << wdseg.frameStart[rInd] << "\t" << wdseg.frameEnd[rInd] << "\t" << wdseg.phone[rInd] << "\t"
				<< wdseg.posterior[rInd] << endl;
	}

}
void CyborgManager::writePhseg(
		string fileName) {
	// ---- Print phseg --- //
	ofstream phfp((m_workingDir + "//" + fileName + ".phseg").c_str());
	if (!phfp.is_open()) {
		cout << "ERROR: " << (m_workingDir + "//" + fileName + ".phseg").c_str()
				<< " is in use by some other program. Please close the program before you run this code" << endl;
		exit(0);
	}
	phfp << "Start\tEnd\tPhone\tPost./Phone" << endl;
	CyborgResultsContainer phseg = getPhseg();
	for (int rInd = 0; rInd < phseg.size; rInd++) {
		phfp << phseg.frameStart[rInd] << "\t" << phseg.frameEnd[rInd] << "\t" << phseg.phone[rInd] << "\t"
				<< phseg.posterior[rInd] << endl;
	}

}
void CyborgManager::writeStseg(
		string fileName) {
	ofstream stfp((m_workingDir + "//" + fileName + ".stseg").c_str());
	if (!stfp.is_open()) {
		cout << "ERROR: " << (m_workingDir + "//" + fileName + ".stseg").c_str()
				<< " is in use by some other program. Please close the program before you run this code" << endl;
		exit(0);
	}

	stfp
			<< "Frame\tPhone\tPhone-State\tState\t Acoustic Scr\tTransition Scr\tNode Cost (AScr+TScr)\twordEndFlag\twordPosition\tword"
			<< endl;
	CyborgResultsContainer stseg = getStseg();
	for (int rInd = 0; rInd < stseg.size; rInd++) {
		stfp << stseg.frameStart[rInd] << "\t" << stseg.phone[rInd] << "\t" << stseg.phoneState[rInd] << "\t"
				<< stseg.state[rInd] << "\t" << stseg.acousticScore[rInd] << "\t" << stseg.transitionScore[rInd] << "\t"
				<< stseg.score[rInd] << "\t" << stseg.wordEnd[rInd] << "\t" << stseg.position[rInd] << "\t"
				<< stseg.word[rInd] << endl;
	}
}

Number CyborgManager::getLyricScore(
		CyborgResultsContainerRef wdseg) {

	int goodFAThreshold = -20000;
	int partiallyGoodFAThreshold = -25000;
	int finalFAScore = 0;
	for (int i = 0; i < wdseg.size; i++) {
		if (wdseg.phone[i].compare("SIL") == 0) {
			continue;
		}

		int vSegScore = wdseg.posterior[i];
		if (vSegScore > goodFAThreshold) //good
				{
			finalFAScore += wdseg.frameEnd[i] - wdseg.frameStart[i];
		} else if ((partiallyGoodFAThreshold <= vSegScore) && (vSegScore < goodFAThreshold)) //partially good 0.5
				{
			finalFAScore += 0.5 * (wdseg.frameEnd[i] - wdseg.frameStart[i]);
		}
//		cout << wdseg.phone[i] <<"\t" <<vSegScore<<endl;
	}

	return finalFAScore / (Number) wdseg.frameEnd[wdseg.size - 1];
}

CyborgManager::~CyborgManager() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
