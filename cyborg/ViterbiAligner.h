/*
 * ViterbiAligner.h
 *
 *  Created on: Jul 24, 2014
 *      Author: sharath
 */

#ifndef VITERBIALIGNER_H_
#define VITERBIALIGNER_H_
#include <iomanip>
#include <algorithm>
#include "AcousticModel.h"
#include "Dictionary.h"
#include "commonUtils.h"
#include "Fst.h"
#include "FstLinkedList.h"
#include "LatticeLinkedList.h"
#include "sort.h"
#include "CyborgResults.h"
#include "../mfcc/Dsp.h"
namespace std {

class ViterbiAligner {

private:
	Number m_varFloor;
	Number m_mixwFloor;
	Number m_tmatFloor;

	Number m_logBase;
	Number m_beamWidth; // DONT USE THIS. USE m_am.m_beamWidth

	string m_workDir;
	string m_fileName;

	int m_backTraceMode;
	int m_nBest;

	AcousticModel m_am;
	Dictionary m_dict;
	FstLinkedList m_fstNetwork;
	CyborgResultsContainer m_stseg, m_phseg, m_wdseg;

	bool m_backwardFA;
	bool m_decoderFlag;

	bool m_neuralNetworkMode;

	ofstream posteriorFid;

public:
	ViterbiAligner() {
	}
	;

	ViterbiAligner(
			Number varFloor,
			Number mixwFloor,
			Number tmatFloor,
			int incompleteMode,
			Number logBase,
			int best,
			Number beamWidth);
	void initAM(
			string amPath);
	void initFST(
			string fsgPath,
			string fsgName);
	void setWorkDir(
			string workDir);
	void loadDictionary(
			string dictPath,
			string fillerDictPath);
	void setBackTraceMode(
			int backTraceMode);

	void setNeuralNetworkMode(
			bool neuralNetworkMode);
	void setBeamWidth(
			Number beamWidth);
	void setBackwardFA(
			bool backwardFA);
	void setDecoderMode(
			bool decoderFlag);

	void createFST(
			string transcription,
			string fstPath,
			string fstName);
	void viterbiAlignerInit(
			string dictPath,
			string fillerDictPath,
			string amPath,
			string audDir,
			string featDir,
			string phsegDir,
			NumberArrayRef floorVals,
			double beamWidth,
			int incompleteMode,
			bool fsgMode,
			string fsgName,
			bool calcFst,
			int nBestResults,
			string thresholdPath);
	void doAlignment(
			Number2DArrayRef mfcc);
	IntArray posteriorScores(
			NumberArrayRef featVec);
	int frameLikelihood(
			int state,
			NumberArrayRef x);
	int transitionLikelihood(
			int s,
			int from,
			int to);
	void calcPhoneFractions(
			vector<LatticeLinkedList*> & NodeListCurrentTimeFrame,
			Number2DArrayRef phoneFrac,
			int i);
	void mergeNodes(
			vector<LatticeLinkedList*> & parentListLevelNplus1,
			vector<LatticeLinkedList*> & parentListLevel);
	void deleteTree(
			vector<LatticeLinkedList*> & NodeListNextTimeFrame);
	int getNBestBacktraceStartNode(
			vector<LatticeLinkedList*> & NodeListCurrentTimeFrame,
			vector<LatticeLinkedList *> & maxCostNode);
	int getBacktraceStartNode(
			vector<LatticeLinkedList*> & NodeListCurrentTimeFrame,
			vector<LatticeLinkedList *> &maxCostNode);
	void backtrace_backwardFA(
			LatticeLinkedList* maxCostNode);
	void backtrace(
			LatticeLinkedList* maxCostNode);
	void reverseResultVector();
	void cleanResultVector();
	CyborgResultsContainer getStseg() {
		return m_stseg;
	}
	;
	CyborgResultsContainer getPhseg() {
		return m_phseg;
	}
	;
	CyborgResultsContainer getWdseg() {
		return m_wdseg;
	}
	;

	~ViterbiAligner() {
	}
	;

};

} /* namespace std */
#endif /* VITERBIALIGNER_H_ */
