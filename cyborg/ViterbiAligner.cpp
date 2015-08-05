/*
 * ViterbiAligner.cpp
 *
 *  Created on: Jul 24, 2014
 *      Author: sharath
 */

#include "ViterbiAligner.h"
#define PROFILE_ENABLED (0)
#define SENONE_PROFILE_ENABLED (0)

namespace std {

ViterbiAligner::ViterbiAligner(
		Number varFloor,
		Number mixwFloor,
		Number tmatFloor,
		int incompleteMode,
		Number logBase,
		int nBest,
		Number beamWidth) :
				m_varFloor(varFloor),
				m_mixwFloor(mixwFloor),
				m_tmatFloor(tmatFloor),
				m_logBase(logBase),
				m_beamWidth(beamWidth),
				m_backTraceMode(incompleteMode),
				m_nBest(nBest),
				m_backwardFA(false),
				m_neuralNetworkMode(false) {
}

void ViterbiAligner::initAM(
		string amPath) {
	if (m_neuralNetworkMode) {
		m_am.acousticModelInit_NN(amPath, m_varFloor, m_mixwFloor, m_tmatFloor, m_logBase, m_beamWidth);
	} else {
		m_am.acousticModelInit(amPath, m_varFloor, m_mixwFloor, m_tmatFloor, m_logBase, m_beamWidth);
	}
}

void ViterbiAligner::setPosteriorPath(
		string posteriorPath,
		string fileId) {
		m_posteriorPath = posteriorPath + "/" + fileId + ".post";
		cout << m_posteriorPath<<endl;
}

void ViterbiAligner::initFST(
		string fstPath,
		string fstName) {
	m_fstNetwork.initFstLinkedList(fstPath, fstName);
}

void ViterbiAligner::setWorkDir(
		string workDir) {
	m_workDir = workDir;
}

void ViterbiAligner::loadDictionary(
		string dictPath,
		string fillerDictPath) {
	m_dict.dictionaryInit(dictPath, fillerDictPath, m_backwardFA);
}

void ViterbiAligner::setBackTraceMode(
		int backTraceMode) {
	m_backTraceMode = backTraceMode;
}

void ViterbiAligner::setNeuralNetworkMode(
		bool neuralNetworkMode) {
	m_neuralNetworkMode = neuralNetworkMode;
}

void ViterbiAligner::setPosteriorMode(
		bool posteriorMode) {
	m_posteriorMode = posteriorMode;
}

void ViterbiAligner::setBeamWidth(
		Number beamWidth) {
	m_am.setBeamWidth(beamWidth);
}

void ViterbiAligner::setBackwardFA(
		bool backwardFA) {
	m_backwardFA = backwardFA;
}

void ViterbiAligner::setDecoderMode(
		bool decoderFlag) {
	m_decoderFlag = decoderFlag;
}

void ViterbiAligner::createFST(
		string transcription,
		string fstPath,
		string fstName) {
	Fst fst(fstPath, fstName, transcription, m_dict, m_backwardFA, m_decoderFlag);
}

// Viterbi tree construction using FST knowledge and likelihood calculation for every MFCC frame.
// User required alignment is backtraced from the tree and results are written in cyborgResults container.
void ViterbiAligner::doAlignment(
		Number2DArrayRef mfcc) {
	//posteriorFid.open((m_workDir+"//"+m_fileName+"_posteriors.txt").c_str());

// Create first node of Viterbi tree
	IntArray senoneScores;
	LatticeLinkedList *root = new LatticeLinkedList();
	root->costOfNode = 0;
	root->parent = NULL;
	root->branchId = 0;
	vector<Triphone> rootPhoneList = m_fstNetwork.getNextPhoneFst(0);
	Triphone rootPhone = rootPhoneList[0];
	root->triPhone = rootPhone.getTriphone();
	root->stateIds = m_am.getStates(root->triPhone + "\t-\t-\t-");
	root->linkedListStateNo = 0;
	root->nextlinkedListStateNo = m_fstNetwork.getNextPhoneFstIndex(root->triPhone, 0);
	root->curStateNo = 0;
	root->curPhoneStateNo = root->stateIds[root->curStateNo + 2];
	root->position = rootPhone.position;

	// Initialize for forward/backward alignment
	int TotalFrames = mfcc.size();
	int from = 1, to = TotalFrames, step = 1;
	// Set posterior array to read
	Number2DArray posterior;
	if (m_backwardFA) {
		//if (m_posteriorMode){
			// Read all Posteriors first and then use as needed
		//	m_am.nn.readPosteriors(m_posteriorPath,posterior);
		//	senoneScores = posteriorScores(posterior[TotalFrames - 1]);
		//} else {
			senoneScores = posteriorScores(mfcc[TotalFrames - 1]);
		//}
		root->observationId = TotalFrames - 1;
		from = TotalFrames - 2;
		to = -1;
		step = -1;
	} else {
		//if (m_posteriorMode){
			//	Read all Posteriors first and then use as needed
		//	m_am.nn.readPosteriors(m_posteriorPath,posterior);
		//	senoneScores = posteriorScores(posterior[0]);
		//} else {
					senoneScores = posteriorScores(mfcc[0]);
		//}
		//senoneScores = posteriorScores(mfcc[0]);
		root->observationId = 0;
	}

	root->languageScore = rootPhone.triphoneWeight;
	root->acousticScore = senoneScores[root->stateIds[root->curStateNo + 2]];
	root->transitionScore = transitionLikelihood(root->stateIds[1], root->curStateNo, root->curStateNo);
	root->costOfNode = root->acousticScore + root->transitionScore;
	root->costOfPath = root->costOfNode;
	root->costOfPathIncludingLw = root->costOfNode + root->languageScore;

	// The pointers NodeListCurrentTimeFrame, NodeListNextTimeFrame keep track of consecutive levels in Viterbi tree.
	vector<LatticeLinkedList*> NodeListCurrentTimeFrame;
	vector<LatticeLinkedList*> NodeListNextTimeFrame;
	NodeListCurrentTimeFrame.push_back(root);
//	Number2DArray phoneFrac(TotalFrames, NumberArray(am.phoneCnt, 0));

// Start building Viterbi tree
	for (int i = from; i != to; i += step) {
		//if (m_posteriorMode){
		//	senoneScores = posteriorScores(posterior[i]);
		//} else {
			senoneScores = posteriorScores(mfcc[i]);
		//}
		for (int j = 0; j < (int) NodeListCurrentTimeFrame.size(); j++) {
			LatticeLinkedList * curTreeNode = NodeListCurrentTimeFrame[j];

			// Each phone has 3 states - 0, 1 and 2. As long as the parent state in the Viterbi tree is not equal to 2,
			// the children nodes will be of the same phone. If the parent is in state 2, then the children will be in state 2
			// and the state 0 of the next phone from FST.

			if (curTreeNode->curStateNo != 2) { // Children nodes of same phone
				LatticeLinkedList *newTreeNode1 = new LatticeLinkedList();
				newTreeNode1->copyFormAnotherNode(*curTreeNode);
				newTreeNode1->observationId = i;
				newTreeNode1->parent = curTreeNode;
				newTreeNode1->languageScore = curTreeNode->languageScore;
				newTreeNode1->acousticScore = senoneScores[newTreeNode1->stateIds[newTreeNode1->curStateNo + 2]];
				newTreeNode1->transitionScore = transitionLikelihood(curTreeNode->stateIds[1], curTreeNode->curStateNo,
						newTreeNode1->curStateNo);
				newTreeNode1->costOfNode = newTreeNode1->acousticScore + newTreeNode1->transitionScore;
				newTreeNode1->costOfPath = curTreeNode->costOfPath + newTreeNode1->costOfNode;
				newTreeNode1->costOfPathIncludingLw = curTreeNode->costOfPathIncludingLw + newTreeNode1->costOfNode
						+ newTreeNode1->languageScore;
				newTreeNode1->position = curTreeNode->position;
				newTreeNode1->newWord = 0;
				newTreeNode1->word = "-";
				newTreeNode1->wordEnd = 0;
				newTreeNode1->lastSIL = curTreeNode->lastSIL;
				newTreeNode1->lastPhone = curTreeNode->lastPhone;

				NodeListNextTimeFrame.push_back(newTreeNode1);
				curTreeNode->numChild++;

				LatticeLinkedList * newTreeNode2 = new LatticeLinkedList();
				newTreeNode2->copyFormAnotherNode(*curTreeNode);
				newTreeNode2->observationId = i;
				newTreeNode2->curStateNo = curTreeNode->curStateNo + 1;
				newTreeNode2->curPhoneStateNo = curTreeNode->curPhoneStateNo + 1;
				newTreeNode2->parent = curTreeNode;
				newTreeNode2->languageScore = curTreeNode->languageScore;
				newTreeNode2->acousticScore = senoneScores[newTreeNode2->stateIds[newTreeNode2->curStateNo + 2]];
				newTreeNode2->transitionScore = transitionLikelihood(curTreeNode->stateIds[1], curTreeNode->curStateNo,
						newTreeNode2->curStateNo);
				newTreeNode2->costOfNode = newTreeNode2->acousticScore + newTreeNode2->transitionScore;
				newTreeNode2->costOfPath = curTreeNode->costOfPath + newTreeNode2->costOfNode;
				newTreeNode2->costOfPathIncludingLw = curTreeNode->costOfPathIncludingLw + newTreeNode2->costOfNode
						+ newTreeNode2->languageScore;
				newTreeNode2->position = curTreeNode->position;
				newTreeNode2->newWord = 0;
				newTreeNode2->word = "-";
				newTreeNode2->wordEnd = 0;
				newTreeNode2->lastSIL = curTreeNode->lastSIL;
				newTreeNode2->lastPhone = curTreeNode->lastPhone;

				NodeListNextTimeFrame.push_back(newTreeNode2);
				curTreeNode->numChild++;
			} else { // Children nodes of current and next phone.
				LatticeLinkedList *newTreeNode1 = new LatticeLinkedList();
				newTreeNode1->copyFormAnotherNode(*curTreeNode);
				newTreeNode1->observationId = i;
				newTreeNode1->parent = curTreeNode;
				newTreeNode1->languageScore = curTreeNode->languageScore;
				newTreeNode1->acousticScore = senoneScores[newTreeNode1->stateIds[newTreeNode1->curStateNo + 2]];
				newTreeNode1->transitionScore = transitionLikelihood(curTreeNode->stateIds[1], curTreeNode->curStateNo,
						newTreeNode1->curStateNo);
				newTreeNode1->costOfNode = newTreeNode1->acousticScore + newTreeNode1->transitionScore;
				newTreeNode1->costOfPath = curTreeNode->costOfPath + newTreeNode1->costOfNode;
				newTreeNode1->costOfPathIncludingLw = curTreeNode->costOfPathIncludingLw + newTreeNode1->costOfNode
						+ newTreeNode1->languageScore;
				newTreeNode1->position = curTreeNode->position;
				newTreeNode1->newWord = 0;
				newTreeNode1->word = "-";
				newTreeNode1->wordEnd = 0;
				newTreeNode1->lastSIL = curTreeNode->lastSIL;
				newTreeNode1->lastPhone = curTreeNode->lastPhone;

				NodeListNextTimeFrame.push_back(newTreeNode1);
				curTreeNode->numChild++;

				if (m_fstNetwork.lastPhoneReached(curTreeNode->nextlinkedListStateNo)) {
					vector<Triphone> nextTriphone = m_fstNetwork.getNextPhoneFst(curTreeNode->nextlinkedListStateNo);
					for (int triPhCnt = 0; triPhCnt < (int) nextTriphone.size(); triPhCnt++) {
						Triphone phone = nextTriphone[triPhCnt];
						if (!phone.getTriphone().compare("") == 0) // if <eps> comes before the last state, then this will return a null string, which has to be ignored.
								{
							LatticeLinkedList *newTreeNode2 = new LatticeLinkedList();
							newTreeNode2->triPhone = phone.getTriphone();
							newTreeNode2->observationId = i;
							newTreeNode2->stateIds = m_am.getStates(newTreeNode2->triPhone + "\t-\t-\t-");

							if (newTreeNode2->stateIds[0] == -1) {
								newTreeNode2->triPhone = m_fstNetwork.getNextPhoneFst(
										curTreeNode->nextlinkedListStateNo)[0].triphone;
								newTreeNode2->stateIds = m_am.getStates(newTreeNode2->triPhone + "\t-\t-\t-");
							}

							newTreeNode2->linkedListStateNo = phone.getlinkedListStateNo();
							newTreeNode2->nextlinkedListStateNo = phone.getNextLinkedListStateNo();
							newTreeNode2->curStateNo = 0;
							newTreeNode2->curPhoneStateNo = newTreeNode2->stateIds[newTreeNode2->curStateNo + 2];

							newTreeNode2->parent = curTreeNode;

							if (nextTriphone[triPhCnt].prevWordFlag == 1) { // handling cases, when word definition is present in #0 state
								newTreeNode2->parent->word = nextTriphone[triPhCnt].prevWord;
							}

							newTreeNode2->languageScore = nextTriphone[triPhCnt].triphoneWeight;
							newTreeNode2->acousticScore = senoneScores[newTreeNode2->stateIds[newTreeNode2->curStateNo
									+ 2]];
							newTreeNode2->transitionScore = transitionLikelihood(curTreeNode->stateIds[1],
									curTreeNode->curStateNo, 3);
							newTreeNode2->costOfNode = newTreeNode2->acousticScore + newTreeNode2->transitionScore;
							newTreeNode2->costOfPath = curTreeNode->costOfPath + newTreeNode2->costOfNode;
							newTreeNode2->costOfPathIncludingLw = curTreeNode->costOfPathIncludingLw
									+ newTreeNode2->costOfNode + newTreeNode2->languageScore;
							newTreeNode2->position = nextTriphone[triPhCnt].position;
							newTreeNode2->newWord = nextTriphone[triPhCnt].newWord;
							newTreeNode2->word = nextTriphone[triPhCnt].word;
							newTreeNode2->wordEnd = nextTriphone[triPhCnt].wordEnd;
							newTreeNode2->lastSIL = nextTriphone[triPhCnt].lastSIL;
							newTreeNode2->lastPhone = nextTriphone[triPhCnt].lastPhone;

							NodeListNextTimeFrame.push_back(newTreeNode2);
							curTreeNode->numChild++;
						}
					}
				}
			}
		}

		// Perform pruning
		mergeNodes(NodeListNextTimeFrame, NodeListCurrentTimeFrame);
//		calcPhoneFractions(NodeListCurrentTimeFrame, phoneFrac, i);
	}

	//------------------- Start backtracing from the best phone ----------------------------//

	vector<LatticeLinkedList*> maxCostNodeArr;
	//	cout << " Last SIL with Best Likelihood " <<maxCostNode->costOfPath << " " << maxCostNode->stateIds[maxCostNode->curStateNo + 2] << " " <<maxCostNode->triPhone<<endl;;

	// N-best results
	if (m_nBest > 1) {

		//find out the N-best nodes to back trace from
		if (getNBestBacktraceStartNode(NodeListCurrentTimeFrame, maxCostNodeArr) == 0) {
			for (int i = 0; i < (int) maxCostNodeArr.size(); i++) {
				backtrace(maxCostNodeArr[i]); // start back tracing
			}
		}
	} else {
		//based on user requirements chose the node to backtrace from
		if (getBacktraceStartNode(NodeListCurrentTimeFrame, maxCostNodeArr) == 0) {
			maxCostNodeArr.push_back(maxCostNodeArr[0]);
			// start back tracing
			if (m_backwardFA) {
				backtrace_backwardFA(maxCostNodeArr[0]); //reverse FA
			} else {
				backtrace(maxCostNodeArr[0]);
			}
		}
	}

	deleteTree(NodeListCurrentTimeFrame); // delete Viterbi tree
	posteriorFid.close();
}

void ViterbiAligner::deleteTree(
		vector<LatticeLinkedList*> & childNodes) {

	vector<LatticeLinkedList*> parentNodes;
	for (int i = 0; i < (int) childNodes.size(); i++) {
		LatticeLinkedList *currentNode = childNodes[i];
		if (currentNode->parent != NULL) {
			currentNode->parent->numChild--;
			if (currentNode->parent->numChild < 0) {
				cout << "WARNING: Backtracking should not have reached here" << endl;
				return;
			}
			if (currentNode->parent->numChild == 0) {
				parentNodes.push_back(currentNode->parent);
			}
//			cout << currentNode->triPhone << " "  ;
			delete currentNode;

		} else {
			delete currentNode;
			return;
		}
	}
	if (parentNodes.size() > 0) {
//		cout<<endl;
		deleteTree(parentNodes);
	}
}

IntArray ViterbiAligner::posteriorScores(
		NumberArrayRef featVec) {
	IntArray senoneScores;
	if (m_neuralNetworkMode) {
		if (m_posteriorMode){
			senoneScores.resize(featVec.size(), 0.0);
			for (int i = 0; i < (int) featVec.size(); i++) {
				senoneScores[i] = (int) ((featVec[i] == 0.0) ? LOGPROB_ZERO : AcousticModel::myLog(featVec[i]));
				//cout << senoneScores[i] << " ";
			}
			//cout << "\n";
		} else {
			// NN POSTERIORS
			NumberArray posteriors;
			m_am.nn.classify(featVec, posteriors);
			senoneScores.resize(posteriors.size(), 0.0);
			for (int i = 0; i < (int) posteriors.size(); i++) {
				//posteriorFid<< posteriors[i]<< " ";
				senoneScores[i] = (int) ((posteriors[i] == 0.0) ? LOGPROB_ZERO : AcousticModel::myLog(posteriors[i]));

			}

			//posteriorFid << endl;
		}
	} else {
		// GMM POSTERIORS
		senoneScores.resize(m_am.SENONES, 0.0);
		int maxSenoneScore = -99999999;
		for (int i = 0; i < m_am.SENONES; i++) {
			senoneScores[i] = frameLikelihood(i, featVec);
			//		cout << "senone " << i << " score " <<  senoneScores[i]<<endl;
			if (senoneScores[i] > maxSenoneScore) {
				maxSenoneScore = senoneScores[i];
			}
		}

		for (int i = 0; i < m_am.SENONES; i++) {
			senoneScores[i] -= maxSenoneScore;
//			posteriorFid<< senoneScores[i]<< " ";

			//		cout << senoneScores[i] << " ";
		}
//		posteriorFid << endl;

		//	cout << endl;
	}
	return senoneScores;
}

void ViterbiAligner::mergeNodes(
		vector<LatticeLinkedList*> & parentListLevelNplus1,
		vector<LatticeLinkedList*> & parentListLevel) {
	int maxCost = parentListLevelNplus1[0]->costOfPathIncludingLw;

	/* Handle similar node cases
	 	         sil,3
	 			  /\
	 	   sil,3     sil,4
	       /\          /\
	  sil,3  sil,4 sil,4 sil,5

	  One of the two sil,4 nodes has to be kept before heading to next level.
	  The sil,4 with lowest path score is saved, and the other is deleted
	 */

	for (int i = 0; i < (int) parentListLevelNplus1.size(); i++) {
		LatticeLinkedList *currentNode = parentListLevelNplus1[i];
		if (currentNode->active) {
			if (currentNode->costOfPathIncludingLw > maxCost) {
				maxCost = currentNode->costOfPathIncludingLw;
			}
			for (int j = i + 1; j < (int) parentListLevelNplus1.size(); j++) {
				LatticeLinkedList *compareNode = parentListLevelNplus1[j];
				if (compareNode->active) {
					if ((currentNode->linkedListStateNo == compareNode->linkedListStateNo)
							&& (currentNode->stateIds[currentNode->curStateNo + 2]
									== compareNode->stateIds[compareNode->curStateNo + 2])) {
						if (currentNode->costOfPathIncludingLw > compareNode->costOfPathIncludingLw) {
							compareNode->parent->numChild--;
							compareNode->active = false;
						} else {
							currentNode->parent->numChild--;
							currentNode->active = false;
							break;
						}
					}
				}
			}
		}
	}

	/* Handle nodes which die before the last level is reached

		   sil,3       sil,4        sil,5
	       /\          / x           /\
	   sil,3  sil,4  (sil,4)    sil,5   a,21
	   /\     / x       x x       / x    x \

	  (sil,4) node which abruptly ends, and all its parents who have zero nodes have to be deleted.
	  The single x cases are handled by the loop after the following loop
	 */
	for (int i = 0; i < (int) parentListLevel.size(); i++) {
		if (parentListLevel[i]->numChild == 0) {
			vector<LatticeLinkedList*> newList;
			newList.push_back(parentListLevel[i]);
			deleteTree(newList);
		}
	}

	parentListLevel.clear();
	int pruneLimit = m_am.m_beamWidth;
	int beam = maxCost + pruneLimit;
	for (int i = 0; i < (int) parentListLevelNplus1.size(); i++) {
		LatticeLinkedList *compareNode = parentListLevelNplus1[i];

		if (compareNode->active) {
			if (compareNode->costOfPathIncludingLw > beam) {
				parentListLevel.push_back(compareNode);
			}
		} else {
			delete compareNode;
		}
	}
	parentListLevelNplus1.clear();
}

int ViterbiAligner::frameLikelihood(
		int state,
		NumberArrayRef x) {
	int gaussian = m_am.GAUSSIAN;
	NumberArray mean, var;
	Number dval, diff, gauscr;
	int score = LOGPROB_ZERO;

	for (int i = 0; i < gaussian; i++) {
		dval = m_am.PRECOMP[state][i];
		for (int j = 0; j < m_am.m_cepLen; j++) {
			diff = x[j] - m_am.MEAN[state][i][j];
			dval -= diff * diff * m_am.VAR[state][i][j];
		}
		gauscr = (m_am.m_logBaseDenominator * dval) + m_am.MIXWT[state][i];
		if (gauscr > score) {
			score = gauscr;
		}
	}
	return score;
}

/**
 * This function returns the transition probability.
 * @param s - state number
 * @param from - transition from
 * @param to - transition to
 *
 */
int ViterbiAligner::transitionLikelihood(
		int s,
		int from,
		int to) {
	int temp;
	//This condition is if the transition is from one HMM
	//to another HMM
	if ((from != to) && (to % 3 == 0))
		temp = m_am.TMAT[s][2][3];
	else
		// this condition is for transition from one state to another
		// of same HMM
		temp = m_am.TMAT[s][from % 3][to % 3];
	return temp;
}

int ViterbiAligner::getNBestBacktraceStartNode(
		vector<LatticeLinkedList*> & NodeListCurrentTimeFrame,
		vector<LatticeLinkedList *> & maxCostNode) {
	//------------------- Backtracing----------------------------//
	// Search for best state in the leaves - two possibilities
	//		a) m_backTraceMode =1 : choose best likelihood - dont bother what alignment
	//		b) m_backTraceMode =0 : choose leaf only if reached SIL
	//		c) m_backTraceMode =2 : choose best likelihood - Force align to given transcript

	if (m_backTraceMode == 1) // N best likelihood
			{
		IntArray silScores, oldSilScores, silInd;

		for (int i = 0; i < (int) NodeListCurrentTimeFrame.size(); i++) {
			silScores.push_back(NodeListCurrentTimeFrame[i]->costOfPathIncludingLw);
			silInd.push_back(i);
		}

		oldSilScores = silScores;
		std::vector<size_t> ind;
		sort(silScores, silScores, ind);

		int vecSize = (int) silScores.size();
		int thEnd = ((vecSize - m_nBest) < 0) ? 0 : (vecSize - m_nBest);
		for (int j = vecSize - 1; j >= thEnd; j--) {
			printf("b[%d] = %d = a[i[%d]] = a[%d] = %d\n", j, silScores[j], j, ind[j], oldSilScores[ind[j]]);
			maxCostNode.push_back(NodeListCurrentTimeFrame[silInd[ind[j]]]);
		}
	} else if (m_backTraceMode == 0) { // N best ending with silence

		IntArray states = m_am.getStates("SIL\t-\t-\t-");
		IntArray silScores, oldSilScores, silInd;

		for (int i = 0; i < (int) NodeListCurrentTimeFrame.size(); i++) {
			if ((NodeListCurrentTimeFrame[i]->lastSIL == 1)
					&& (NodeListCurrentTimeFrame[i]->stateIds[NodeListCurrentTimeFrame[i]->curStateNo + 2] == states[4])) // b)FORCE TO ALIGN LAST PHONE WITH SIL
					{
				silScores.push_back(NodeListCurrentTimeFrame[i]->costOfPathIncludingLw);
				silInd.push_back(i);
			}
		}

		oldSilScores = silScores;
		std::vector<size_t> ind;
		sort(silScores, silScores, ind);
		int vecSize = (int) silScores.size();
		int thEnd = ((vecSize - m_nBest) < 0) ? 0 : (vecSize - m_nBest);
		for (int j = vecSize - 1; j >= thEnd; j--) {
			printf("b[%d] = %d = a[i[%d]] = a[%d] = %d\n", j, silScores[j], j, ind[j], oldSilScores[ind[j]]);
			maxCostNode.push_back(NodeListCurrentTimeFrame[silInd[ind[j]]]);
		}

		if (maxCostNode.size() == 0) {
			cout << "ERROR: Alignment failed. Did not reach silence at the end" << endl;
			return -1;
		}
	} else { // choose best likelihood - Force align to given transcript

		IntArray states = m_am.getStates("SIL\t-\t-\t-");
		IntArray silScores, oldSilScores, silInd;

		for (int i = 0; i < (int) NodeListCurrentTimeFrame.size(); i++) {
			if ((NodeListCurrentTimeFrame[i]->lastPhone == 1) && (NodeListCurrentTimeFrame[i]->curStateNo == 2)) // b)FORCE TO ALIGN LAST PHONE WITH SIL
					{
				silScores.push_back(NodeListCurrentTimeFrame[i]->costOfPathIncludingLw);
				silInd.push_back(i);
			}
		}

		oldSilScores = silScores;
		std::vector<size_t> ind;
		sort(silScores, silScores, ind);
		int vecSize = (int) silScores.size();
		int thEnd = ((vecSize - m_nBest) < 0) ? 0 : (vecSize - m_nBest);
		for (int j = vecSize - 1; j >= thEnd; j--) {
			maxCostNode.push_back(NodeListCurrentTimeFrame[silInd[ind[j]]]);
		}

		if (maxCostNode.size() == 0) {
			cout << "ERROR: Alignment failed. Did not reach silence at the end" << endl;
			return -1;
		}
	}
	return 0;
}

int ViterbiAligner::getBacktraceStartNode(
		vector<LatticeLinkedList*> & NodeListCurrentTimeFrame,
		vector<LatticeLinkedList *> & maxCostNode) {
	//------------------- Backtracing----------------------------//
	// Search for best state in the leaves - two possibilities
	//		a) m_incompleteMode =1 : choose best likelihood - dont bother what alignment
	//		b) m_incompleteMode =0 : choose leaf only if reached SIL
	//		c) m_incompleteMode =2 : choose best likelihood - Force align to given transcript

	if (m_backTraceMode == 1) { // a) BEST LIKELIHOOD

		LatticeLinkedList *node = NodeListCurrentTimeFrame[0];
		maxCostNode.push_back(node);
		int maxCost = node->costOfPathIncludingLw;

		for (int i = 1; i < (int) NodeListCurrentTimeFrame.size(); i++) {
			node = NodeListCurrentTimeFrame[i];
			if (maxCost < node->costOfPathIncludingLw) {
				maxCost = node->costOfPathIncludingLw;
				maxCostNode[0] = node;
			}
		}
	} else if (m_backTraceMode == 0) { //b) FORCE TO ALIGN LAST PHONE WITH SIL

		int maxCost = -999999999;
		IntArray states = m_am.getStates("SIL\t-\t-\t-");
		maxCostNode.resize(1);
		for (int i = 0; i < (int) NodeListCurrentTimeFrame.size(); i++) {
			if ((NodeListCurrentTimeFrame[i]->lastSIL == 1)
					&& (NodeListCurrentTimeFrame[i]->stateIds[NodeListCurrentTimeFrame[i]->curStateNo + 2] == states[4])) // b)FORCE TO ALIGN LAST PHONE WITH SIL
					{
				if (maxCost < NodeListCurrentTimeFrame[i]->costOfPathIncludingLw) {
					maxCost = NodeListCurrentTimeFrame[i]->costOfPathIncludingLw;
					maxCostNode[0] = NodeListCurrentTimeFrame[i];
				}
			}
		}

		if (maxCostNode[0] == NULL) {
			cout << "ERROR: Alignment failed. Did not reach silence at the end" << endl;
			return -1;

		}
	} else { //choose best likelihood - Force align to given transcript

		int maxCost = -999999999;
		LatticeLinkedList *node = NULL;
		maxCostNode.resize(1);

		for (int i = 0; i < (int) NodeListCurrentTimeFrame.size(); i++) {
			node = NodeListCurrentTimeFrame[i];
			if ((node->lastPhone == 1) && (node->curStateNo == 2)) {
				if (maxCost < node->costOfPathIncludingLw) {
					maxCost = node->costOfPathIncludingLw;
					maxCostNode[0] = node;
				}
			}
		}
		if (maxCostNode[0] == NULL) {
			cout << "ERROR: Alignment failed. Did not reach silence at the end" << endl;
			return -1;
		}
	}
	return 0;
}

void ViterbiAligner::cleanResultVector() {

	m_wdseg.phone.clear();
	m_wdseg.word.clear();
	m_wdseg.frameStart.clear();
	m_wdseg.frameEnd.clear();
	m_wdseg.acousticScore.clear();
	m_wdseg.transitionScore.clear();
	m_wdseg.score.clear();
	m_wdseg.posterior.clear();
	m_wdseg.state.clear();
	m_wdseg.position.clear();
	m_wdseg.wordEnd.clear();
	m_wdseg.size = 0;

	m_stseg.phone.clear();
	m_stseg.word.clear();
	m_stseg.frameStart.clear();
	m_stseg.frameEnd.clear();
	m_stseg.acousticScore.clear();
	m_stseg.transitionScore.clear();
	m_stseg.score.clear();
	m_stseg.posterior.clear();
	m_stseg.state.clear();
	m_stseg.position.clear();
	m_stseg.wordEnd.clear();
	m_stseg.size = 0;

	m_phseg.phone.clear();
	m_phseg.word.clear();
	m_phseg.frameStart.clear();
	m_phseg.frameEnd.clear();
	m_phseg.acousticScore.clear();
	m_phseg.transitionScore.clear();
	m_phseg.score.clear();
	m_phseg.posterior.clear();
	m_phseg.state.clear();
	m_phseg.position.clear();
	m_phseg.wordEnd.clear();
	m_phseg.size = 0;

}

void ViterbiAligner::reverseResultVector() {

	std::reverse(m_wdseg.phone.begin(), m_wdseg.phone.end());
	std::reverse(m_wdseg.word.begin(), m_wdseg.word.end());
	std::reverse(m_wdseg.frameStart.begin(), m_wdseg.frameStart.end());
	std::reverse(m_wdseg.frameEnd.begin(), m_wdseg.frameEnd.end());
	std::reverse(m_wdseg.acousticScore.begin(), m_wdseg.acousticScore.end());
	std::reverse(m_wdseg.transitionScore.begin(), m_wdseg.transitionScore.end());
	std::reverse(m_wdseg.score.begin(), m_wdseg.score.end());
	std::reverse(m_wdseg.posterior.begin(), m_wdseg.posterior.end());
	std::reverse(m_wdseg.state.begin(), m_wdseg.state.end());
	std::reverse(m_wdseg.position.begin(), m_wdseg.position.end());
	std::reverse(m_wdseg.wordEnd.begin(), m_wdseg.wordEnd.end());
	m_wdseg.size = m_wdseg.frameStart.size();

	std::reverse(m_phseg.phone.begin(), m_phseg.phone.end());
	std::reverse(m_phseg.word.begin(), m_phseg.word.end());
	std::reverse(m_phseg.frameStart.begin(), m_phseg.frameStart.end());
	std::reverse(m_phseg.frameEnd.begin(), m_phseg.frameEnd.end());
	std::reverse(m_phseg.acousticScore.begin(), m_phseg.acousticScore.end());
	std::reverse(m_phseg.transitionScore.begin(), m_phseg.transitionScore.end());
	std::reverse(m_phseg.score.begin(), m_phseg.score.end());
	std::reverse(m_phseg.posterior.begin(), m_phseg.posterior.end());
	std::reverse(m_phseg.state.begin(), m_phseg.state.end());
	std::reverse(m_phseg.position.begin(), m_phseg.position.end());
	std::reverse(m_phseg.wordEnd.begin(), m_phseg.wordEnd.end());
	m_phseg.size = m_phseg.frameStart.size();

	std::reverse(m_stseg.phone.begin(), m_stseg.phone.end());
	std::reverse(m_stseg.word.begin(), m_stseg.word.end());
	std::reverse(m_stseg.frameStart.begin(), m_stseg.frameStart.end());
	std::reverse(m_stseg.frameEnd.begin(), m_stseg.frameEnd.end());
	std::reverse(m_stseg.acousticScore.begin(), m_stseg.acousticScore.end());
	std::reverse(m_stseg.transitionScore.begin(), m_stseg.transitionScore.end());
	std::reverse(m_stseg.score.begin(), m_stseg.score.end());
	std::reverse(m_stseg.posterior.begin(), m_stseg.posterior.end());
	std::reverse(m_stseg.state.begin(), m_stseg.state.end());
	std::reverse(m_stseg.phoneState.begin(), m_stseg.phoneState.end());
	std::reverse(m_stseg.position.begin(), m_stseg.position.end());
	std::reverse(m_stseg.wordEnd.begin(), m_stseg.wordEnd.end());
	m_stseg.size = m_stseg.frameStart.size();

}

void ViterbiAligner::backtrace(
		LatticeLinkedList* maxCostNode) {
	// Start collecting word, phone and state segmentations

	int endFrame = maxCostNode->observationId;
	int endCost = maxCostNode->costOfPath;

	int wordEndFrame = maxCostNode->observationId;
	int wordEndCost = maxCostNode->costOfPath;
	bool wordFoundFlag = false;
	bool skipFlag = false;
	string wordFound = "";

	int costOfPhone, normCostOfPhone;
	while (maxCostNode->parent != NULL) {
		if ((maxCostNode->curStateNo == 0) && (maxCostNode->parent->curStateNo == 2)) {
			costOfPhone = endCost - maxCostNode->parent->costOfPath;
			normCostOfPhone = costOfPhone / (Number) (endFrame - maxCostNode->parent->observationId);
			m_phseg.frameStart.push_back(maxCostNode->parent->observationId + 1);
			m_phseg.frameEnd.push_back(endFrame);
			m_phseg.score.push_back(costOfPhone);
			m_phseg.phone.push_back(maxCostNode->triPhone);
			m_phseg.posterior.push_back(normCostOfPhone);

			endFrame = maxCostNode->parent->observationId;
			endCost = maxCostNode->parent->costOfPath;
		}

		if (maxCostNode->wordEnd) { // Handle Wordseg

			costOfPhone = wordEndCost - maxCostNode->parent->costOfPath;
			m_wdseg.frameStart.push_back(maxCostNode->parent->observationId + 1);
			m_wdseg.frameEnd.push_back(wordEndFrame);
			m_wdseg.score.push_back(costOfPhone);
			m_wdseg.posterior.push_back(costOfPhone / (wordEndFrame - maxCostNode->parent->observationId));

			if (maxCostNode->word.compare("-") != 0) {
				m_wdseg.phone.push_back(maxCostNode->word);
				skipFlag = true;
			} else if (wordFoundFlag == true) {
				m_wdseg.phone.push_back(wordFound);
			} else {
				m_wdseg.phone.push_back("SIL");
			}
			wordFoundFlag = false;
			wordEndFrame = maxCostNode->parent->observationId;
			wordEndCost = maxCostNode->parent->costOfPath;

		}

		if ((maxCostNode->word.compare("-") != 0) && (wordFoundFlag == false) && (skipFlag == false)) // this is done coz the word definition may not align with #1
				{ // it was observed if two words "char - c aa r" and "che - c e" were present in the vocab, the word definition of char is placed on "aa" and for che on "e".
				  //This is because, "c" state is shared by char and che. and hence this check had to be done
			wordFound = maxCostNode->word;
			wordFoundFlag = true;
		}
		skipFlag = false;
		m_stseg.frameStart.push_back(maxCostNode->observationId);
		m_stseg.phone.push_back(maxCostNode->triPhone);
		m_stseg.state.push_back(maxCostNode->curStateNo);
		m_stseg.phoneState.push_back(maxCostNode->curPhoneStateNo);
		m_stseg.acousticScore.push_back(maxCostNode->acousticScore);
		m_stseg.transitionScore.push_back(maxCostNode->transitionScore);
		m_stseg.score.push_back(maxCostNode->costOfNode);
		m_stseg.position.push_back(maxCostNode->position);
		m_stseg.wordEnd.push_back(maxCostNode->wordEnd);
		m_stseg.word.push_back(maxCostNode->word);

		maxCostNode = maxCostNode->parent;
	}

	// Reached first state
	m_wdseg.frameStart.push_back(maxCostNode->observationId);
	m_wdseg.frameEnd.push_back(wordEndFrame);
	m_wdseg.score.push_back(wordEndCost);
	m_wdseg.phone.push_back("SIL");
	m_wdseg.posterior.push_back(wordEndCost / (wordEndFrame - maxCostNode->observationId + 1));

	m_stseg.frameStart.push_back(maxCostNode->observationId);
	m_stseg.phone.push_back(maxCostNode->triPhone);
	m_stseg.state.push_back(maxCostNode->curStateNo);
	m_stseg.phoneState.push_back(maxCostNode->curPhoneStateNo);
	m_stseg.acousticScore.push_back(maxCostNode->acousticScore);
	m_stseg.transitionScore.push_back(maxCostNode->transitionScore);
	m_stseg.score.push_back(maxCostNode->costOfNode);
	m_stseg.position.push_back(maxCostNode->position);
	m_stseg.wordEnd.push_back(maxCostNode->wordEnd);
	m_stseg.word.push_back(maxCostNode->word);

	m_phseg.frameStart.push_back(maxCostNode->observationId);
	m_phseg.frameEnd.push_back(endFrame);
	m_phseg.score.push_back(endCost);
	m_phseg.phone.push_back(maxCostNode->triPhone);
	m_phseg.posterior.push_back(endCost / (Number) (endFrame - maxCostNode->observationId));

	reverseResultVector();
}

void ViterbiAligner::backtrace_backwardFA(
		LatticeLinkedList* maxCostNode) {
	// Start collecting word, phone and state segmentations

	int endFrame = maxCostNode->observationId;
	int endCost = maxCostNode->costOfPath;

	int wordEndFrame = maxCostNode->observationId;
	int wordEndCost = maxCostNode->costOfPath;
	bool wordFoundFlag = false;
	bool skipFlag = false;
	string wordFound = "";

	int costOfPhone, normCostOfPhone;
	while (maxCostNode->parent != NULL) {
		if ((maxCostNode->curStateNo == 0) && (maxCostNode->parent->curStateNo == 2)) {
			costOfPhone = endCost - maxCostNode->parent->costOfPath;
			normCostOfPhone = costOfPhone / (Number) (maxCostNode->parent->observationId - endFrame);
			m_phseg.frameEnd.push_back(maxCostNode->parent->observationId - 1);
			m_phseg.frameStart.push_back(endFrame);
			m_phseg.score.push_back(costOfPhone);
			m_phseg.phone.push_back(maxCostNode->triPhone);
			m_phseg.posterior.push_back(normCostOfPhone);

			endFrame = maxCostNode->parent->observationId;
			endCost = maxCostNode->parent->costOfPath;
		}

		if (maxCostNode->wordEnd) { // Handle Wordseg

			costOfPhone = wordEndCost - maxCostNode->parent->costOfPath;
			m_wdseg.frameEnd.push_back(maxCostNode->parent->observationId - 1);
			m_wdseg.frameStart.push_back(wordEndFrame);
			m_wdseg.score.push_back(costOfPhone);
			m_wdseg.posterior.push_back(costOfPhone / (maxCostNode->parent->observationId - wordEndFrame));

			if (maxCostNode->word.compare("-") != 0) {
				m_wdseg.phone.push_back(maxCostNode->word);
				skipFlag = true;
			} else if (wordFoundFlag == true) {
				m_wdseg.phone.push_back(wordFound);
			} else {
				m_wdseg.phone.push_back("SIL");
			}
			wordFoundFlag = false;
			wordEndFrame = maxCostNode->parent->observationId;
			wordEndCost = maxCostNode->parent->costOfPath;

		}

		if ((maxCostNode->word.compare("-") != 0) && (wordFoundFlag == false) && (skipFlag == false)) // this is done coz the word definition may not align with #1
				{ // it was observed if two words "char - c aa r" and "che - c e" were present in the vocab, the word definition of char is placed on "aa" and for che on "e".
				  //This is because, "c" state is shared by char and che. and hence this check had to be done
			wordFound = maxCostNode->word;
			wordFoundFlag = true;
		}
		skipFlag = false;
		m_stseg.frameStart.push_back(maxCostNode->observationId);
		m_stseg.phone.push_back(maxCostNode->triPhone);
		m_stseg.state.push_back(maxCostNode->curStateNo);
		m_stseg.phoneState.push_back(maxCostNode->curPhoneStateNo);
		m_stseg.acousticScore.push_back(maxCostNode->acousticScore);
		m_stseg.transitionScore.push_back(maxCostNode->transitionScore);
		m_stseg.score.push_back(maxCostNode->costOfNode);
		m_stseg.position.push_back(maxCostNode->position);
		m_stseg.wordEnd.push_back(maxCostNode->wordEnd);
		m_stseg.word.push_back(maxCostNode->word);

		maxCostNode = maxCostNode->parent;
	}

	// Reached first state
	m_wdseg.frameEnd.push_back(maxCostNode->observationId);
	m_wdseg.frameStart.push_back(wordEndFrame);
	m_wdseg.score.push_back(wordEndCost);
	m_wdseg.phone.push_back("SIL");
	m_wdseg.posterior.push_back(wordEndCost / (maxCostNode->observationId - wordEndFrame));

	m_stseg.frameStart.push_back(maxCostNode->observationId);
	m_stseg.phone.push_back(maxCostNode->triPhone);
	m_stseg.state.push_back(maxCostNode->curStateNo);
	m_stseg.phoneState.push_back(maxCostNode->curPhoneStateNo);
	m_stseg.acousticScore.push_back(maxCostNode->acousticScore);
	m_stseg.transitionScore.push_back(maxCostNode->transitionScore);
	m_stseg.score.push_back(maxCostNode->costOfNode);
	m_stseg.position.push_back(maxCostNode->position);
	m_stseg.wordEnd.push_back(maxCostNode->wordEnd);
	m_stseg.word.push_back(maxCostNode->word);

	m_phseg.frameEnd.push_back(maxCostNode->observationId);
	m_phseg.frameStart.push_back(endFrame);
	m_phseg.score.push_back(endCost);
	m_phseg.phone.push_back(maxCostNode->triPhone);
	m_phseg.posterior.push_back(endCost / (Number) (maxCostNode->observationId - endFrame));

	m_wdseg.size = m_wdseg.frameStart.size();
	m_phseg.size = m_phseg.frameStart.size();
	m_stseg.size = m_stseg.frameStart.size();

}



void ViterbiAligner::calcPhoneFractions(
		vector<LatticeLinkedList*> & NodeListCurrentTimeFrame,
		Number2DArrayRef phoneFrac,
		int frame) {
	IntArray NodeListCurrentTimeFrame_phoneFrac(m_am.phoneCnt, 0);
//	int silCnt=0;
	for (int i = 0; i < (int) NodeListCurrentTimeFrame.size(); i++) {
		LatticeLinkedList *currentNode = NodeListCurrentTimeFrame[i];
		for (int j = 0; j < m_am.phoneCnt; j++) {
			if (currentNode->triPhone.compare(m_am.phoneList[j]) == 0) {
//				if(currentNode->triPhone.compare("SIL")!=0)
				{
					NodeListCurrentTimeFrame_phoneFrac[j]++;
				}
//				else
//				{
//					silCnt++;
//				}

			}
		}

	}
//	cout<<endl;
//	for(int j = 0; j < am.phoneCnt; j++)
//	{
//		cout <<NodeListCurrentTimeFrame_phoneFrac[j] << " ";
//	}
//	cout<<endl;

	for (int j = 0; j < m_am.phoneCnt; j++) {
//		if ((NodeListCurrentTimeFrame.size()-silCnt) > 0 )
//		{
//			phoneFrac[frame][j]=NodeListCurrentTimeFrame_phoneFrac[j]/(Number)(NodeListCurrentTimeFrame.size()-silCnt);
//		}
		phoneFrac[frame][j] = NodeListCurrentTimeFrame_phoneFrac[j] / (Number) NodeListCurrentTimeFrame.size();
//		cout << phoneFrac[frame][j] << " ";
	}
//	cout<<endl;
}

} /* namespace std */
