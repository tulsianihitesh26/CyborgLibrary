/*
 * FstLinkedList.cpp
 *
 *  Created on: Jul 24, 2014
 *      Author: sharath
 */

#include "FstLinkedList.h"

namespace std {

void FstLinkedList::initFstLinkedList(
		string fstPath,
		string fileName) {
	listAllFstNodes.clear();
	finalFstLinkedList.clear();
	string line;
	StringArray words;
	ifstream fileIdFile((fstPath + "//" + fileName + ".fst").c_str());
	if (!fileIdFile.is_open()) {
		cout << "FST file missing : " << (fstPath + "//" + fileName + ".fst").c_str() << endl;
		exit(0);
	}
	IntArray lastPhoneArr;

	while (getline(fileIdFile, line)) {
		words = Dictionary::string2words(line);
		FstNode curFstNode;
		if (words.size() == 5) //FST contains weights
				{
			curFstNode.start = atoi(words[0].c_str());
			curFstNode.end = atoi(words[1].c_str());
			curFstNode.in = words[2];
			curFstNode.op = words[3];
			curFstNode.wt = AcousticModel::myLog(atof(words[4].c_str()));
			listAllFstNodes.push_back(curFstNode);
//			cout <<curFstNode.start << " " << curFstNode.end<<" " << curFstNode.in<<" " << curFstNode.op <<endl;
		} else if (words.size() == 4) //FST does not contain weights
				{
			curFstNode.start = atoi(words[0].c_str());
			curFstNode.end = atoi(words[1].c_str());
			curFstNode.in = words[2];
			curFstNode.op = words[3];
			listAllFstNodes.push_back(curFstNode);
//			cout <<curFstNode.start << " " << curFstNode.end<<" " << curFstNode.in<<" " << curFstNode.op <<endl;
		} else if (words.size() == 2) // last state has weight
				{
			curFstNode.start = atoi(words[0].c_str());
			curFstNode.end = -1;
			curFstNode.wt = atof(words[1].c_str());
			listAllFstNodes.push_back(curFstNode);
			lastPhoneArr.push_back(curFstNode.start);
//			cout <<curFstNode.start << " " << curFstNode.end<<endl;
		} else if (words.size() == 1) // last state has no weight
				{
			curFstNode.start = atoi(words[0].c_str());
			curFstNode.end = -1;
			listAllFstNodes.push_back(curFstNode);
			lastPhoneArr.push_back(curFstNode.start);
//			cout <<curFstNode.start << " " << curFstNode.end<<endl;
		} else {
			cout << "ERROR : FST network should not have reached here : " << line << endl;
			exit(0);
		}

	}

//	cout << listAllFstNodes.size() <<endl;

	// if the lastPhoneArr points to #0, then make it point the previous phone
	for (int i = 0; i < (int) listAllFstNodes.size(); i++) {
		FstNode *curFstNode = &listAllFstNodes[i];

		for (int j = 0; j < (int) lastPhoneArr.size(); j++) {
			if (curFstNode->end == lastPhoneArr[j]) {
				if (curFstNode->in.compare("#0") == 0) {
					lastPhoneArr[j] = curFstNode->start;
				}
			}
		}
	}

	int currentState = 0;
	finalFstLinkedList.push_back(listAllFstNodes[0]);
	for (int i = 1; i < (int) listAllFstNodes.size(); i++) {
		FstNode *curFstNode = &listAllFstNodes[i];

		for (int j = 0; j < (int) lastPhoneArr.size(); j++) {
			if (curFstNode->end == lastPhoneArr[j]) {
				curFstNode->lastPhone = 1;
				if (curFstNode->in.compare("SIL") == 0) {
					curFstNode->lastSIL = 1;
				}
			}
		}

		if (curFstNode->start > currentState) {
			finalFstLinkedList.push_back(*curFstNode);
			currentState = curFstNode->start;

		} else {
			FstNode *endOfLinkedList = &finalFstLinkedList[finalFstLinkedList.size() - 1];
			while (endOfLinkedList->nextFstNode != NULL) {
				endOfLinkedList = endOfLinkedList->nextFstNode;
			}
			endOfLinkedList->nextFstNode = curFstNode;

		}
//		cout <<curFstNode->start << " " << curFstNode->end<<" " << curFstNode->in<<" " << curFstNode->op << " "<<curFstNode->lastPhone <<" "<<curFstNode->lastSIL <<endl;
	}
}

vector<Triphone> FstLinkedList::getNextPhoneFst(
		int fstStateNo) {
	vector<Triphone> nextPhoneFst;
	FstNode *curFstNode = &finalFstLinkedList[fstStateNo];

	while (curFstNode != NULL) {
		if (curFstNode->end != -1) // condition to check end of FST, last state of FST has only start. in and end are missing
				{

			if ((curFstNode->in.find("#") != std::string::npos)
					|| (curFstNode->in.find("<eps>") != std::string::npos)) {
				vector<Triphone> newList = getNextPhoneFst(curFstNode->end);
				for (int i = 0; i < (int) newList.size(); i++) {
					Triphone * newListTriphone = &newList[i];
					if (curFstNode->in.find("#1") != std::string::npos) {
						newListTriphone->newWord = 1;
						newListTriphone->wordEnd = 1;
						newListTriphone->position = 0;
						if (checkNewWord(curFstNode->op) == 1) {
							newListTriphone->word = curFstNode->op;
						}
					} else if (curFstNode->in.find("#0") != std::string::npos) {
						newListTriphone->position = 2;
						newListTriphone->wordEnd = 1;
						newListTriphone->newWord = 1;
						if (checkNewWord(curFstNode->op) == 1) { // this word is actually of the previous state.
																 // raise flag and let viterbiAligner handle it.
							newListTriphone->prevWordFlag = 1;
							newListTriphone->prevWord = curFstNode->op;
						}
					} else {
						newListTriphone->position = 3;
					}
				}
				nextPhoneFst.insert(nextPhoneFst.end(), newList.begin(), newList.end());
			} else {
				Triphone currentPhone;
				currentPhone.triphone = curFstNode->in;
				currentPhone.triphoneWeight = curFstNode->wt;
				currentPhone.lastSIL = curFstNode->lastSIL;
				currentPhone.lastPhone = curFstNode->lastPhone;
				currentPhone.linkedListStateNo = curFstNode->start;
				currentPhone.nextlinkedListStateNo = curFstNode->end;
				if (checkNewWord(curFstNode->op) == 1) {
					currentPhone.word = curFstNode->op;
				}
				nextPhoneFst.push_back(currentPhone);
			}
		} else {
			// FST ends here, last state of FST
			Triphone currentPhone;
			currentPhone.triphone = curFstNode->in;
			currentPhone.triphoneWeight = curFstNode->wt;
			currentPhone.lastSIL = curFstNode->lastSIL;
			currentPhone.lastPhone = curFstNode->lastPhone;
			currentPhone.linkedListStateNo = curFstNode->start;
			currentPhone.nextlinkedListStateNo = curFstNode->end;
			nextPhoneFst.push_back(currentPhone);
		}
		curFstNode = curFstNode->nextFstNode;
	}
	return nextPhoneFst;
}

vector<Triphone> FstLinkedList::getNextTriphoneFst(
		int fstStateNo,
		string prevTriphn) {
	vector<Triphone> nextPhoneFst;
	vector<Triphone> centerPhoneList = getNextPhoneFst(fstStateNo);
	int TriphonePosition = 1;				//0-beginning, 1-intermediate, 2-end
	if (prevTriphn.find("+") != std::string::npos) {
		prevTriphn = "SIL";
	}

	for (int i = 0; i < (int) centerPhoneList.size(); i++) {
		Triphone centerPhone = centerPhoneList[i];
		TriphonePosition = 1;
		if (centerPhone.position == 0) {
			TriphonePosition = 0;
		}
		if (!((centerPhone.triphone.find("+") != std::string::npos)
				|| (centerPhone.triphone.find("SIL") != std::string::npos))) {
			vector<Triphone> rightContextList = getNextPhoneFst(centerPhone.nextlinkedListStateNo);
			for (int j = 0; j < (int) rightContextList.size(); j++) {
				Triphone rightContext = rightContextList[i];
				Triphone currentPhone;
				if (rightContext.position == 2) {
					TriphonePosition = 2;
				}
				if (rightContext.nextlinkedListStateNo != -1) {
					if (!(rightContext.triphone.find("+") != std::string::npos)) {
						currentPhone.triphone = (centerPhone.triphone + "\t" + Dictionary::string2words(prevTriphn)[0]
								+ "\t" + rightContext.triphone);
						currentPhone.linkedListStateNo = centerPhone.linkedListStateNo;
						currentPhone.nextlinkedListStateNo = centerPhone.nextlinkedListStateNo;
					} else {
						currentPhone.triphone = (centerPhone.triphone + "\t" + Dictionary::string2words(prevTriphn)[0]
								+ "\t" + "SIL");
						currentPhone.linkedListStateNo = centerPhone.linkedListStateNo;
						currentPhone.nextlinkedListStateNo = centerPhone.nextlinkedListStateNo;
					}
				} else {
					currentPhone.triphone = (centerPhone.triphone + "\t" + Dictionary::string2words(prevTriphn)[0]
							+ "\t" + "SIL");
					currentPhone.linkedListStateNo = centerPhone.linkedListStateNo;
					currentPhone.nextlinkedListStateNo = centerPhone.nextlinkedListStateNo;
				}
				if (TriphonePosition == 0) {
					currentPhone.triphone = currentPhone.triphone + "\tb";
				} else {
					if (TriphonePosition == 1) {
						currentPhone.triphone = currentPhone.triphone + "\ti";
					} else {
						currentPhone.triphone = currentPhone.triphone + "\te";
					}
				}
				nextPhoneFst.push_back(currentPhone);
			}
		} else {
			centerPhone.triphone = centerPhone.triphone + "\t-\t-\t-";
			nextPhoneFst.push_back(centerPhone);
		}
	}
	return nextPhoneFst;
}

bool FstLinkedList::lastPhoneReached(
		int linkedListStateNo) {
	FstNode curFstNode = finalFstLinkedList[linkedListStateNo];
	if (curFstNode.end != -1) {
		return true;
	} else {
		return false;
	}
}

int FstLinkedList::getNextPhoneFstIndex(
		string phone,
		int curLinkedListIndex) {
	FstNode *curFstNode = &finalFstLinkedList[curLinkedListIndex];
	while (curFstNode != NULL) {
		if (phone.compare(curFstNode->in) == 0) {
			return curFstNode->end;
		} else {
			curFstNode = curFstNode->nextFstNode;
		}
	}
	return -1;
}

int FstLinkedList::checkNewWord(
		string word) {
	if ((word.compare("<eps>") == 0) || (word.compare("<sil>") == 0) || (word.compare("<s>") == 0)
			|| (word.compare("</s>") == 0)) {
		return 0;
	} else {
		return 1;
	}
}
} /* namespace std */
