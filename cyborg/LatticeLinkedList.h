/*
 * LatticeLinkedList.h
 *
 *  Created on: Jul 24, 2014
 *      Author: sharath
 */

#ifndef LATTICELINKEDLIST_H_
#define LATTICELINKEDLIST_H_
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "commonUtils.h"

namespace std {

class Triphone {
public:
	string triphone;
	int linkedListStateNo;
	int nextlinkedListStateNo;
	int position;				//0-beginning, 1-intermediate, 2-end
	int newWord;
	int wordEnd;
	int triphoneWeight;
	string word;
	string prevWord;
	int prevWordFlag;
	int lastSIL;
	int lastPhone;

	Triphone() {
		position = 1;
		newWord = 0;
		wordEnd = 0;
		lastSIL = 0;
		lastPhone = 0;
		word = "-";
		triphoneWeight = 0;
		prevWordFlag = 0;
	}
	;

	string getTriphone() {
		return triphone;
	}

	string getWord() {
		return word;
	}
	int getlinkedListStateNo() {
		return linkedListStateNo;
	}

	int getNextLinkedListStateNo() {
		return nextlinkedListStateNo;
	}
};

class LatticeLinkedList {
public:
	LatticeLinkedList();
	//Reference to parent TreeNode
	LatticeLinkedList *parent;

	/* stateInfo[0] : State ID index
	 * stateInfo[1] : Observation id
	 * stateInfo[2] : State tmat value
	 */
	//int stateInfo[];
	IntArray stateIds;
	int tmatId;
	int observationId;
	int curStateNo;
	int curPhoneStateNo;
	int linkedListStateNo;
	int nextlinkedListStateNo;
	int position;
	int newWord;
	int wordEnd;
	string word;
	int branchId;
	int lastSIL;
	int lastPhone;
	//Corresponding triphone
	string triPhone;

	//Stores languageScore of each frame
	int languageScore;

	//Stores acousticScore of each frame
	int acousticScore;
	int antiPhoneAcousticScore;

	//Stores transitionScore of each frame
	int transitionScore;

	//Stores likelihood of each frame
	int costOfNode;
	int antiPhoneCostOfNode;

	//Stores likelihood score of each path
	int costOfPath;
	int costOfPathIncludingLw;
	int antiPhoneCostOfPath;

	//Stores likelihood score of each path
	int numChild;

	bool active;

	//Stores list of children in left to right order.
	vector<LatticeLinkedList> children;

	//Stores branchIds in case
	IntArray siblingsBranchIds;
	IntArray getSiblingsBranchIds();
	void setSiblingsBranchIds(
			IntArray siblingsBranchIds);
	int getBranchId();
	void setBranchId(
			int newBranchId);
	bool isActive();
	void setActive(
			bool newActive);
	LatticeLinkedList *getParent();
	void setParent(
			LatticeLinkedList * newParent);
	string getTriPhone();
	void setTriPhone(
			string newTriPhone);
	vector<LatticeLinkedList> getChildren();
	void setChildren(
			vector<LatticeLinkedList> newChildren);
	Number getCostOfNode();
	void setCostOfNode(
			Number newCostOfNode);
	Number getCostOfPath();
	void setCostOfPath(
			Number newCostOfPath);
	void setChild(
			LatticeLinkedList child);
	void copyFormAnotherNode(
			LatticeLinkedList & curLatticeLinkedList);

};

} /* namespace std */
#endif /* LATTICELINKEDLIST_H_ */
