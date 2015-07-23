/*
 * FstLinkedList.h
 *
 *  Created on: Jul 24, 2014
 *      Author: sharath
 */

#ifndef FSTLINKEDLIST_H_
#define FSTLINKEDLIST_H_
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "commonUtils.h"
#include "AcousticModel.h"
#include "Dictionary.h"
#include "LatticeLinkedList.h"
namespace std {

class FstNode {
public:
	int start;
	int end;
	string in;
	string op;
	int wt;
	FstNode *nextFstNode;
	int lastSIL;
	int lastPhone;

	FstNode() {
		nextFstNode = NULL;
		lastSIL = 0;
		lastPhone = 0;
		wt = 0;
	}
	;

	int getStart() {
		return start;
	}

	void setStart(
			int newStart) {
		start = newStart;
	}

	int getEnd() {
		return end;
	}

	void setEnd(
			int newEnd) {
		end = newEnd;
	}

	string getIn() {
		return in;
	}

	void setIn(
			string newIn) {
		in = newIn;
	}

	string getOp() {
		return op;
	}

	void setOp(
			string newOp) {
		op = newOp;
	}
};

class FstLinkedList {
private:
	vector<FstNode> listAllFstNodes;
	vector<FstNode> finalFstLinkedList;
public:
	FstLinkedList() {
	}
	;
	void initFstLinkedList(
			string fstPath,
			string fileName);
	vector<Triphone> getNextPhoneFst(
			int fstStateNo);
	vector<Triphone> getNextTriphoneFst(
			int fstStateNo,
			string prevTriphn);
	bool lastPhoneReached(
			int linkedListStateNo);
	int getNextPhoneFstIndex(
			string phone,
			int curLinkedListIndex);
	int checkNewWord(
			string word);

	~FstLinkedList() {
	}
	;
};

} /* namespace std */
#endif /* FSTLINKEDLIST_H_ */
