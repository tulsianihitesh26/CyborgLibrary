/*
 * LatticeLinkedList.cpp
 *
 *  Created on: Jul 24, 2014
 *      Author: sharath
 */

#include "LatticeLinkedList.h"

namespace std {

LatticeLinkedList::LatticeLinkedList() {
	active = true;
	numChild = 0;
	position = 2;
	newWord = 0;
	wordEnd = 0;
	word = "-";
}

IntArray LatticeLinkedList::getSiblingsBranchIds() {
	return siblingsBranchIds;
}

void LatticeLinkedList::setSiblingsBranchIds(
		IntArray newSiblingsBranchIds) {
	siblingsBranchIds = newSiblingsBranchIds;
}

int LatticeLinkedList::getBranchId() {
	return branchId;
}

void LatticeLinkedList::setBranchId(
		int newBranchId) {
	branchId = newBranchId;
}

bool LatticeLinkedList::isActive() {
	return active;
}
void LatticeLinkedList::setActive(
		bool newActive) {
	active = newActive;
}

LatticeLinkedList* LatticeLinkedList::getParent() {
	return parent;
}

void LatticeLinkedList::setParent(
		LatticeLinkedList * newParent) {
	parent = newParent;
}

string LatticeLinkedList::getTriPhone() {
	return triPhone;
}

void LatticeLinkedList::setTriPhone(
		string newTriPhone) {
	triPhone = newTriPhone;
}

vector<LatticeLinkedList> LatticeLinkedList::getChildren() {
	return children;
}
void LatticeLinkedList::setChildren(
		vector<LatticeLinkedList> newChildren) {
	children = newChildren;
}

Number LatticeLinkedList::getCostOfNode() {
	return costOfNode;
}

void LatticeLinkedList::setCostOfNode(
		Number newCostOfNode) {
	costOfNode = newCostOfNode;
}

Number LatticeLinkedList::getCostOfPath() {
	return costOfPathIncludingLw;
}

void LatticeLinkedList::setCostOfPath(
		Number newCostOfPath) {
	costOfPathIncludingLw = newCostOfPath;
}

void LatticeLinkedList::setChild(
		LatticeLinkedList child) {
	children.push_back(child);
}

void LatticeLinkedList::copyFormAnotherNode(
		LatticeLinkedList & curLatticeLinkedList) {
	curPhoneStateNo = curLatticeLinkedList.curPhoneStateNo;
	curStateNo = curLatticeLinkedList.curStateNo;
	observationId = curLatticeLinkedList.observationId;
	stateIds = curLatticeLinkedList.stateIds;
	tmatId = curLatticeLinkedList.tmatId;
	linkedListStateNo = curLatticeLinkedList.linkedListStateNo;
	nextlinkedListStateNo = curLatticeLinkedList.nextlinkedListStateNo;
	triPhone = curLatticeLinkedList.triPhone;
//	triPhone = new string();
//	triPhone = string.format("%s",curLatticeLinkedList.triPhone);
}

} /* namespace std */
