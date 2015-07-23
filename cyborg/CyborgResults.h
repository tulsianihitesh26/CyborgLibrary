/*
 * CyborgResults.h
 *
 *  Created on: Feb 3, 2015
 *      Author: sharath
 */

#ifndef CYBORGRESULTS_H_
#define CYBORGRESULTS_H_

struct CyborgResultsContainer {
	StringArray phone;
	StringArray word;
	IntArray frameStart;
	IntArray frameEnd;
	IntArray acousticScore;
	IntArray transitionScore;
	IntArray score;
	IntArray posterior;
	IntArray state;
	IntArray phoneState;
	IntArray position;
	IntArray wordEnd;
	int size;
};

typedef CyborgResultsContainer & CyborgResultsContainerRef;

#endif /* CYBORGRESULTS_H_ */
