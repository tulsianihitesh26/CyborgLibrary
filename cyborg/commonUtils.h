#ifndef COMMONUTILS_H_
#define COMMONUTILS_H_

#include <iostream>
#include <vector>
#include "Log.h"
using namespace std;

#define M_PI 3.14159265359
#define LOGPROB_ZERO -939524096	/** Integer version of log of zero Approx -infinity!! */

typedef float Number;
typedef vector<Number> NumberArray;
typedef vector<Number> & NumberArrayRef;
typedef vector<Number> const & NumberArrayConstRef;
typedef vector<vector<Number> > Number2DArray;
typedef vector<vector<Number> > & Number2DArrayRef;
typedef vector<vector<vector<Number> > > Number3DArray;
typedef vector<vector<vector<Number> > > & Number3DArrayRef;
typedef vector<int> IntArray;
typedef IntArray & IntArrayRef;
typedef vector<vector<int> > Int2DArray;
typedef vector<vector<int> > & Int2DArrayRef;
typedef vector<vector<vector<int> > > Int3DArray;
typedef vector<vector<vector<int> > > & Int3DArrayRef;

typedef vector<string> StringArray;
typedef vector<string>& StringArrayRef;
typedef vector<Number> const & PitchesConstRef;

#endif
