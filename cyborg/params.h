/*
 * params.h
 *
 *  Created on: Feb 3, 2015
 *      Author: sharath
 */

#ifndef PARAMS_H_
#define PARAMS_H_

//#define SYS_FS 16000
#define SYS_FS 8000

#if (SYS_FS==16000)
#define UPPER_F			(6855.4976)
#define LOWER_F			(133.33334)
#define NUM_FILT		(40)
#define FS				(SYS_FS)
#define NFFT			(512)
#else
#define UPPER_F			(3500)
#define LOWER_F			(133)
#define NUM_FILT		(31)
#define FS				(SYS_FS)
#define NFFT			(256)
#endif

#define VAR_FLOOR		(0.0001)
#define MIXWT_FLOOR		(0.0000001)
#define TMAT_FLOOR		(0.0001)
#define BEAM_WIDTH		(1e-38) //Beam selecting active HMMs (relative to  best) in each frame [0(widest)..1(narrowest)]
//For widest beam, we cant have 0, because log(0) =Infinity.
// But the lowest float value supported 1e-38, so we use it.

#define INCOMPLETE_MODE	(2) //0 - force SIL at the end, 1 - best likelihood whatever alignment, 2- force align to transcript
#define N_BEST			(1) //Has to be greater than 1. 1- prints the best alignment, N- prints the N best alignments
#define LOGBASE			(1.0003)//loge=2.71828182846;
#define WINLEN_SEC		(0.025625)
#define DCT_TYPE		(0) //legacy=0, dct=1,htk=2

#endif /* PARAMS_H_ */
