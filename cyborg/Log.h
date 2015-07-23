/*
 * Log.h
 *
 *  Created on: 14-Jul-2014
 *      Author: Sujeet
 */

#ifndef LOG_H_
#define LOG_H_

#define JNI_DEBUG

#ifdef ANDROID_BUILD
#include <jni.h>
#include <GLES2/gl2.h>
#ifdef JNI_DEBUG
#include <android/log.h>
#define  LOG_TAG    "libradioke"
#define  LOGV(...)  __android_log_print(ANDROID_LOG_VERBOSE,LOG_TAG,__VA_ARGS__)
#define  LOGD(...)  __android_log_print(ANDROID_LOG_DEBUG,LOG_TAG,__VA_ARGS__)
#define  LOGI(...)  __android_log_print(ANDROID_LOG_INFO,LOG_TAG,__VA_ARGS__)
#define  LOGE(...)  __android_log_print(ANDROID_LOG_ERROR,LOG_TAG,__VA_ARGS__)
#else
#define  LOGV(...)
#define  LOGD(...)
#define  LOGI(...)
#define  LOGE(...)
#endif

#else
#include <stdio.h>
//#include "jni.h"
#define  LOGV(...)  fprintf(stdout,__VA_ARGS__)
#define  LOGD(...)  fprintf(stdout,__VA_ARGS__)
#define  LOGI(...)  fprintf(stdout,__VA_ARGS__)
#define  LOGE(...)  fprintf(stderr,__VA_ARGS__)
#endif /*ANDROID_BUILD*/

#endif /* LOG_H_ */
