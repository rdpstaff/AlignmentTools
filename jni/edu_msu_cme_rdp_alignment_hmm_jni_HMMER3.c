/*
 * Copyright (C) 2012 Michigan State University <rdpstaff at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "edu_msu_cme_rdp_alignment_hmm_jni_HMMER3.h"
#include "hmmer_wrapper.h"
#include "easel.h"

static jclass hitClass;
static jfieldID nameFieldId;
static jfieldID bitsFieldId;
static jfieldID alignedSeqFieldId;
static jfieldID hmmStartFieldId;
static jfieldID hmmEndFieldId;
static jfieldID seqStartFieldId;
static jfieldID seqEndFieldId;
static jmethodID hitConstructor;

void JNU_ThrowByName(JNIEnv *env, const char *name, const char *msg)
{
  jclass cls = (*env)->FindClass(env, name);
  /* if cls is NULL, an exception has already been thrown */
  if (cls != NULL) {
    (*env)->ThrowNew(env, cls, msg);
  }
  /* free the local ref */
  (*env)->DeleteLocalRef(env, cls);
}

JNIEXPORT jint JNICALL Java_edu_msu_cme_rdp_alignment_hmm_jni_HMMER3_initHmmer(JNIEnv *env, jobject jobj, jstring f_name) {
  jboolean iscopy;
  const char *hmmfile = (*env)->GetStringUTFChars(env, f_name, &iscopy);
  int result = init_hmmer_wrapper(hmmfile);
  (*env)->ReleaseStringUTFChars(env, f_name, hmmfile);

  if(result <= 0) {
    if(result == MODELS_NOT_PRESSED) {
      JNU_ThrowByName(env, "java.io.IOException", "HMM file not pressed!");
    } else if(result == -eslENOTFOUND) {
      JNU_ThrowByName(env, "java.io.IOException", "HMM file not found!");
    } else if(result == -eslEFORMAT) {
      JNU_ThrowByName(env, "java.io.IOException", "HMM file is not a recognized format!");
    } else if(result == 0) {
      JNU_ThrowByName(env, "java.io.IOException", "Didn't read any models from the file");
    } else {
      JNU_ThrowByName(env, "java.io.IOException", "Unknown error encountered when reading hmm file");
    }
    return 0;
  }

  jclass cls = (*env)->FindClass(env, "edu/msu/cme/rdp/alignment/hmm/jni/HMMER3Hit");
  if(cls == NULL) {
    JNU_ThrowByName(env, "java/lang/ClassNotFoundException", "edu/msu/cme/rdp/alignment/hmm/jni/HMMER3Hit");
    return 0;
  }

  hitClass = (*env)->NewGlobalRef(env, cls);
  (*env)->DeleteLocalRef(env, cls);

  hitConstructor = (*env)->GetMethodID(env, hitClass, "<init>", "()V");
  if(hitConstructor == NULL) {
    JNU_ThrowByName(env, "java/lang/NoSuchMethodException", "<init> ()V on class HMMER3$HMMER3Hit");
    return 0;
  }

  nameFieldId = (*env)->GetFieldID(env, hitClass, "modelName", "Ljava/lang/String;");
  if(nameFieldId == NULL) {
    JNU_ThrowByName(env, "java/lang/NoSuchFieldException", "Field 'modelName' (type String) does not exist in class HMMER3Hit");
    return 0;
  }

  bitsFieldId = (*env)->GetFieldID(env, hitClass, "bits", "D");
  if(bitsFieldId == NULL) {
    JNU_ThrowByName(env, "java/lang/NoSuchFieldException", "Field 'bits' (type double) does not exist in class HMMER3Hit");
    return 0;
  }

  alignedSeqFieldId = (*env)->GetFieldID(env, hitClass, "alignedSeq", "Ljava/lang/String;");
  if(alignedSeqFieldId == NULL) {
    JNU_ThrowByName(env, "java/lang/NoSuchFieldException", "Field 'alignedSeq' (type String) does not exist in class HMMER3Hit");
    return 0;
  }

  hmmStartFieldId = (*env)->GetFieldID(env, hitClass, "hmmStart", "I");
  if(hmmStartFieldId == NULL) {
    JNU_ThrowByName(env, "java/lang/NoSuchFieldException", "Field 'hmmStart' (type int) does not exist in class HMMER3Hit");
    return 0;
  }

  hmmEndFieldId = (*env)->GetFieldID(env, hitClass, "hmmEnd", "I");
  if(hmmEndFieldId == NULL) {
    JNU_ThrowByName(env, "java/lang/NoSuchFieldException", "Field 'hmmEnd' (type int) does not exist in class HMMER3Hit");
    return 0;
  }

  seqStartFieldId = (*env)->GetFieldID(env, hitClass, "seqStart", "I");
  if(seqStartFieldId == NULL) {
    JNU_ThrowByName(env, "java/lang/NoSuchFieldException", "Field 'seqStart' (type int) does not exist in class HMMER3Hit");
    return 0;
  }

  seqEndFieldId = (*env)->GetFieldID(env, hitClass, "seqEnd", "I");
  if(seqEndFieldId == NULL) {
    JNU_ThrowByName(env, "java/lang/NoSuchFieldException", "Field 'seqEnd' (type int) does not exist in class HMMER3Hit");
    return 0;
  }

  return result;
}

/*
 * Class:     edu_msu_cme_rdp_alignment_hmm_jni_HMMER3
 * Method:    destroyHmmer
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_edu_msu_cme_rdp_alignment_hmm_jni_HMMER3_destroyHmmer(JNIEnv *env, jobject jobj) {
  destroy_hmmer_wrapper();

  if(hitClass != NULL) {
    (*env)->DeleteGlobalRef(env, hitClass);
  }
}

/*
 * Class:     edu_msu_cme_rdp_alignment_hmm_jni_HMMER3
 * Method:    hmmer3
 * Signature: (Ljava/lang/String;)Ledu/msu/cme/rdp/alignment/hmm/jni/HMMER3/HMMER3Hit;
 */
JNIEXPORT jobjectArray JNICALL Java_edu_msu_cme_rdp_alignment_hmm_jni_HMMER3_hmmer3(JNIEnv *env, jobject jobj, jstring jseq) {
  jboolean iscopy;
  const char* seq = (*env)->GetStringUTFChars(env, jseq, &iscopy);
  run_hmmer_pipeline(seq);

  int index;
  jobjectArray ret = (jobjectArray)(*env)->NewObjectArray(env, num_results, hitClass, NULL);

  for(index = 0;index < num_results;index++) {
    WRAPPER_RESULT* r = wrapper_results[index];

    jobject hit = (*env)->NewObject(env, hitClass, hitConstructor);
    jstring name = (*env)->NewStringUTF(env, r->model_name);
    jstring alignedSeq  = (*env)->NewStringUTF(env, r->aligned_seq);

    (*env)->SetObjectField(env, hit, nameFieldId, name);
    (*env)->SetObjectField(env, hit, alignedSeqFieldId, alignedSeq);
    (*env)->SetDoubleField(env, hit, bitsFieldId, r->bits);
    (*env)->SetIntField(env, hit, hmmStartFieldId, r->hmm_start);
    (*env)->SetIntField(env, hit, hmmEndFieldId, r->hmm_stop);
    (*env)->SetIntField(env, hit, seqStartFieldId, r->seq_start);
    (*env)->SetIntField(env, hit, seqEndFieldId, r->seq_stop);

    (*env)->SetObjectArrayElement(env, ret, index, hit);
  }

  (*env)->ReleaseStringUTFChars(env, jseq, seq);

  return ret;
}

