#ifndef __HMMER3_WRAPPER__
#define __HMMER3_WRAPPER__

#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hmmer.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_gumbel.h"
#include "esl_exponential.h"

#define INIT_NUM_MODELS 1
#define INC_NUM_MODELS 10

#define MODELS_NOT_PRESSED -100

typedef struct {
  char* model_name;
  double bits;
  char* aligned_seq;

  int matches;
  int inserts;
  int deletes;
  int n;
  int c;

  int seq_start, seq_stop;
  int hmm_start, hmm_stop;

  size_t seq_alloc;
  size_t name_alloc;
} WRAPPER_RESULT;

WRAPPER_RESULT** wrapper_results;
int num_results;

//Setup hmmer3, loads ALL models in to memory
//Must be called before using hmmscan_seq
int init_hmmer_wrapper(const char* hmmfile);

//Scan the models to see which this seq may belong to, top hits
//will contain the results of the scan
void run_hmmer_pipeline(const char* seq);

//Free all hmmer components, be kind, rewind
void destroy_hmmer_wrapper();

#endif
