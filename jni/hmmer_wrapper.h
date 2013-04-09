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
