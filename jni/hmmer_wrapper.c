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

#include "hmmer_wrapper.h"

static P7_OPROFILE** models = NULL;
static P7_PROFILE** gmodels = NULL;
static P7_BG* bg = NULL;
static P7_HMMFILE* hmm_fp = NULL;
static ESL_ALPHABET* abc = NULL;
static P7_OMX* oxf = NULL;
static P7_OMX* oxb = NULL;
static P7_GMX* gxf = NULL;
static P7_GMX* gxb = NULL;
static P7_TRACE* tr = NULL;

const double f1 = 0.02;
const double f2 = 1e-3;
const double f3 = 1e-5;

static int num_models = 0;
static int models_capacity = 0;
static int hmmer_error = 0;

void hmmer_error_handler(int errcode, int use_errno, char *sourcefile, int sourceline, char *format, va_list argp) {
  vfprintf(stderr, format, argp);
  fprintf(stderr, "\n");
  hmmer_error = 1;
  return;
}

WRAPPER_RESULT* new_result() {
  WRAPPER_RESULT* ret = (WRAPPER_RESULT*)malloc(sizeof(WRAPPER_RESULT));
  ret->model_name = (char*)malloc(sizeof(char) * 40);
  ret->name_alloc = 40;

  ret->aligned_seq = (char*)malloc(sizeof(char) * 200);
  ret->seq_alloc = 200;

  ret->matches = 0;
  ret->deletes = 0;
  ret->inserts = 0;
  ret->n = 0;
  ret->c = 0;

  return ret;
}

void reuse_result(WRAPPER_RESULT* result, int N, const char* model_name) {
  if(result->seq_alloc < N + 1) {
    result->aligned_seq = (char*)realloc(result->aligned_seq, sizeof(char) * (N + 1));
    result->seq_alloc = N + 1;
  }

  if(model_name != NULL) {
    int name_len = strlen(model_name);
    if(name_len > result->name_alloc) {
      result->model_name = (char*)realloc(result->model_name, sizeof(char) * name_len);
      result->name_alloc = name_len;
    }
    
    strcpy(result->model_name, model_name);
  }
  result->matches = 0;
  result->deletes = 0;
  result->inserts = 0;
  result->n = 0;
  result->c = 0;
}

void destroy_result(WRAPPER_RESULT* result) {
  if(result && result->aligned_seq) {
    free(result->aligned_seq);
  }

  if(result && result->model_name) {
    free(result->model_name);
  }

  if(result) {
    free(result);
  }
}

int init_hmmer_wrapper(const char* hmmfile) {
  p7_FLogsumInit();
  impl_Init();
  esl_exception_SetHandler(hmmer_error_handler);

  int status, index;
  
  P7_HMM* model = NULL;
  P7_PROFILE* gm = NULL;
  P7_OPROFILE* om = NULL;

  status = p7_hmmfile_Open(hmmfile, p7_HMMDBENV, &hmm_fp);

  if(status != eslOK) {
    if(hmm_fp) {
      p7_hmmfile_Close(hmm_fp);
    }

    if(!hmm_fp->is_pressed) {
      return MODELS_NOT_PRESSED;
    }

    return -status;
  }

  model = NULL;

  while((status = p7_hmmfile_Read(hmm_fp, &abc, &model)) == eslOK) {
    if(bg == NULL) {
      bg = p7_bg_Create(abc);
      p7_bg_SetLength(bg, 400);
    }

    gm = p7_profile_Create(model->M, abc);
    om = p7_oprofile_Create(model->M, abc);

    p7_ProfileConfig(model, bg, gm, 400, p7_UNILOCAL);
    p7_oprofile_Convert(gm, om);
  
    /*while((status = p7_oprofile_ReadMSV(hmm_fp, &abc, &om)) == eslOK) {
      p7_oprofile_ReadRest(hmm_fp, om);*/

    if(num_models >= models_capacity) {
      models_capacity += INC_NUM_MODELS;
      models = (P7_OPROFILE**)realloc(models, sizeof(P7_OPROFILE*) * models_capacity);
      gmodels = (P7_PROFILE**)realloc(gmodels, sizeof(P7_PROFILE*) * models_capacity);
    }

    models[num_models] = om;
    gmodels[num_models] = gm;

    p7_hmm_Destroy(model);

    model = NULL;
    om = NULL;
    gm = NULL;
    num_models++;
  }

  if(models == 0) {
    return 0;
  }

  tr = p7_trace_CreateWithPP();

  oxf = p7_omx_Create(400, 400, 400);
  oxb = p7_omx_Create(400, 400, 400);
  wrapper_results = (WRAPPER_RESULT**)malloc(sizeof(WRAPPER_RESULT*) * num_models);
  for(index = 0;index < num_models;index++) {
    wrapper_results[index] = new_result();
  }

  num_results = 0;

  return num_models;
}

void trace_into(P7_TRACE* tr, WRAPPER_RESULT* result, ESL_SQ* sq, ESL_ALPHABET* abc, int M) {
  int k = 0, i = 0;
  int z1, z2, z;

  for (z1 = 0;z1 < tr->N; z1++) if (tr->st[z1] == p7T_M) break;
  result->seq_start = tr->i[z1];
  result->hmm_start = tr->k[z1];

  for (z=z1; z < tr->N; z++) if (tr->st[z] == p7T_M) z2 = z;
  result->seq_stop = tr->i[z2];
  result->hmm_stop = tr->k[z2];

  /*for(i = 1;i <= result->seq_start;i++) {
    result->aligned_seq[k++] = (char)tolower((int)abc->sym[sq->dsq[i]]);
  }
  for(i = 1;i <= result->hmm_start;i++) {
    result->aligned_seq[k++] = '-';
    }*/

  int numB = 0, numE = 0;
  for(z = 2;z < tr->N;z++) {
    int xi = sq->dsq[tr->i[z]];

    char sym = abc->sym[xi];
    switch(tr->st[z]) {
    case p7T_M:
      if(xi == eslDSQ_SENTINEL) {
	continue;
      }

      result->matches++;
      result->aligned_seq[k++] = sym;
      break;
    case p7T_D:
      result->deletes++;
      result->aligned_seq[k++] = '-';
      break;
    case p7T_I:
      if(xi == eslDSQ_SENTINEL) {
	continue;
      }

      result->inserts++;
      result->aligned_seq[k++] = (char)tolower((int)sym);
      break;
    case p7T_N:
      if(xi == eslDSQ_SENTINEL) {
	continue;
      }
      result->n++;
      result->aligned_seq[k++] = (char)tolower((int)sym);
      break;
    case p7T_C:
      if(xi == eslDSQ_SENTINEL) {
	continue;
      }
      result->c++;
      result->aligned_seq[k++] = (char)tolower((int)sym);
      break;
    case p7T_B:
      if(!numB) {
	for(i = 1;i < result->hmm_start;i++) {
	  result->aligned_seq[k++] = '-';
	}
      }
      numB++;
      break;
    case p7T_E:
      if(!numE) {
	for(i = result->hmm_stop;i < M;i++) {
	  result->aligned_seq[k++] = '-';
	}
      }
      numE++;
      break;
    }
  }

  /*for(i = result->hmm_stop;i < M;i++) {
    result->aligned_seq[k++] = '-';
  }
  for(i = result->seq_stop + 1;i <= sq->n;i++) {
    result->aligned_seq[k++] = (char)tolower((int)abc->sym[sq->dsq[i]]);
    }*/
  
  result->aligned_seq[k] = '\0';
}

void run_hmmer_pipeline(const char* seq) {
  int index, i, status;
  ESL_SQ* sq = esl_sq_CreateFrom(NULL, seq, NULL, NULL, NULL);
  P7_OPROFILE *om = NULL;
  P7_PROFILE *gm = NULL;
  float usc, vfsc, fwdsc;   /* filter scores                           */
  float filtersc;           /* HMM null filter score                   */
  float nullsc;             /* null model score                        */
  float seqbias;
  float seq_score;          /* the corrected per-seq bit score */
  double P;
  WRAPPER_RESULT* result;

  num_results = 0;
  if(sq->n == 0) {
    esl_sq_Destroy(sq);
    return;
  }

  esl_sq_Digitize(abc, sq);  

  int n = 0;
  float oasc;

  for(index = 0;index < num_models;index++) {
    om = models[index];

    p7_omx_Reuse(oxf);
    p7_omx_Reuse(oxb);

    p7_omx_GrowTo(oxf, om->M, sq->n, sq->n);
    p7_omx_GrowTo(oxb, om->M, sq->n, sq->n);

    p7_oprofile_ReconfigLength(om, sq->n);

    p7_bg_SetFilter(bg, om->M, om->compo);
    p7_bg_SetLength(bg, sq->n);

    //Calibrate null model
    p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

    //MSV Filter
    p7_MSVFilter(sq->dsq, sq->n, om, oxf, &usc);
    seq_score = (usc - nullsc) / eslCONST_LOG2;
    P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
    if (P > f1) continue;

    //Bias filter (model compo)
    p7_bg_FilterScore(bg, sq->dsq, sq->n, &filtersc);
    seq_score = (usc - filtersc) / eslCONST_LOG2;
    P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
    if (P > f1) continue;

    //Viterbi filter (Only do if P value from Bias is high)
    if(P > f2) {
      p7_ViterbiFilter(sq->dsq, sq->n, om, oxf, &vfsc);
      seq_score = (vfsc - filtersc) / eslCONST_LOG2;
      P = esl_gumbel_surv(seq_score,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
      if (P > f2) continue;
    }

    //Get the real probability (forward)
    p7_Forward(sq->dsq, sq->n, om, oxf, &fwdsc);
    seq_score = (fwdsc - filtersc) / eslCONST_LOG2;
    P = esl_exp_surv(seq_score,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
    if(hmmer_error) {
      fprintf(stderr, "HMM: %s, seq: %s", om->name, seq);
      hmmer_error = 0;
      continue;
    }
    if (P > f3) continue;

    //Real hit, go in to posterior decoding and alignment
    p7_omx_Reuse(oxb);
    p7_trace_Reuse(tr);

    p7_Backward(sq->dsq, sq->n, om, oxf, oxb, NULL);

    status = p7_Decoding(om, oxf, oxb, oxb);

    if(status == eslOK) {
      //And then trace the result
      p7_OptimalAccuracy(om, oxb, oxf, &oasc);
      p7_OATrace(om, oxb, oxf, tr);
    } else if(status == eslERANGE) {
      fprintf(stderr, "Decoding overflow on model %s\n", om->name);
      gm = gmodels[index];
      if(gxf == NULL) {
	gxf = p7_gmx_Create(gm->M, sq->n);
	gxb = p7_gmx_Create(gm->M, sq->n);
      } else {
	p7_gmx_GrowTo(gxf, gm->M, sq->n);
	p7_gmx_GrowTo(gxb, gm->M, sq->n);
      }

      p7_ReconfigLength(gm, sq->n);

      p7_GForward (sq->dsq, sq->n, gm, gxf, &fwdsc);
      p7_GBackward(sq->dsq, sq->n, gm, gxb, NULL);

      p7_GDecoding(gm, gxf, gxb, gxb);
      p7_GOptimalAccuracy(gm, gxb, gxf, &oasc);
      p7_GOATrace        (gm, gxb, gxf, tr);

      p7_gmx_Reuse(gxf);
      p7_gmx_Reuse(gxb);
    }

    if(hmmer_error) {
      fprintf(stderr, "HMM: %s, seq: %s", om->name, seq);
      hmmer_error = 0;
      continue;
    }

    result = wrapper_results[num_results];
    reuse_result(result, tr->N + om->M, om->name); //We're way overallocating here, but it's hard to know at this point how much space we'll need for the alignment (plus leading and trailing gaps)
    trace_into(tr, result, sq, abc, om->M);
    result->bits = seq_score;
    num_results++;
  }

  esl_sq_Destroy(sq);
}

void destroy_hmmer_wrapper() {
  int index;
  if(models != NULL) {
    for(index = 0;index < num_models;index++) {
      p7_oprofile_Destroy(models[index]);
      p7_profile_Destroy(gmodels[index]);
    }
    free(models);
    free(gmodels);
  }

  if(wrapper_results != NULL) {
    for(index = 0;index < num_models;index++) {
      destroy_result(wrapper_results[index]);
    }
    free(wrapper_results);
  }
  
  if(bg != NULL) {
    p7_bg_Destroy(bg);
  }

  if(hmm_fp != NULL) {
    p7_hmmfile_Close(hmm_fp);
  }

  if(oxf) {
    p7_omx_Destroy(oxf);
  }
  if(oxb) {
    p7_omx_Destroy(oxb);
  }

  if(gxf) {
    p7_gmx_Destroy(gxf);
  }
  if(gxb) {
    p7_gmx_Destroy(gxb);
  }

  if(abc) {
    esl_alphabet_Destroy(abc);
  }
  if(tr) {
    p7_trace_Destroy(tr);
  }
}

int main(int argc, char** argv) {
  if(argc != 3) {
    fprintf(stderr, "USAGE: wrapper <hmm> <seqs>\n");
    exit(1);
  }

  const char* hmm_file = argv[1];
  char* line = malloc(sizeof(char) * 4096);
  FILE* seqfile = fopen(argv[2], "r");
  int index;

  if(seqfile == NULL) {
    perror("Failed to open seq file");
    free(line);
    exit(1);
  } else {
    int n_models = init_hmmer_wrapper(hmm_file);
    if(n_models < 1) {
      fprintf(stderr, "Couldn't open file %s", hmm_file);
      exit(1);
    }

    printf("Created %d models\n", n_models);
    int lineno = 0;
    while(fgets(line, 4096, seqfile) != NULL ) {
      line[strlen(line) - 1] = '\0';
      lineno++;
      printf("%d\t%s\n", lineno, line);

      run_hmmer_pipeline(line);
      for(index = 0;index < num_results;index++) {
	WRAPPER_RESULT* r = wrapper_results[index];
	printf("Model %s, score=%f, matches=%d, inserts=%d, deletes=%d, n=%d, s=%d\n", r->model_name, r->bits, r->matches, r->inserts, r->deletes, r->n, r->c);
	printf("hmm: %d - %d, seq: %d - %d\n", r->hmm_start, r->hmm_stop, r->seq_start, r->seq_stop);
	printf("Alignment: %s\n", r->aligned_seq);
      }

      printf("Pipeline finished\n");
    }

    destroy_hmmer_wrapper();
    fclose(seqfile);
    free(line);
  }
}
