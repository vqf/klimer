#ifndef KMERREAD_H_INCLUDED
#define KMERREAD_H_INCLUDED

#include "kmer.h"
#include "traceIdList.h"
#include "stdio.h"

typedef enum { absent, present, fragmented} testSeq;


char getLastBase(kmerHolder** khp, uint32_t pos){
  kmerHolder* kh = *khp;
  uint8_t nb = kh->nBases;
  uint8_t md = (uint8_t) (pos % nb);
  char result = kh->ms->bi->valBase[md];
  return result;
}

kcLL* followTrace(kmerHolder** khp, uint32_t pos, LISTTYPE tid, LISTTYPE posInTrace){
  if (pos == NOKMER) return NULL;
  kcLL* result = NULL;
  memstruct* ms = (*khp)->ms;
  tIdList* thisTrace = NULL;
  kmerConnector* kc = getKcWithTId(&ms, pos, tid, posInTrace, &thisTrace);
  if (!kc){
    D_(0, "Error at kc %lu, id %lu, pos %lu\n", (LUI) pos, (LUI) tid, (LUI) posInTrace);
    SEQ(&ms, ms->kmerArray[pos]);
    X_;
  }
  while (kc && !isLastInTrace(&thisTrace, posInTrace)){
    E_(1,
      char* seq2 = calloc(12, sizeof(char));
      pos2seq(&ms, kc->dest, seq2);
      printKmerConnector(kc, seq2);
      printTIdList(kc->idflags);
      free(seq2);
      //getc(stdin)
    );
    kcpush(&result, &kc, tid, posInTrace);
    kc = nextKc(&ms, &kc, &thisTrace, posInTrace);
    if (!kc){
      D_(0, "Trace %lu interrupted\n", (LUI) tid);
    }
    posInTrace++;
  }
  if (kc) kcpush(&result, &kc, tid, posInTrace);
  return result;
}

kcLL* nextTrace(kmerHolder** khp){
  kmerHolder* kh = *khp;
  memstruct* ms = kh->ms;
  kcLL* result = NULL;
  uint32_t i = 0;//ms->status->current;
  while (i < ms->nPos){
    kmerConnector* kc = ms->kmerArray[i];
    while (kc && kc->n > 0){
      traceVessel* l = traceFirst(&kc->idflags);
      while (l){
        tIdList* start = l->tidl;
        if (!IS(start, IN_USE)){
          SET(start, IN_USE);
          kcpush(&result, &kc, start->trace.n, 0);
          result->cPos = i;
          D_(1, "Found starting trace at %lu\n", (LUI) i);
          ms->status->current = i;
          kcLL* f = followTrace(khp, i, start->trace.n, 1);
          result->next = f;
          return result;
        }
        l = l->next;
      }
      destroyTraceVessel(&l);
      kc = kc->next;
    }
    i++;
  }
  return result;
}

char* getTraceSeq(kmerHolder** khp, kcLL** kcp){
  kmerHolder* kh = *khp;
  kcLL* kcll = *kcp;
  kcLL* tmp = kcll;
  uint32_t l = 0;
  SCALAR(tmp, l);
  uint32_t cl = (uint32_t) (kh->kmerSize + l + 1);
  char* result = (char*) calloc((size_t) cl, sizeof(char));
  pos2seq(&kh->ms, kcll->cPos, result);
  //DIE("%s\n", result);
  uint32_t j = (uint32_t) ((kh->kmerSize));
  while (kcll){
    result[j] = getLastBase(khp, kcll->kc->dest);
    j++;
    kcll = kcll->next;
  }
  return result;
}


seqCollection* recursivelyFollowTrace(kmerHolder** khp, kmerConnector** kcp, LISTTYPE tid, LISTTYPE posInTrace){
  seqCollection* result = newSeqCollection();
  memstruct* ms = (*khp)->ms;
  kmerConnector* kc = *kcp;
  kcpush(&result->trace, &kc, tid, posInTrace);
  tIdList* thisTrace = _getTrace(&kc->idflags, tid, posInTrace);
  kc = nextKc(&ms, &kc, &thisTrace, posInTrace);
  posInTrace++;
  if (!kc) return result;
  while (kc && thisTrace){
    kcpush(&result->trace, &kc, tid, posInTrace);
    kc = nextKc(&ms, &kc, &thisTrace, posInTrace);
    posInTrace++;
  }
  return result;
}

seqCollection* allTraces(kmerHolder** khp){
  seqCollection* result = newSeqCollection();
  seqCollection* pointer = result;
  memstruct* ms = (*khp)->ms;
  uint32_t i = 0;
  for (i = 0; i < ms->nPos; i++){
    D_(3, "Reading %lu\n", (LUI) i);
    kmerConnector* tc = ms->kmerArray[i];
    while (tc){
      tIdList* ti = tc->idflags;
      while (ti){
        if (IS(ti, FIRST_IN_TRACE)){
          D_(1, "Found starting trace at %lu\n", (LUI) i);
          seqCollection* newseqs = recursivelyFollowTrace(khp, &tc, ti->trace.n, 1);
          pointer->next = newseqs;
          pointer = newseqs;
        }
        ti = ti->next;
      }
      tc = tc->next;
    }
  }
  return result;
}

uint16_t getTotalNumberOfReads(kcLL** kclp){
  uint16_t result = 0;
  kcLL* tmp = *kclp;
  while (tmp){
    result += tmp->kc->n;
    tmp = tmp->next;
  }
  return result;
}

seqCollection* bestSeqCollection(seqCollection** scp){
  seqCollection* sc = *scp;
  seqCollection* result = sc;
  uint16_t top = getTotalNumberOfReads(&sc->trace);
  while (sc && sc->next){
    sc = sc->next;
    uint16_t n = getTotalNumberOfReads(&sc->trace);
    if (n > top){
      top = n;
      result = sc;
    }
  }
  return result;
}

void printSeq(kmerHolder** khp, seqCollection** scp){
  //kmerHolder* kh = *khp;
  seqCollection* sc = *scp;
  if (sc->trace){
    char* seq = getTraceSeq(khp, &sc->trace);
    printf("%s\n", seq);
    free(seq);
    //printKcLL(kh->ms, sc->trace);
  }
}

void printSeqCollection(kmerHolder** khp, seqCollection** scp){
  seqCollection* sc = *scp;
  uint32_t i = 0;
  if (!sc) printf("No traces\n");
  while (sc){
    i++;
    if (sc->delme) D_(0, "delme\n");
    printf(">trace%"PRIu32"_pos%"PRIu32"\n", sc->trace->ntid, sc->posInTrace);
    printSeq(khp, &sc);
    sc = sc->down;
  }
}

void doNothing(kmerHolder** khp, uint32_t t){
  return;
}

seqCollection* searchSeq(kmerHolder** khp, char* seq){
  kmerHolder* kh = *khp;
  memstruct* ms = kh->ms;
  uint32_t lseq = (uint32_t) strlen(seq);
  seqCollection* result = NULL;
  D_(1, "Resetting trace\n");
  resetTrace(khp);
  uint32_t prevPos = NOKMER;
  for (uint32_t i = 0; i < lseq; i++){
    char b = seq[i];
    uint32_t cPos = updateKmer(khp, &b, doNothing);
    if (cPos != NOKMER && prevPos != NOKMER){
      kmerConnector* kc = ms->kmerArray[prevPos];
      while (kc && kc->dest != cPos){
        kc = kc->next;
      }
      if (!kc) return NULL;
      if (result){
        seqCollection* tsc = result;
        while (tsc){
          kcLL* tr = tsc->trace;
          tsc->posInTrace++;
          tIdList* x = _getTrace(&kc->idflags, tr->ntid, tsc->posInTrace);
          if (x){
            kcpush(&tr, &kc, x->trace.n, 0);
          }
          else{
            tsc->delme = true;
          }
          tsc = tsc->down;
        }
      }
      else{
        tIdList* tmp = kc->idflags;
        while (tmp){
          tIdList* eachPos = tmp->posInTrace;
          while (eachPos){
            kcLL* tkcll = NULL;
            kcpush(&tkcll, &kc, tmp->trace.n, cPos);
            pushSeq(&result, &tkcll, eachPos->trace.n);
            eachPos = eachPos->next;
          }
          tmp = tmp->next;
        }
      }
      prevPos = cPos;
    }
    else{
      prevPos = cPos;
    }
  }
  pruneSeqCollection(&result);
  return result;
}


#endif // KMERREAD_H_INCLUDED
