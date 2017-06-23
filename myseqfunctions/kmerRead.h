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
          result->posInTrace = i;
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
  pos2seq(&kh->ms, kcll->posInTrace, result);
  uint32_t j = (uint32_t) ((kh->kmerSize) - 1);
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
  kmerHolder* kh = *khp;
  seqCollection* sc = *scp;
  if (sc->trace){
    /*char* seq = getTraceSeq(khp, &sc->trace);
    printf("%s\n", seq);
    free(seq);*/
    printKcLL(kh->ms, sc->trace);
  }
}

void printSeqCollection(kmerHolder** khp, seqCollection** scp){
  seqCollection* sc = *scp;
  uint32_t i = 0;
  while (sc){
    i++;
    printf(">trace%lu\n", (LUI) i);
    printSeq(khp, &sc);
    sc = sc->down;
  }
}

seqCollection* searchSeq(kmerHolder** khp, char* seq){
  kmerHolder* kh = *khp;
  uint32_t lseq = (uint32_t) strlen(seq);
  seqCollection* result = NULL;
  D_(1, "Resetting trace\n");
  kh->cVal = 0;
  kh->posInMemBuffer = 0;
  for (uint32_t i = 0; i < lseq; i++){
    char b = seq[i];
    uint8_t baseVal = base(&kh->ms, &b);
    D_(2, "Base %c, val %u\n", (int) b, baseVal);
    if (!isBase(&kh->ms, seq)){
      D_(0, "Sequence contains non-standard character %c. Results may not be reliable\n", b);
      return result;
    }
    uint32_t ipos = kh->posInMemBuffer - (uint32_t) kh->kmerSize;
    if (kh->posInMemBuffer >= (uint32_t) kh->kmerSize){ //kmer filled
      uint32_t order = power(kh->ms->bi->nbases, kh->kmerSize - 1);
      D_(2, "%" PRIu32 "\n", ipos);
      unsigned char oldval = kh->buffer[ipos];
      D_(2, "cVal: %" PRIu32 "\n", kh->cVal);
      kmerConnector* kc = kh->ms->kmerArray[kh->cVal];
      if (result){
        seqCollection* tmpsc = result;
        while (tmpsc){
          if (!tmpsc->delme){
            kcLL* kcT = tmpsc->trace;
            LISTTYPE i = kcT->ntid;
            LISTTYPE p = kcT->posInTrace;
            kcT->posInTrace++;
            tIdList* point = _getTrace(&kcT->kc->idflags, i, p);
            if (point){
              kcpush(&kcT, &kc, i, p);
            }
            else{ //Delete sc
              tmpsc->delme = true;
            }
          }
          tmpsc = tmpsc->down;
        }
      }
      else{
        if (!kc) return NULL;
        tIdList* tmpid = kc->idflags;
        while (tmpid){
          LISTTYPE ntid = tmpid->trace.n;
          D_(2, "%"PRIu32"\n", ntid);
          tIdList* tmppos = tmpid->posInTrace;
          while (tmppos){
            kcLL* starting = NULL;
            LISTTYPE pos = tmppos->trace.n;
            kcpush(&starting, &kc, ntid, pos);
            pushSeq(&result, &starting);
            tmppos = tmppos->next;
          }
          tmpid = tmpid->next;
        }
      }
      kh->cVal -=  order * (uint32_t) oldval;
    }
    kh->cVal *= (uint32_t) kh->ms->bi->nbases;
    kh->cVal += (uint32_t) baseVal;
    kh->buffer[kh->posInMemBuffer] = baseVal;
    kh->posInMemBuffer++;
    if (kh->posInMemBuffer >= (KMERBUFFERSIZE - 1)){ // Rewind
      //printf("Rewinding\n");
      memmove(kh->buffer, kh->buffer + ipos + 1, sizeof(char) * (kh->kmerSize));
      kh->posInMemBuffer = (uint32_t) kh->kmerSize;
    }
  }
  return result;
}


#endif // KMERREAD_H_INCLUDED
