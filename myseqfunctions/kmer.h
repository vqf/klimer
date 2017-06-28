#ifndef KMER_H_INCLUDED
#define KMER_H_INCLUDED


#define ENCLOSE(a) printf("\n---before, line %d\n", __LINE__); a; printf("\nafter, line %d---\n", __LINE__);
#define SEQ(a, b) char* _seq = (char*) calloc(12, sizeof(char)); pos2seq(a, b->dest, _seq); printf("%s\n", _seq); free(_seq);

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <inttypes.h>
#include "neatHtml.h"
#include "traceIdList.h"
#include "seqCollection.h"

#define KMERBUFFERSIZE 0xFFFF
#define NOKMER 0xFFFFFFFF
#define LUI long unsigned int
#define HISTMAX 0xFFFFF

#ifndef BOOL
#define BOOL
typedef enum { false, true } bool;
#endif /* BOOL */


// Set, unset and check flags in traces
// a -> tId object
// b -> flag to operate

uint32_t _UID = 0;

//#define addPosInTrace(a, b) printf("Line %d, cpos %lu\n", __LINE__, (LUI) cposInTrace); addPosInTrace(a, b)


// FLAGS

//

typedef struct tmpKcVessel{
  tId* dest;
  kmerConnector* loop;
  struct tmpKcVessel* next;
} tmpKcVessel;

typedef struct bs{
  uint8_t* baseVal; //Fast conversion of bases to values (see initMemory)
  uint8_t* valBase; //Fast conversion from values to bases
  uint8_t nbases; // Usually 4, or 5 with N
} baseInfo;



typedef struct st{
  traceLL* inUse;       // Which lists are in use
  tSet* traceSet; // Keep tab of trace IDs the current read belongs to
  LISTTYPE cId;      // Current id for kmerConnector
  LISTTYPE posInRead;
  kcLL* trace;      // Keeps track of the related kmers in the current trace
  uint32_t firsInTrace;  // Can keep tab of current when following trace
  uint32_t current; // Last kmer read
  bool start;       // No kc has still been filled
  bool isFirst;     // Is it the first in a trace?
  uint8_t addExistingTraceStatus; // See _existingTrac
} msStatus;


typedef struct ms{ //Used to store all kmer combinations
  baseInfo* bi;
  uint32_t nPos; // Number of positions
  uint32_t nBytes; // Number of bytes in kmerArray
  msStatus* status;
  uint8_t kmerSize;
  kmerConnector** kmerArray;
} memstruct;


typedef struct kh{  //Used to read kmers
  unsigned char* buffer;
  memstruct* ms;
  uint32_t posInMemBuffer;
  uint8_t kmerSize;
  uint8_t nBases;
  uint32_t cVal;
  bool canonic;  // Has it been standardized?
} kmerHolder;


void resetTrace(kmerHolder**);
void resetKmer(kmerHolder**);
void summarize (memstruct*);
void resetKcLL(kcLL**);
void printKmerConnector(kmerConnector*, char*);
void printKmerConnectors(memstruct*, uint32_t);
void printKc(memstruct*);
void printRead(memstruct*);
void printKcLL(memstruct*, kcLL*);
void destroyTIdList(tIdList**);
kmerConnector* newKmerConnector(uint32_t);
kmerConnector* getConnector(memstruct**, uint32_t, uint32_t);
void delConnector(kmerConnector**);
void cleanTraceStatus(kmerHolder**);
kcLL* followTrace(kmerHolder**, uint32_t, LISTTYPE, LISTTYPE);
void printSeqCollection(kmerHolder**, seqCollection**);


// SIGNAL HANDLING
kmerHolder* khDebug = NULL;
kmerConnector* kcDebug = NULL;
LISTTYPE iDebug = 0;
LISTTYPE posDebug = 0;
char* getTraceSeq(kmerHolder**, kcLL**);
kmerConnector* nextKc(memstruct**, kmerConnector**, tIdList**, LISTTYPE);
static void sig_handler(int signo) {
    printf("Interactive attention signal caught.\n");
    printf("Trace %lu\n", (LUI) iDebug);
    DEBUG++;
    STOP("Debug: %d\n", DEBUG);
    /*
    kmerConnector* kcptr = kcDebug;
    kcLL* toseq = NULL;
    while (kcptr){
      kcpush(&toseq, &kcptr, iDebug, posDebug);
      tIdList* tmp = _getTrace(&kcptr->idflags, iDebug, posDebug);
      kcptr = nextKc(&khDebug->ms, &kcptr, &tmp, posDebug);
      posDebug++;
    }
    char* seq = getTraceSeq(&khDebug, &toseq);
    printf("%s\n", seq);
    printTIdList(kcDebug->idflags); printf("\n");*/
}

// SIGNAL HANDLING

uint32_t lengthKcll(kcLL** lp){
  kcLL* l = *lp;
  uint32_t result = 0;
  while(l){
    result++;
    l = l->next;
  }
  return result;
}

tmpKcVessel* newTmpKcVessel(tId* dst, kmerConnector** lp, tmpKcVessel** nxt){
  tmpKcVessel* tonext = NULL;
  if (nxt) tonext = *nxt;
  tmpKcVessel* result = (tmpKcVessel*) calloc(1, sizeof(tmpKcVessel));
  result->dest = dst;
  result->loop = *lp;
  result->next = tonext;
  return result;
}

void destroytmpKcVessel(tmpKcVessel** kvp){
  tmpKcVessel* kv = *kvp;
  while (kv){
    tmpKcVessel* tmp = kv->next;
    free(kv);
    kv = tmp;
  }
}

msStatus* initStatus(){
  msStatus* result = (msStatus*) calloc(1, sizeof(msStatus));
  result->trace = NULL;
  result->addExistingTraceStatus = 0;
  result->inUse = NULL;
  result->traceSet = NULL;
  result->firsInTrace = NOKMER;
  result->current = NOKMER;
  result->cId = 0;
  result->posInRead = 0;
  return result;
}

memstruct* initFourBase(){  /*
  * Init MM_baseVal for fast base value retrieval
  * MM_baseVal contains 0xFF elements, all except
  * the corresponding ACGT with a value of 0xFF.
  * Lowercase and uppercase ACGT have 0x00, 0x01, ...
  */

  memstruct* result = (memstruct*) malloc(sizeof(memstruct));
  result->bi = (baseInfo*) calloc(1, sizeof(baseInfo));
  result->bi->baseVal = (uint8_t*) malloc(sizeof(uint8_t) * (0x100));
  result->bi->valBase = (uint8_t*) calloc(5, sizeof(uint8_t));
  assert(result->bi->baseVal);
  assert(result->bi->valBase);
  result->bi->valBase[0x00] = 0x41;
  result->bi->valBase[0x01] = 0x43;
  result->bi->valBase[0x02] = 0x47;
  result->bi->valBase[0x03] = 0x54;
  int i = 0;
  for(i = 0; i < 0xFF; i++){
    result->bi->baseVal[i] = 0xFF;
  }
  result->bi->baseVal[0x61] = 0x00; result->bi->baseVal[0x41] = 0x00; //aA
  result->bi->baseVal[0x63] = 0x01; result->bi->baseVal[0x43] = 0x01; //cC
  result->bi->baseVal[0x67] = 0x02; result->bi->baseVal[0x47] = 0x02; //gG
  result->bi->baseVal[0x74] = 0x03; result->bi->baseVal[0x54] = 0x03; //tT
  result->bi->nbases = 4;
  // Init other fields
  result->status = initStatus();
  return result;
}

memstruct* initFiveBase(){
  /*
  * Init MM_baseVal for fast base value retrieval
  * MM_baseVal contains 0xFF elements, all except
  * the corresponding ACGTN with a value of 0xFF.
  * Lowercase and uppercase ACGTN have 0x00, 0x01, ...
  */

  memstruct* result = (memstruct*) malloc(sizeof(memstruct));
  result->bi = (baseInfo* ) calloc(1, sizeof(baseInfo));
  result->bi->baseVal = (uint8_t*) malloc(sizeof(uint8_t) * (0x100));
  result->bi->valBase = (uint8_t*) calloc(5, sizeof(uint8_t));
  result->bi->valBase[0x00] = 0x41;
  result->bi->valBase[0x01] = 0x43;
  result->bi->valBase[0x02] = 0x47;
  result->bi->valBase[0x03] = 0x54;
  result->bi->valBase[0x04] = 0x4E;
  int i = 0;
  for(i = 0; i < 0xFF; i++){
    result->bi->baseVal[i] = 0xFF;
  }
  result->bi->baseVal[0x100] = 0x00;
  result->bi->baseVal[0x61] = 0x00; result->bi->baseVal[0x41] = 0x00; //aA
  result->bi->baseVal[0x63] = 0x01; result->bi->baseVal[0x43] = 0x01; //cC
  result->bi->baseVal[0x67] = 0x02; result->bi->baseVal[0x47] = 0x02; //gG
  result->bi->baseVal[0x74] = 0x03; result->bi->baseVal[0x54] = 0x03; //tT
  result->bi->baseVal[0x6E] = 0x04; result->bi->baseVal[0x4E] = 0x04; //nN
  result->bi->nbases = 5;
  result->status = initStatus();
  return result;
}

void destroyCircular(void** toDest){
  void* p = *toDest;
  if (p){
    kcLL* k = (kcLL*) p;
    resetKcLL(&k);
  }
}

void destroyMs(memstruct** msp){
  memstruct* ms = *msp;
  resetKcLL(&ms->status->trace);
  destroyTSet(&ms->status->traceSet);
  destroyTraceLL(&ms->status->inUse);
  free(ms->status);
  free(ms->bi->valBase);
  free(ms->bi->baseVal);
  free(ms->bi);
  // The big one: kmerArray and connected kcs
  kmerConnector* tmp = NULL;
  //fprintf(stderr, "Cleaning ms array, located at %p...\n", *(ms->kmerArray));
  for (int i = 0; i < ms->nPos; i++){
    tmp = ms->kmerArray[i];
    delConnector(&tmp);
  }
  free(ms->kmerArray);
  free(ms);
  ms = NULL;
}

void destroyKh(kmerHolder** khp){
  kmerHolder* kh = *khp;
  free(kh->buffer);
  destroyMs(&kh->ms);
  free(*khp);
}

bool isBase (memstruct** msp, char* b){
  memstruct* ms = *msp;
  uint8_t val = ms->bi->baseVal[(uint8_t) b[0]];
  bool result = false;
  if (val >= 0 && val < ms->bi->nbases) result = true;
  return result;
}

uint8_t base (memstruct** msp, char* b){
  memstruct* ms = *msp;
  uint8_t result = (uint8_t) ms->bi->baseVal[(uint8_t) b[0]];
  return result;
}

uint32_t power(uint8_t base, uint8_t exp) {
  if (exp == 0)
    return 1;
  else if (exp % 2)
    return base * power(base, exp - 1);
  else {
    int temp = power(base, exp / 2);
    return temp * temp;
  }
}

void setKmerArray(kmerHolder** kp, uint8_t wlength){
  /*
   *Initializes and zeroes a block of memory with the values for each
   *k-mer.
  */
  memstruct* ms = (*kp)->ms;
  ms->kmerSize = wlength;
  ms->nPos = power(ms->bi->nbases, wlength);
  cleanTraceStatus(kp);
  uint64_t nb = sizeof(kmerConnector) * ms->nPos;
  ms->nBytes = nb;
  //DIE("%" PRIu32 " times %" PRIu32 " equals %" PRIu64 "\n", (uint32_t) (sizeof(kmerConnector)), (uint32_t) ms->nPos, nb);
  D_(1, "Allocating %" PRIu64 " bytes...\n", nb);
  ms->kmerArray = (kmerConnector**) calloc(ms->nPos, sizeof(kmerConnector*));
  D_(1, "Done\n");
}

uint32_t seq2pos (memstruct* ms, char* seq){
  /*
   * The length of the seq is set to ms kmer length
   */
  uint32_t pos = NOKMER;
  uint16_t sl = strlen(seq);
  uint16_t i;
  for (i = 0; i < sl; i++){
    char* tb = seq + i;
    if (isBase(&ms, tb)){
      if (pos == NOKMER) pos = 0;
      pos = ms->bi->nbases * pos + base(&ms, tb);
    }
    else{
      fprintf(stderr, "Unknown base %s\n", tb);
    }
  }
  return pos;
}

void pos2seq (memstruct** msp, uint32_t pos, char* result){
  memstruct* ms = *msp;
  int8_t i = ms->kmerSize - 1;
  int8_t j = 0;
  while (pos > 0){
    uint8_t r = (uint8_t) (pos % ms->bi->nbases);
    pos = pos / ms->bi->nbases;
    result[i] = ms->bi->valBase[r];
    i--;
  }
  for (j = i; j >= 0; j--){
    result[j] = ms->bi->valBase[0x00];
  }
}

uint16_t _getVal(memstruct* ms, uint32_t pos){
  uint16_t result = ms->kmerArray[pos]->n;
  return result;
}

uint16_t getVal(memstruct* ms, char* kmer){
  uint32_t kpos = seq2pos(ms, kmer);
  uint16_t result = _getVal(ms, kpos);
  return result;
}

/*
 *Simple addCount, no relationships
 */
void addCount(kmerHolder** kp, uint32_t pos){
  memstruct *ms = (*kp)->ms;
  kmerConnector* result = ms->kmerArray[pos];
  if (pos < 0 || pos > ms->nPos){
    //fprintf(stderr, "Too much\n");
    return;
  }
  if (result){
    if (result->n <= 0xFFFF) result->n++;
  }
  else{
    ms->kmerArray[pos] = newKmerConnector(0);
    ms->kmerArray[pos]->n++;
  }
}

void scEliminateDuplicates(kmerHolder** khp, seqCollection** scp){
  seqCollection* tsc = *scp;
  while (tsc){
    tsc->seq = getTraceSeq(khp, &tsc->trace);
    tsc = tsc->down;
  }
  seqCollection* sc = *scp;
  while (sc){
    seqCollection* f = sc->down;
    while (f){
      E_(1, printf("%s\t%s\n", sc->seq, f->seq));
      if (EQ(sc->seq, f->seq)){
        f->delme = true;
      }
      f = f->down;
    }
    sc = sc->down;
  }
}

seqCollection* exploreSeqs(kmerHolder** khp, uint32_t start, uint32_t l){
  kmerHolder* kh = *khp;
  memstruct* ms  = kh->ms;
  seqCollection* result = NULL;
  kmerConnector* kc = ms->kmerArray[start];
  if (!kc) return result;
  uint32_t i = 0;
  while (i < (l - kh->kmerSize)){
    if (result){
      seqCollection* tsc = result;
      while (tsc){
        kcLL* tr = tsc->trace;
        kmerConnector* pkc = tr->last->kc;
        tIdList* tt = _getTrace(&pkc->idflags, tr->ntid, tsc->posInTrace);
        if (!tt){
          printf("%" PRIu32 "\t%" PRIu32 "\n", tr->ntid, tsc->posInTrace);
          printTIdList(pkc->idflags);
          X_
        }
        kmerConnector* tkc = nextKc(&ms, &pkc, &tt, tsc->posInTrace);
        if (tkc){
          kcpush(&tr, &tkc, tr->ntid, 0);
          tsc->posInTrace++;
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
          E_(2, printf("%" PRIu32 "\t%" PRIu32 "\n", tmp->trace.n, eachPos->trace.n));
          kcLL* tkcll = NULL;
          kcpush(&tkcll, &kc, tmp->trace.n, start);
          pushSeq(&result, &tkcll, eachPos->trace.n);
          eachPos = eachPos->next;
        }
        tmp = tmp->next;
      }
    }
    i++;
  }
  scEliminateDuplicates(khp, &result);
  pruneSeqCollection(&result);
  return result;
}

uint32_t nReads(kcLL** kclp){
  kcLL* kcl = *kclp;
  uint32_t result = (uint32_t) kcl->kc->n;
  while (kcl){
    tIdList* tid = _getTrace(&kcl->kc->idflags, kcl->ntid, 0);
    uint32_t n = (uint32_t) kcl->kc->n;
    if (n < result) result = n;
    kcl = kcl->next;
  }
  return result;
}

void printHistogram (kmerHolder** khp, uint8_t kl){
  kmerHolder* kh = *khp;
  memstruct* ms = kh->ms;
  uint16_t hmax = 0;
  uint32_t* hist = (uint32_t*) calloc(HISTMAX, sizeof(uint32_t));
  for (uint32_t i = 0; i < ms->nPos; i++){
    kmerConnector* kc = ms->kmerArray[i];
    if (kc){
      seqCollection* fromHere = exploreSeqs(khp, i, kl);
      //printSeqCollection(khp, &fromHere);
      while (fromHere){
        kcLL* tr = fromHere->trace;
        uint32_t v = nReads(&tr);
        printf("%s\t%"PRIu32"\n", fromHere->seq, v);
        fromHere = fromHere->down;
      }
      X_
      uint16_t n = kc->n;
      if (n < HISTMAX){
        if (n > hmax) hmax = n;
        hist[n]++;
      }
    }
  }
  for (uint32_t i = 1; i <= hmax; i++){
    printf("%" PRIu32 "\t%" PRIu32 "\n", i, hist[i]);
  }
  free(hist);
}

/*
 *With relationships
 */


void nextId(memstruct** msp){
  memstruct* ms = *msp;
  if (ms->status->cId == MAXLISTTYPE){
    D_(0, "%s\n", "Overflow in cId");
    return;
  }
  ms->status->cId++;
}


/*DEBUG*/
int measureTrace(kcLL** llp){
  kcLL* ll = *llp;
  int result = 0;
  while (ll->next){
    result++;
    ll = ll->next;
  }
  return result;
}
void debugTrace(kcLL** llp){
  kcLL* ll = *llp;
  if (!ll) return;
  printf("---\n%.8x -> \n---\n", ll->kc->uid);
  while (ll->next){
    ll = ll->next;
    printf("---\n%.8x -> \n---\n", ll->kc->uid);
  }
}
/**/

kmerConnector* getKcWithTId(memstruct** msp, uint32_t pos, LISTTYPE tid, LISTTYPE posInTrace, tIdList** tp){
  //Returns the kc and sets tp to the corresponding trace
  memstruct* ms = *msp;
  if (pos == NOKMER) return NULL;
  kmerConnector* kc = ms->kmerArray[pos];
  kmerConnector* result = NULL;
  while (kc && !result){
    tIdList* tt = _getTrace(&kc->idflags, tid, posInTrace);
    if (tt){
      result = kc;
      *tp = tt;
    }
    kc = kc->next;
  }
  return result;
}



kmerConnector* nextKc(memstruct** msp, kmerConnector** kcp, tIdList** p, LISTTYPE pos){
  kmerConnector* kc = *kcp;
  D_(2, "From %.8x\n", (unsigned int) kc->uid);
  kmerConnector* result = NULL;
  result = getKcWithTId(msp, kc->dest, (*p)->trace.n, pos + 1, p);
  return result;
}
//#define nextKc(a,b,c,d) nextKc(a,b,c,d); printf("%d\n", __LINE__);

bool isLastInTrace(tIdList** tp, LISTTYPE pos){
  tIdList* t = *tp;
  tIdList* myl = isInTIdList(&t->posInTrace, pos);
  if (IS(myl, LAST_IN_TRACE)){
    return true;
  }
  return false;
}




void cleanTraceStatus(kmerHolder** khp){
  memstruct* ms = (*khp)->ms;
  //ms->status->cId = 0;
  ms->status->posInRead = 0;
  destroyTSet(&ms->status->traceSet);
  ms->status->addExistingTraceStatus = 0;
  ms->status->start   = true;
  ms->status->isFirst = true;
  ms->status->current = NOKMER;
  resetKmer(khp);
  resetKcLL(&ms->status->trace);
  ms->status->trace = NULL;
}


void shiftPosInTrace(memstruct** msp, kmerConnector** kcp, tIdList** tp, LISTTYPE n){
  // returns number of kcs shifted
  tIdList* t = *tp;
  tIdList* tmp = t;
  D_(2, "Shifting trace %lu by %lu\n", (LUI) t->trace.n, (LUI) n);
  kmerConnector* kc = *kcp;
  LISTTYPE cPos = 1;
  if (t){
    LISTTYPE i = t->trace.n;
    traceVessel* posToShift = NULL;
    while (kc){
      tIdList* toShift = tmp;
      if (toShift){
        tIdList* ps = _getPos(&toShift, i, cPos);
        pushTraceInVessel(&posToShift, &ps);
      }
      else{
        D_(0, "Error, lacking trace %lu, pos %lu\n", (LUI) i, (LUI) cPos);
      }
      kc = nextKc(msp, &kc, &tmp, cPos);
      cPos++;
    }
    traceVessel* tmp = posToShift;
    while (tmp){
      tmp->tidl->trace.n += n;
      tmp = tmp->next;
    }
    destroyTraceVessel(&posToShift);
  }
}


void resetTrace(kmerHolder** kp){
  D_(2, "Resetting trace\n");
  memstruct* ms = (*kp)->ms;
  kcLL** tmpp = &ms->status->trace;
  kcLL* tmp = *tmpp;
  bool newTrace = true;
  uint32_t tlen = lengthKcll(tmpp);
  D_(2, "%"PRIu32"\t%"PRIu8"\n", tlen, (*kp)->kmerSize);
  if (tlen < ((uint32_t) (*kp)->kmerSize)){
    cleanTraceStatus(kp);
    return;
  }
  if (ms->status->traceSet){
    newTrace = false;
    // If several tSets refer to the same trace, only one
    // will be incorporated.
    tIdList* used = NULL;
    tSet* tmpTSet = ms->status->traceSet;
    while (tmpTSet){
      if (isInTIdListNotInUse(&used, tmpTSet->traceId)){
        tmpTSet->conflict = true;
      }
      else{
        insertInTIdList(&used, tmpTSet->traceId);
      }
      tmpTSet = tmpTSet->next;
    }
    destroyTIdList(&used);
  }
  LISTTYPE cposInTrace = 1;  // Used for newTrace
  if (newTrace){
    nextId(&ms);
    D_(1, "This is a new trace (%lu)\n", (LUI) ms->status->cId);
    tIdList* fpos = NULL;
    tIdList* lpos = NULL;
    while (tmp){
      tIdList* tl = addTrace(&tmp->kc->idflags, ms->status->cId);
      tIdList* cPosTid = addPosInTrace(&tl, cposInTrace);
      if (fpos){
        lpos = cPosTid;
      }
      else{
        fpos = cPosTid;
      }
      cposInTrace++;
      tmp = tmp->next;
    }
    if (ms->status->trace){
      setAsFirst(&ms->status->trace->kc->idflags, ms->status->cId, 1);
      setAsLast(&ms->status->trace->last->kc->idflags, ms->status->cId, cposInTrace - 1);
      //SET(fpos, FIRST_IN_TRACE); SET(lpos, LAST_IN_TRACE); // Avoid unnecessary calls to _getTrace
      SET(fpos, FREAD); SET(lpos, LREAD);
    }
  }
  else{
    tSet* eachTSet = ms->status->traceSet;
    E_(2, printTSet(ms->status->traceSet));
    while (eachTSet){
      cposInTrace = 1;
      tIdList* rFst = NULL;
      tIdList* rLst = NULL;
      if (!eachTSet->conflict){
        tmp = ms->status->trace;
        LISTTYPE cId = eachTSet->traceId;
        LISTTYPE posInRead = 1; // Only used for type 2 and 3
        if (eachTSet->type == 0){
          LISTTYPE nId = eachTSet->traceId;
          cposInTrace = eachTSet->startingPos;
          D_(1, "Inside trace (%" PRIu32 ")\n", nId);
          while (tmp){
            tIdList* tl = addTrace(&tmp->kc->idflags, nId);
            tIdList* cPos = addPosInTrace(&tl, cposInTrace);
            if (rFst){
              rLst = cPos;
            }
            else{
              rFst = cPos;
            }
            cposInTrace++;
            tmp = tmp->next;
          }
        }
        cposInTrace = eachTSet->startingPos;
        if (eachTSet->type == 1 || eachTSet->type == 3){ // Extending up
          D_(1, "Extending up trace %" PRIu32 "\n", eachTSet->traceId);
          LISTTYPE counter = 1;
          while (counter < eachTSet->upShift){
            counter++;
            tmp = tmp->next;
          }
          tIdList* fst = _getTrace(&tmp->kc->idflags, cId, 1);
          if (!fst){
            D_(0, "Warning: There should be a starting trace here, but I cannot find it\n");
            printTIdList(tmp->kc->idflags);
            X_
          }
          if (isFirst(&fst)){
            unsetAsFirst(&fst, cId, 1);
            shiftPosInTrace(&ms, &tmp->kc, &fst, eachTSet->upShift - 1);
          }
          else{
            char* s = getTraceSeq(kp, tmpp);
            printf ("%s\n", s);
            D_(0, "Warning, there should be a first_in_trace in trace %" PRIu32 \
               " pos %" PRIu32 "\n", eachTSet->traceId, cposInTrace);
            printTIdList(tmp->kc->idflags);
          }
          tmp = ms->status->trace;
          while (cposInTrace < eachTSet->upShift){
            tIdList* tl = addTrace(&tmp->kc->idflags, cId);
            tIdList* cPos = addPosInTrace(&tl, cposInTrace);
            if (rFst){
              rLst = cPos;
            }
            else{
              rFst = cPos;
            }
            cposInTrace++;
            posInRead++;
            tmp = tmp->next;
          }
          setAsFirst(&ms->status->trace->kc->idflags, cId, 1);
          if (eachTSet->type == 3){
            eachTSet->type = 2; //Defer to type 2
          }
        }
        if (eachTSet->type == 2){
          while(posInRead < eachTSet->dnShift){
            tIdList* tl = addTrace(&tmp->kc->idflags, cId);
            tIdList* cPos = addPosInTrace(&tl, cposInTrace);
            if (rFst){
              rLst = cPos;
            }
            else{
              rFst = cPos;
            }
            cposInTrace++;
            posInRead++;
            tmp = tmp->next;
          }
          tIdList* lst = _getTrace(&tmp->kc->idflags, cId, cposInTrace );
          if (!lst){
            D_(0, "Cannot find previous last\n"); X_;
          }
          if (isLast(&lst)){
            unsetAsLast(&lst, cId, cposInTrace);
          }
          else{
            D_(0, "Warning, this one should be last!\n");
          }
          while(tmp){
            tIdList* tl = addTrace(&tmp->kc->idflags, cId);
            tIdList* cPos = addPosInTrace(&tl, cposInTrace);
            if (rFst){
              rLst = cPos;
            }
            else{
              rFst = cPos;
            }
            cposInTrace++;
            tmp = tmp->next;
          }
          setAsLast(&ms->status->trace->last->kc->idflags, cId, cposInTrace - 1);
        }
        if (eachTSet->type > 3){
          D_(0, "At this point, type must be 0, 1, 2 or 3, not %u\n", eachTSet->type);
          X_;
        }
        // finish it
        while(tmp){
          tIdList* tl = addTrace(&tmp->kc->idflags, cId);
          tIdList* cPos = addPosInTrace(&tl, cposInTrace);
          if (rFst){
            rLst = cPos;
          }
          else{
            rFst = cPos;
          }
          cposInTrace++;
          tmp = tmp->next;
        }
      }
      if (rFst) SET(rFst, FREAD); if (rLst) SET(rLst, LREAD);
      eachTSet = eachTSet->next;
    }
  }
  cleanTraceStatus(kp);
}
//#define resetTrace(a) D_(0, "Calling resetTrace from %s, line %d\n", __FUNCTION__, __LINE__); resetTrace(a);




/*
 * kmerConnector functions
 */

uint32_t rnd32(){
  uint32_t r = 0;
  int h = RAND_MAX >> 1;
  for (int i = 0; i < 32; i++){
    int n = rand();
    int b = 0;
    if (n > h){
      b = 1;
    }
    r = r << 1;
    r += b;
  }
  //printf("yo %.8x\n", r);
  return r;
}

kmerConnector* newKmerConnector(uint32_t to){
  kmerConnector* result = (kmerConnector*) calloc(1, sizeof(kmerConnector));
  result->dest    = to;
  result->uid     = GETUID;
  result->n       = 0;
  result->idflags = NULL;
  result->next    = NULL;
  return result;
}

void delConnector(kmerConnector** kcp){
  kmerConnector* tmp = *kcp;
  while (tmp){
    kmerConnector* nxt = tmp->next;
    destroyTIdList(&tmp->idflags);
    free(tmp);
    tmp = nxt;
  }
  *kcp = NULL;
}


LISTTYPE _getLastPos(tIdList** lp){
  LISTTYPE result = 0;
  tIdList* l = *lp;
  if (!l) return result;
  tIdList* eachPos = l->posInTrace;
  while (eachPos){
    if (IS(eachPos, LAST_IN_TRACE)){
      result = eachPos->trace.n;
    }
    eachPos = eachPos->next;
  }
  return result;
}

void _existingTrace(memstruct** msp, kmerConnector** kcp, uint32_t from){
  memstruct* ms = *msp;
  kmerConnector* kc = *kcp;
  ms->status->posInRead++;
  // Remove incompatible traces
  tSet* eachTSet = ms->status->traceSet;
  traceVessel* lst = NULL;
  if (eachTSet) D_(2, "-Existing trace\n======\n");
  if (kc->idflags) lst = traceLast(&kc->idflags);
  while (eachTSet){
    bool getNext = true;
    if (!eachTSet->fixed){
      D_(2, "tSet [%lu, %lu]\n", (LUI) eachTSet->traceId, (LUI) eachTSet->currentPos);
      if (eachTSet->type != 2){
        tIdList* thisTrace = NULL;
        if (kc->idflags) thisTrace = _getTrace(&kc->idflags, eachTSet->traceId, eachTSet->currentPos+1);
        if (thisTrace){
          if (lst){
            D_(2, "Compatible\n");
            tIdList* islst = _getTrace(&lst->tidl, eachTSet->traceId, eachTSet->currentPos+1);
            LISTTYPE lastPos = _getLastPos(&islst);
            if (islst && lastPos == (eachTSet->currentPos + 1)){
              eachTSet->type += 2;
              eachTSet->dnShift = ms->status->posInRead;
              eachTSet->fixed = true;
              D_(2, "Ends trace at pos %lu\n", (LUI) ms->status->posInRead);
            }
          }
          eachTSet->currentPos++;
        }
        else{
          D_(2, "Incompatible with id %lu pos %lu\n", (LUI) eachTSet->traceId, (LUI) (eachTSet->currentPos));
          tSet* delme = eachTSet;
          eachTSet = eachTSet->next;
          exciseTSet(&ms->status->traceSet, &delme);
          E_(2, printTSet(eachTSet));
          getNext = false;
        }
      }
    }
    if (getNext){
      eachTSet = eachTSet->next;
    }
  }
  destroyTraceVessel(&lst);
  // Add traces if necessary
  if (kc->idflags){
    tIdList* eachTrace = kc->idflags;
    while (eachTrace){
      if (ms->status->isFirst){
        D_(2, "Inside, not extending");
        tIdList* dPos = eachTrace->posInTrace;
        while (dPos){
          D_(2, "Trace %lu, pos %lu\n", (LUI) eachTrace->trace.n, (LUI) dPos->trace.n);
          unshiftTSet(&ms->status->traceSet, eachTrace->trace.n, dPos->trace.n);
          ms->status->traceSet->type = 0;
          dPos = dPos->next;
        }
      }
      else if (isFirst(&eachTrace)){
        D_(2, "Might extend up trace %lu\n", (LUI) eachTrace->trace.n);
        unshiftTSet(&ms->status->traceSet, eachTrace->trace.n, 1);
        ms->status->traceSet->type = 1;
        ms->status->traceSet->upShift = ms->status->posInRead;
      }
      eachTrace = eachTrace->next;
    }
  }
  E_(2, printTSet(ms->status->traceSet));
  D_(2, "\n======\n");
  if (ms->status->isFirst) ms->status->isFirst = false;
}


kmerConnector* getConnector(memstruct** msp, uint32_t from, uint32_t to){
  memstruct *ms = *msp;
  kmerConnector* result = ms->kmerArray[from];
  D_(2, "Getting connector from %lu to %lu\n", (LUI) from, (LUI) to);
  if (!result){ // This connector did not exist
    D_(2, "This is a new one\n");
    ms->kmerArray[from] = newKmerConnector(to);
    result = ms->kmerArray[from];
    _existingTrace(msp, &result, from);
  }
  else{
    while (result->dest != to && result->next){
      result = result->next;
    }
    if (result->dest == to){
      D_(2, "Oh, you mean this one\n");
      _existingTrace(msp, &result, from);
    }
    else{ // This connector did not exist
      kmerConnector* nkc = newKmerConnector(to);
      result->next = nkc;
      result = result->next;
      D_(2, "New destination\n");
      _existingTrace(msp, &result, from);
    }
  }
  //ms->status->cId = maxInList(ms->status->traceSet);
  kcpush(&ms->status->trace, &result, 0, 0);
  return result;
}

void addRelationship(kmerHolder** kp, uint32_t to){
  memstruct* ms = (*kp)->ms;
  if (ms->status->start){
    ms->status->start = false;
    ms->status->current = to;
    ms->status->firsInTrace = to;
    return;
  }
  uint32_t from = ms->status->current;
  if (from == NOKMER){
    resetTrace(kp);
  }
  else{
    kmerConnector* kc = getConnector(&ms, from, to);
    //printKmerConnectors(ms, from);
    kc->n++;
  }
  ms->status->current = to;
}

void printKmerConnector(kmerConnector* kc, char* seq){
  printf("Uid: %.8x\nFlags: %u\n", kc->uid, kc->flags);
  printf("Dest: %s\nN: %d\nIds: ", seq, kc->n);
}

void printRead(memstruct* ms){
  kcLL* tmp = ms->status->trace;
  while (tmp){
    printf("^^^\n");
    kmerConnector* kc = tmp->kc;
    char* seq = (char*) calloc(ms->kmerSize + 1, sizeof(char));
    pos2seq(&ms, kc->dest, seq);
    printKmerConnector(kc, seq);
    printTIdList(kc->idflags);
    printf("\n");
    free(seq);
    tmp = tmp->next;
    printf("vvv\n");
  }
}

void printKmerConnectors(memstruct* ms, uint32_t pos){
  kmerConnector* kc = ms->kmerArray[pos];
  while (kc){
    char* seq = (char*) calloc(ms->kmerSize + 1, sizeof(char));
    pos2seq(&ms, kc->dest, seq);
    printKmerConnector(kc, seq);
    printTIdList(kc->idflags);
    printf("\n");
    free(seq);
    kc = kc->next;
  }
}

void printKc (memstruct* ms){
  kcLL* ll = ms->status->trace;
  char seq[12];
  char scurr[12];
  if (ll){
    pos2seq(&ms, ll->kc->dest, seq);
    pos2seq(&ms, ms->status->current, scurr);
    printf("-----KmerConnector------\n");
    printf("Current: %s\n", scurr);
    printf("Dest: %s\n", seq);
    if (!ll){
      printf("NULL\n");
      return;
    }
    while (ll){
      printKmerConnector(ll->kc, "");
      ll = ll->next;
    }
    printf("\n-----END KmerConnector------\n");
  }
}



void printKcLL(memstruct* ms, kcLL* k){
  while (k){
    char* seq = (char*) calloc(ms->kmerSize + 1, sizeof(char));
    pos2seq(&ms, k->kc->dest, seq);
    printKmerConnector(k->kc, seq);
    free(seq);
    k = k->next;
  }
}

void printMs(memstruct* ms){
  printf("---Memstruct---\n");

  printf("---End Memstruct---\n");
}

/*
 * kmerConnector functions //
 */



void summarize (memstruct* ms){
  uint32_t i = 0;
  for (i = 0; i < ms->nPos; i++){
    if (ms->kmerArray[i] && ms->kmerArray[i]->n > 0){
      char* seq = (char*) calloc(ms->kmerSize + 1, sizeof(char));
      pos2seq(&ms, i, seq);
      printf("----\n%s\n", seq);
      free(seq);
      seq = NULL;
      printKmerConnectors(ms, i);
    }
  }
}

uint32_t getStart(kmerHolder** khp, LISTTYPE nid){
  kmerHolder* kh = *khp;
  memstruct* ms = kh->ms;
  uint32_t result = NOKMER;
  for (uint32_t i = 0; i < ms->nPos; i++){
    kmerConnector* p = ms->kmerArray[i];
    if (p){
      tIdList* ok = _getTrace(&p->idflags, nid, 1);
      if (ok){
        result = i;
        return result;
      }
    }
  }
  return result;
}

uint32_t nTraces(kmerHolder** khp, bool printAnalysis){
  kmerHolder* kh = *khp;
  memstruct* ms = kh->ms;
  uint32_t MAX = 0xFFFF;
  uint32_t* nt = (uint32_t*) calloc(MAX, sizeof(uint32_t));
  uint32_t ntmax = 0;
  uint32_t result = 0;
  char* maxseq = (char*) calloc(kh->kmerSize + 1, sizeof(char));
  tIdList* cmplx = 0;
  for (uint32_t i = 0; i < ms->nPos; i++){
    kmerConnector* kc = ms->kmerArray[i];
    while (kc){
      if (kc->idflags){
        uint32_t nids = tIdListLength(&kc->idflags);
        if (nids > 0) nt[nids]++;
        if (nids > ntmax) ntmax = nids;
        if (nids == ntmax) cmplx = kc->idflags;
        pos2seq(&kh->ms, i, maxseq);
      }
      traceVessel* t = traceFirst(&kc->idflags);
      if (t){
        result++;
        t = t->next;
      }
      destroyTraceVessel(&t);
      kc = kc->next;
    }
  }
  if (printAnalysis){
    for (uint32_t i = 1; i < ntmax; i++){
      if (nt[i] > 0) printf("%" PRIu32 "\t%" PRIu32 "\n", i, nt[i]);
    }
  }
  free(nt);
  // Show most complex
  D_(0, "The winner is %s\n", maxseq);
  free(maxseq);
  tIdList* tmp = cmplx;
  while (tmp){
    LISTTYPE tid = tmp->trace.n;
    printf(">Tid %" PRIu32 "\n", tid);
    uint32_t start = getStart(khp, tid);
    char* sseq = (char*) calloc(kh->kmerSize + 1, sizeof(char));
    pos2seq(&kh->ms, start, sseq);
    D_(0, ">Start at %" PRIu32 ", seq %s\n", tid, sseq);
    free(sseq);
    D_(0, ">Tid %" PRIu32 "\n", tid);
    D_(0, "Start in %" PRIu32 "\n", start);
    if (start != NOKMER){
      kcLL* tseq = followTrace(khp, start, tid, 1);
      char* seq = getTraceSeq(khp, &tseq);
      printf("%s\n", seq);
    }
    else{
      printf("Error\n");
    }
    tmp = tmp->next;
  }
  return result;
}

void drawMs (memstruct* ms, char* fname){
  FILE* fp = fopen(fname, "w");
  if (fp){

  }
  else{
    fprintf(stderr, "Could not open %s\n", fname);
  }
  fclose(fp);
}

void resetKmer(kmerHolder** kp){
  kmerHolder* k = *kp;
  k->posInMemBuffer = 0;
  k->cVal = 0;
}

uint32_t updateKmer(kmerHolder** kp, char* b, void (*callback)(kmerHolder**, uint32_t)){
  kmerHolder* k = *kp;
  uint32_t result = NOKMER;
  uint8_t baseVal = base(&k->ms, b);
  D_(2, "Base %c, val %u\n", (int) b[0], baseVal);
  if (!isBase(&k->ms, b)){
    D_(2, "Non-standard base %d\n", (int) b[0]);
    if (b[0] > 0x30) resetTrace(kp);
    return NOKMER;
  }
  uint32_t ipos = k->posInMemBuffer - (uint32_t) k->kmerSize;
  if (k->posInMemBuffer >= (uint32_t) k->kmerSize){ //kmer filled
    uint32_t order = power(k->ms->bi->nbases, k->kmerSize - 1);
    unsigned char oldval = k->buffer[ipos];
    //printf("%d\t", ipos);
    k->cVal -=  order * (uint32_t) oldval;
  }
  k->cVal *= (uint32_t) k->ms->bi->nbases;
  k->cVal += (uint32_t) baseVal;
  k->buffer[k->posInMemBuffer] = baseVal;
  k->posInMemBuffer++;
  if (k->posInMemBuffer >= (KMERBUFFERSIZE - 1)){ // Rewind
    //printf("Rewinding\n");
    memmove(k->buffer, k->buffer + ipos + 1, sizeof(char) * (k->kmerSize));
    k->posInMemBuffer = (uint32_t) k->kmerSize;
  }
  if (k->posInMemBuffer >= (uint32_t) k->kmerSize){
    result = k->cVal;
  }
  //printf("%d-%c-%d\n", k->posInMemBuffer, b[0], result);

  if (result != NOKMER) callback(&k, result); //Side effects
  if (k->canonic) k->canonic = false;
  return result;
}

kmerHolder* initKmer(uint8_t kmerSize, uint8_t nBases){
  kmerHolder* result = (kmerHolder*) malloc(sizeof(kmerHolder));
  //
  if (nBases == 4){
    result->ms = initFourBase();
    result->nBases = 4;
  }
  else if (nBases == 5){
    result->ms = initFiveBase();
    result->nBases = 5;
  }
  else{
    fprintf(stderr, "At this point, only four (ACGT) or five (ACGTN) bases are accepted");
  }
  result->kmerSize = kmerSize;
  if (kmerSize >= KMERBUFFERSIZE){
    fprintf(stderr, "kmer size too large, max is %d", KMERBUFFERSIZE);
    exit(1);
  }
  result->buffer = (unsigned char*) calloc(KMERBUFFERSIZE, sizeof(unsigned char));
  setKmerArray(&result, kmerSize);
  result->posInMemBuffer = 0;
  result->cVal = 0;
  result->canonic = false;
  return result;
}

void mergeKh(kmerHolder** khp1, kmerHolder** khp2){
  //TODO: merge idflags
  kmerHolder* kh1 = *khp1;
  memstruct* ms1 = kh1->ms;
  kmerHolder* kh2 = *khp2;
  if (kh1->kmerSize != kh2->kmerSize){
    D_(0, "Incompatible kmerHolders\n");
    return;
  }
  for (uint32_t i = 0; i < kh2->ms->nPos; i++){
    kmerConnector* kc2 = kh2->ms->kmerArray[i];
    while (kc2){
      uint32_t from = i;
      uint32_t to = kc2->dest;
      kmerConnector* kc1 = getConnector(&ms1, from, to);
      kc1->flags = kc1->flags | kc2->flags;
      kc1->n += kc2->n;
      kc2 = kc2->next;
    }
  }
}

void filterMergeKh(kmerHolder** khp1, kmerHolder** khp2, uint16_t filter){
  //TODO: merge idflags
  kmerHolder* kh1 = *khp1;
  memstruct* ms1 = kh1->ms;
  kmerHolder* kh2 = *khp2;
  if (kh1->kmerSize != kh2->kmerSize){
    D_(0, "Incompatible kmerHolders\n");
    return;
  }
  for (uint32_t i = 0; i < kh2->ms->nPos; i++){
    kmerConnector* kc2 = kh2->ms->kmerArray[i];
    while (kc2 && kc2->n < filter){
      uint32_t from = i;
      uint32_t to = kc2->dest;
      kmerConnector* kc1 = getConnector(&ms1, from, to);
      kc1->flags = kc1->flags | kc2->flags;
      kc1->n += kc2->n;
      kc2 = kc2->next;
    }
  }
}

void setTraceDelme(memstruct** msp, kmerConnector** kcp, LISTTYPE i){
  kmerConnector* kcptr = *kcp;
  LISTTYPE cPos = 1;
  tIdList* p = _getTrace(&kcptr->idflags, i, cPos);
  while (kcptr){
    SET(p, DELME);
    kcptr = nextKc(msp, &kcptr, &p, cPos);
    cPos++;
  }
}

void _consolidate(kmerHolder** khp, kmerConnector** kcp, LISTTYPE i,
                  LISTTYPE n, LISTTYPE delta, kcLL** toDelp){
  // Delete trace n if it is compatible with i
  kmerHolder* kh = *khp;
  memstruct* ms = kh->ms;
  kmerConnector* kc = *kcp;
  LISTTYPE cPosI = delta;
  LISTTYPE cPosN = 1; // Position in trace n. Corresponds to pos delta in i
  kmerConnector* kcptr = kc;
  bool doConsol = true;
  bool doGoon   = true;
  tIdList* iptr = _getTrace(&kcptr->idflags, i, cPosI);
  tIdList* tmpiptr = iptr;
  if (i == n){
    printTIdList(tmpiptr); printf("\n\n");
    tIdList* nptr = _getTrace(&kcptr->idflags, n, cPosN);
    printTIdList(nptr); printf("\n");
    DIE("yo");
  }
  while (kcptr && doGoon){ // Check compatibility
    // This follows trace i and looks for trace n
    // The last time n is found, it should be LAST_IN_TRACE
    tIdList* nptr = _getTrace(&kcptr->idflags, n, cPosN);
    if (nptr){
      if (isPosLast(&nptr, cPosN)){
        tIdList* thisPosN = _getPos(&kcptr->idflags, n, cPosN);
        if (thisPosN && IS(thisPosN, LAST_IN_TRACE)){
          doGoon = false;
        }
      }
      // Could check here whether there is a nextKc in n
    }
    else{
      D_(2, "Cannot consolidate\n");
      doConsol = false;
      doGoon = false;
    }
    kcptr = nextKc(&ms, &kcptr, &tmpiptr, cPosI);
    cPosI++; cPosN++;
  }
  // Consolidate if compatible
  tIdList* consolTrace = _getTrace(&kc->idflags, n, 1);
  bool NalreadyDeleted = IS(consolTrace, DELME);
  if (doConsol){
    if (signal(SIGINT, sig_handler) == SIG_ERR)
    printf("\nCan't catch SIGINT\n");
    D_(1, "Trace %lu will be eliminated in favor of trace %lu\n", (LUI) n, (LUI) i);
    D_(1, "Starts at pos %lu\n", (LUI) delta);
    cPosI = delta;
    cPosN = 1;
    kcptr = kc;
    bool goon = true;
    kmerConnector* bridge = kcptr;
    tIdList* todel = NULL;
    while (kcptr && goon){ // Here they will be simply marked for deletion
      todel = _getTrace(&kcptr->idflags, n, cPosN);
      if (todel){
        if (!NalreadyDeleted) kcpush(toDelp, &kcptr, n, cPosN);
      }
      else{
        goon = false;
      }
      bridge = kcptr;
      kcptr = nextKc(&ms, &kcptr, &iptr, cPosI);
      cPosI++; cPosN++;
      D_(2, "CposI: %lu, CposN: %lu\n", (LUI) cPosI, (LUI) cPosN);
    }
    if (goon){ // Here, the remaining traces will be extended, and the old traces marked for deletion
      if (bridge) kcptr = nextKc(&ms, &bridge, &todel, cPosN - 1);
      tIdList* oldTrace = NULL;
      if (kcptr){
        // We go on, and therefore bridge should not be last in trace
        unsetAsLast(&bridge->idflags, i, cPosI - 1);
        D_(2, "Unsetting last in trace %lu\n", (LUI) i);
        oldTrace = _getTrace(&kcptr->idflags, n, cPosN);
      }
      while (kcptr){
        // The whole trace will be deleted afterwards
        if (!NalreadyDeleted) kcpush(toDelp, &kcptr, n, cPosN);
        tIdList* newTrace = addTrace(&kcptr->idflags, i);
        tIdList* oldPos = _getTrace(&oldTrace->posInTrace, cPosN, 0);
        tIdList* newPos = addPosInTrace(&newTrace, cPosI);
        newPos->trace.flag = oldPos->trace.flag;
        newPos->trace.nReads = oldPos->trace.nReads;
        newTrace->trace.flag = oldTrace->trace.flag;
        newTrace->trace.nReads = oldTrace->trace.nReads;
        SET(newTrace, IN_USE);
        kcptr = nextKc(&ms, &kcptr, &oldTrace, cPosN);
        cPosI++; cPosN++;
        D_(2, "CposI: %lu, CposN: %lu\n", (LUI) cPosI, (LUI) cPosN);
      }
    }
    setTraceDelme(&ms, &kc, n);
  }
}

void setTraceInUse(memstruct** msp, kmerConnector** kcp, tIdList** cTid){
  kmerConnector* kcptr = *kcp;
  LISTTYPE cPos = 1;
  tIdList* tid = *cTid;
  tIdList* p = tid;
  while (kcptr){
    SET(p, IN_USE);
    kcptr = nextKc(msp, &kcptr, &p, cPos);
    cPos++;
  }
}

void unsetTraceInUse(memstruct** msp, kmerConnector** kcp, tIdList** cTid){
  kmerConnector* kcptr = *kcp;
  tIdList* tid = *cTid;
  tIdList* p = tid;
  LISTTYPE cPos = 1;
  while (kcptr){
    UNSET(p, IN_USE);
    kcptr = nextKc(msp, &kcptr, &p, cPos);
    cPos++;
  }
}

void _canonizeThis(kmerHolder** khp, kmerConnector** kcp, LISTTYPE i,
                   kcLL** toDelp){
  kmerHolder* kh = *khp;
  memstruct* ms = kh->ms;
  kmerConnector* kc = *kcp;
  LISTTYPE cPos = 1;
  kmerConnector* kcptr = kc;
  // Set all in use
  tIdList* tid = _getTrace(&kc->idflags, i, cPos);
  tIdList* start = tid;
  if (!IS(tid, IN_USE)){
    setTraceInUse(&ms, kcp, &tid);
    D_(1, "Setting trace %" PRIu32 " in use\n", tid->trace.n);
    E_(2, printTIdList(start));
    while(kcptr){
      if (kcptr->idflags){
        tIdList* idptr = _getTrace(&kcptr->idflags, i, cPos);
        tIdList* p = _getTrace(&idptr->posInTrace, cPos, 0);
        if (!IS(p, DELME)){
          if (idptr){
            traceVessel* firsts = traceFirst(&kcptr->idflags);
            traceVessel* fstptr = firsts;
            E_(2, printTraceLL((traceLL*) firsts));
            while (fstptr && fstptr->tidl){
              traceVessel* nxt = fstptr->next;
              if (IS(fstptr->tidl, IN_USE)){
                D_(1, "Will not check trace %" PRIu32 ": In use\n", fstptr->tidl->trace.n);
              }
              else{
                D_(1, "Canonizing trace %" PRIu32 " with trace %" PRIu32 "\n", i, fstptr->tidl->trace.n);
                if (DEBUG >= 1){
                  khDebug = *khp;
                  kcDebug = kc;
                  iDebug  = fstptr->tidl->trace.n;
                  posDebug = 1;
                }
                if (!IS(fstptr->tidl, CANON) && !IS(fstptr->tidl, IN_USE)){
                  _canonizeThis(khp, &kcptr, fstptr->tidl->trace.n, toDelp);
                }
                if (!IS(fstptr->tidl, DELME) && !IS(fstptr->tidl, IN_USE)){
                  _consolidate(khp, &kcptr, i, fstptr->tidl->trace.n, cPos, toDelp);
                }
                SET(fstptr->tidl, CANON);
              }
              fstptr = nxt;
            }
            destroyTraceVessel(&firsts);
          }
        }
      }
      kcptr = nextKc(&ms, &kcptr, &tid, cPos);
      cPos++;
    }
    unsetTraceInUse(&ms, kcp, &start);
    D_(1, "Unetting trace %" PRIu32 " in use\n", start->trace.n);
  }
}

void _canonize(kmerHolder** khp){
  kmerHolder* kh = *khp;
  memstruct* ms  = kh->ms;
  uint32_t i = 0;//ms->status->current;
  // First pass
  D_(1, "Starting canonization\n");
  kcLL* toDel = NULL;
  for (i = 0; i < ms->nPos; i++){
    kmerConnector* kc = ms->kmerArray[i];
    while (kc && kc->n > 0){
      tIdList* l = kc->idflags;
      while (l){
        if (isFirst(&l) && !IS(l, CANON) && !IS(l, DELME)){
          D_(2, "Checking trace %lu\n", (LUI) l->trace.n)
          _canonizeThis(khp, &kc, l->trace.n, &toDel);
        }
        l = l->next;
      }
      kc = kc->next;
    }
  }
  //Unset IN_USE
  D_(1, "Unset canon flag\n");
  for (i = 0; i < ms->nPos; i++){
    kmerConnector* kc = ms->kmerArray[i];
    while (kc && kc->n > 0){
      tIdList* l = kc->idflags;
      while (l){
        if (isFirst(&l)){
          UNSET(l, CANON);
        }
        l = l->next;
      }
      kc = kc->next;
    }
  }
  D_(1, "Deleting merged traces\n");
  kcLL* tmp = toDel;
  while (tmp){
    if (tmp->kc){
      E_(1, printf("%" PRIu32 "\n", tmp->ntid); printTIdList(tmp->kc->idflags); printf("\n");)
      delTIdFromList(&tmp->kc->idflags, tmp->ntid, 0);
    }
    tmp = tmp->next;
  }
  resetKcLL(&toDel);
}



#endif // KMER_H_INCLUDED
