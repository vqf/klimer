#ifndef KMER_H_INCLUDED
#define KMER_H_INCLUDED

#ifndef DEBUG
#define DEBUG 1
#endif /* DEBUG */

#define ENCLOSE(a) printf("\n---before, line %d\n", __LINE__); a; printf("\nafter, line %d---\n", __LINE__);
#define D_(a, ...) if (DEBUG > a) printf(__VA_ARGS__);
#define P_(a) printf("%s at %d: %p\n", #a, __LINE__, a);

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "neatHtml.h"
#include "traceIdList.h"

#define KMERBUFFERSIZE 0xFFFF
#define NOKMER 0xFFFFFFFF
#define LUI long unsigned int

#ifndef BOOL
#define BOOL
typedef enum { false, true } bool;
#endif /* BOOL */

#ifndef IN_USE
#define IN_USE  0x01
#endif /* IN_USE */

// Set, unset and check flags in traces
// a -> tId object
// b -> flag to operate

#define SETKC(a, b) (((a)->flags = (a)->flags | b))
#define UNSETKC(a, b) ((a)->flags = (a)->flags ^ b)
#define ISKC(a, b) ((a)->flags & b)



// FLAGS

//


typedef struct kC{
  uint32_t dest;
  uint32_t uid;  // For debugging purposes, may lose later
  uint16_t n; // Number of events supporting connection
  uint8_t flags; // Connector-specific flags (in_use, ...)
  tIdList* idflags; // Info about each trace (id, isFirst, isLast, inUse...)
  struct kC* next; // Another kC, same from, different dest
} kmerConnector;

typedef struct mll{
  kmerConnector* kc;
  struct mll* next;
  struct mll* last;
} kcLL;

typedef struct bs{
  uint8_t* baseVal; //Fast conversion of bases to values (see initMemory)
  uint8_t* valBase; //Fast conversion from values to bases
  uint8_t nbases; // Usually 4, or 5 with N
} baseInfo;

typedef struct st{
  tIdList* extendMe;      // Which tIds can be extended?
  tIdList* traceSet; // Keep tab of trace IDs the current read belongs to
  LISTTYPE cId;      // Current id for kmerConnector
  kcLL* trace;      // Keeps track of the related kmers in the current trace
  uint32_t current; // Last kmer read
  bool start;       // No kc has still been filled
  bool isFirst;     // Is it the first in a trace?
  bool addExistingTrace;    // Should trace ids be added? (If not, cId is inserted)
  uint8_t addExistingTraceStatus; // See _existingTrace
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
  uint32_t cVal;
} kmerHolder;


void resetTrace(kmerHolder**);
void resetKmer(kmerHolder**);
void summarize (memstruct*);
void resetKcLL(kcLL**);
void printKmerConnector(kmerConnector*, char*);
void printKmerConnectors(memstruct*, uint32_t);
void printKc(memstruct*);
void printKcLL(memstruct*, kcLL*);
void destroyTIdList(tIdList**, void (*callback)(void**));
void delConnector(kmerConnector**);


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
  result->status = (msStatus*) calloc(1, sizeof(msStatus));
  result->status->trace = NULL;
  result->status->addExistingTrace = true;
  result->status->addExistingTraceStatus = 0;
  result->status->cId = 0;
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
  result->status->trace = NULL;
  result->status->addExistingTrace = true;
  result->status->addExistingTraceStatus = 0;
  result->status->extendMe = NULL;
  result->status->cId = 0;
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
  destroyTIdList(&ms->status->traceSet, destroyCircular);
  destroyTIdList(&ms->status->extendMe, destroyCircular);
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

bool isBase (memstruct* ms, char* b){
  uint8_t val = ms->bi->baseVal[(uint8_t) b[0]];
  bool result = false;
  if (val >= 0 && val < ms->bi->nbases) result = true;
  return result;
}

uint8_t base (memstruct* ms, char* b){
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
  resetTrace(kp);
  uint32_t nb = sizeof(kmerConnector) * ms->nPos;
  ms->nBytes = nb;
  if (DEBUG){
    fprintf(stderr, "Allocating %lu bytes...\n", (long unsigned int) nb);
  }
  ms->kmerArray = (kmerConnector**) calloc(ms->nPos, sizeof(kmerConnector*));
  if (DEBUG){
    fprintf(stderr, "Done\n");
  }
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
    if (isBase(ms, tb)){
      if (pos == NOKMER) pos = 0;
      pos = ms->bi->nbases * pos + base(ms, tb);
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
  if (pos < 0 || pos > ms->nPos){
    //fprintf(stderr, "Too much\n");
    return;
  }
  //fprintf(stderr, "Add to %lu\n", (long unsigned int) pos);
  if (ms->kmerArray[pos]->n <= 0xFFFF) ms->kmerArray[pos]->n++;
}

/*
 *With relationships
 */


void nextId(memstruct* ms){
  if (ms->status->cId == MAXLISTTYPE){
    if (DEBUG){
      fprintf(stderr, "%s\n", "Overflow in cId");
    }
    return;
  }
  ms->status->cId++;
}


kcLL* newKcPointer(kmerConnector** newkcp){
  kmerConnector* newkc = *newkcp;
  kcLL* result = (kcLL*) calloc(1, sizeof(kcLL));
  result->kc = newkc;
  result->next = NULL;
  result->last = result;
  return result;
}

void resetKcLL(kcLL** llp){
  kcLL* ll = *llp;
  while (ll){
    kcLL* todel = ll;
    ll = ll->next;
    free(todel);
  }
  llp = NULL;
}

void kcpush (kcLL** llp, kmerConnector** newkcp){
  kcLL* ll = *llp;
  if (ll){
    kcLL* tmp = ll->last;
    if (!ll->last) tmp = ll;
    tmp->next = newKcPointer(newkcp);
    ll->last = tmp->next;
  }
  else{
    ll = newKcPointer(newkcp);
    ll->last = ll;
  }
  *llp = ll;
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

bool isLastInTrace(memstruct* ms, kmerConnector* kc){
  bool result = true;
  kmerConnector* nxt = ms->kmerArray[kc->dest];
  if (nxt && nxt->dest){
    result = false;
  }
  return result;
}

bool isInUse(kmerConnector** t){
  bool result = ISKC(*t, IN_USE);
  return result;
}

void resetTrace(kmerHolder** kp){
  memstruct* ms = (*kp)->ms;
  kcLL** tmpp = &ms->status->trace;
  kcLL* tmp = *tmpp;
  bool markFirst = false;
  D_(1, "Resetting trace with status %d\n", ms->status->addExistingTraceStatus);
  if (ms->status->addExistingTraceStatus == 0 || ms->status->addExistingTraceStatus == 3){
    ms->status->addExistingTrace = false;
  }
  if (tmpp && tmp){
    // If there is no idflags, this is the first in a trace
    if (!tmp->kc->idflags){
      markFirst = true;
    }
    kcLL* kcCirc = NULL; // Which kcs in this trace loop
    traceLL* tCirc = newTraceLL(); // Which traces loop
    while(tmp){
      if (ms->status->addExistingTrace){
        mergeTIdLists(&tmp->kc->idflags, ms->status->traceSet, destroyCircular);
        unsetAsFirst(&tmp->kc->idflags, ms->status->traceSet);
        unsetAsLast(&tmp->kc->idflags, ms->status->traceSet);
      }
      else{
        //ENCLOSE(printTIdList(tmp->kc->idflags); printf("cID: %d\n", ms->status->cId));
        insertInTIdList(&tmp->kc->idflags, ms->status->cId, destroyCircular);
        insertInTIdList(&ms->status->traceSet, ms->status->cId, destroyCircular);
      }
      // Check for trace loops
      tIdList* cTraces = circTraces(&tmp->kc->idflags, destroyCircular);
      if (tmp->next && cTraces){
        pushTrace(&tCirc, cTraces, destroyCircular);
        traceVessel* tloop = _getTraces(&tmp->kc->idflags, cTraces);
        traceVessel* ptr = tloop;
        while (ptr->tidl){
          if (tmp->kc->idflags->trace.circular){
            kcpush((kcLL**) &ptr->tidl->trace.circular, &tmp->next->kc);
          }
          else{
            ptr->tidl->trace.circular = newKcPointer(&tmp->next->kc);
          }
          ptr = ptr->next;
        }
        destroyTraceVessel(&tloop);
      }
      destroyTIdList(&cTraces, destroyCircular);
      kcpush(&kcCirc, &tmp->kc);
      //ENCLOSE(printTIdList(tmp->kc->idflags); printTIdList(ms->status->traceSet));
      setAsUsed(&tmp->kc->idflags, ms->status->traceSet);
      //ENCLOSE(printTIdList(tmp->kc->idflags));
      tmp = tmp->next;
      nextTraceLL(&tCirc);
    }
    if (markFirst){
      tIdList* which = ms->status->traceSet;
      if (ms->status->addExistingTraceStatus == 0 || ms->status->addExistingTraceStatus == 3){
        which = ms->status->trace->kc->idflags;
      }
      setAsFirst(&ms->status->trace->kc->idflags, which);
    }
    //This is the last one in a read. Is it the last one in a trace?
    if (isLastInTrace(ms, ms->status->trace->last->kc)){
      tIdList* which = ms->status->traceSet;
      if (ms->status->addExistingTraceStatus == 0 || ms->status->addExistingTraceStatus == 3){
        which = ms->status->trace->last->kc->idflags;
      }
      setAsLast(&ms->status->trace->last->kc->idflags, which);
    }
    /* Set circ bits*/
    kcLL* kctmp    = kcCirc;
    traceLL* ttmp  = tCirc;
    while (kctmp){
      setCircular(&kctmp->kc->idflags, ttmp->tidl);
      D_(1, "Setting %d as Circ\n", ttmp->tidl->trace.n);
      kctmp = kctmp->next;
      ttmp  = ttmp->next;
    }
    /* Unset used bits */
    tmp = ms->status->trace;
    while (tmp){
      unsetInUse(&tmp->kc->idflags, ms->status->traceSet);
      tmp = tmp->next;
    }
    destroyTraceLL(&tCirc, destroyCircular);
    resetKcLL(&kcCirc);
  }
  /* Clean for next trace */
  ms->status->cId = 0;
  destroyTIdList(&ms->status->traceSet, destroyCircular);
  free(ms->status->traceSet);
  destroyTIdList(&ms->status->extendMe, destroyCircular);
  free(ms->status->extendMe);
  ms->status->traceSet = NULL;
  ms->status->addExistingTrace = true;
  ms->status->addExistingTraceStatus = 0;
  ms->status->start   = true;
  ms->status->isFirst = true;
  ms->status->current = NOKMER;
  resetKmer(kp);
  resetKcLL(&ms->status->trace);
  ms->status->trace = NULL;
}




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
  result->uid     = rnd32();
  result->n       = 0;
  result->idflags = NULL;
  result->next    = NULL;
  return result;
}

void delConnector(kmerConnector** kcp){
  kmerConnector* tmp = *kcp;
  while (tmp){
    kmerConnector* nxt = tmp->next;
    destroyTIdList(&tmp->idflags, destroyCircular);
    free(tmp);
    tmp = nxt;
  }
  *kcp = NULL;
}


/*
 * addExistingTrace is used to decide whether a new (cId) or existing (traceSet) trace is added.
 * addExistingTraceStatus can take several values that describe the current situation:
 *   0 -> No existing Ids have been added to traceSet. The situation is uncertain. If the trace
 * ends with this value, addExistingTrace should be false. If an existing id is found later, the
 * value should change to 1.
 *   1 -> Existing, non-empty Ids have been added to traceSet. If the trace ends with this value,
 * addExistingTrace should be true. If traceSet gets empty, its value should get to 2.
 * addExistingTrace can still be true if the last kc is the last in a trace (exists extendMe).
 * Otherwise, addExistingTrace should be false.
 *   2 -> There is a traceSet that is being extended. If the trace ends with this value,
 * addExistingTrace should be true. If a new traceId is found, its value should get to 3, and
 * addExistingTrace should be definitely false.
 *   3 -> A new trace (cId) should be added. In fact, this is equivalent to 0, but this value
 * is kept for debugging purposes
 *   4 -> A new trace (cId) should be added. This is final and happens when a read connects
 * two incompatible traces directly.
 */

void _existingTrace(memstruct** msp, kmerConnector** kcp){
  memstruct* ms = *msp;
  kmerConnector* kc = *kcp;
  tIdList* fst = traceFirst(kc->idflags, destroyCircular);
  if (kc->idflags){
    LISTTYPE m = maxInList(kc->idflags);
    if (ms->status->cId <= m){
      ms->status->cId = m + 1;
    }
  }
  if (ms->status->addExistingTraceStatus == 0 || ms->status->addExistingTraceStatus == 3){
    if (kc->idflags){
      if (fst){
        ms->status->traceSet = fst;
        ms->status->addExistingTraceStatus = 1;
        D_(1, "Found an existing trace: 0 to 1");
      }
      else if (ms->status->isFirst){
        ms->status->isFirst = false;
        ms->status->traceSet = copyTIdList(kc->idflags, destroyCircular);
        ms->status->addExistingTraceStatus = 1;
        D_(1, "Coincidence: 0 to 1");
      }
      else{
        ms->status->addExistingTraceStatus = 3;
        D_(1, "Found existing trace, not first: 0 to 3");
      }
    }
  }
  else if (ms->status->addExistingTrace){
    if (ms->status->addExistingTraceStatus == 1){
      intersectTIdLists(&ms->status->traceSet, kc->idflags, destroyCircular);
      if (!ms->status->traceSet){
        if (ms->status->extendMe){
          ms->status->addExistingTraceStatus = 2;
          destroyTIdList(&ms->status->traceSet, destroyCircular);
          ms->status->traceSet = ms->status->extendMe;
          D_(1, "Extending existing traces: 1 to 2");
        }
        else{
          ms->status->addExistingTraceStatus = 3;
          D_(1, "No traces to extend: 1 to 3");
        }
      }
    }
    else if (ms->status->addExistingTraceStatus == 2){
      if (kc->idflags){
        ms->status->addExistingTraceStatus = 3;
        D_(1, "Conflict: 2 to 3");
      }
    }
  }
  if (DEBUG) printf("TraceStatus: %d\n", ms->status->addExistingTraceStatus);
}


kmerConnector* getConnector(memstruct** msp, uint32_t from, uint32_t to){
  memstruct *ms = *msp;
  kmerConnector* result = ms->kmerArray[from];
  if (!result){ // This connector did not exist
    ms->kmerArray[from] = newKmerConnector(to);
    result = ms->kmerArray[from];
    _existingTrace(msp, &result);
  }
  else{
    while (result->dest != to && result->next){
      result = result->next;
    }
    if (result->dest == to){
      _existingTrace(msp, &result);
    }
    else{ // This connector did not exist
      kmerConnector* nkc = newKmerConnector(to);
      result->next = nkc;
      result = result->next;
      _existingTrace(msp, &result);
    }
  }
  //ms->status->cId = maxInList(ms->status->traceSet);
  ms->status->extendMe = traceLast(result->idflags, destroyCircular);
  kcpush(&ms->status->trace, &result);
  return result;
}

void addRelationship(kmerHolder** kp, uint32_t to){
  memstruct* ms = (*kp)->ms;
  if (ms->status->start){
    ms->status->start = false;
    ms->status->current = to;
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
  printf("Uid: %.8x\n", kc->uid);
  printf("Dest: %s\nN: %d\nIds: ", seq, kc->n);
}

void printKmerConnectors(memstruct* ms, uint32_t pos){
  kmerConnector* kc = ms->kmerArray[pos];
  while (kc){
    char* seq = (char*) calloc(ms->kmerSize + 1, sizeof(char));
    pos2seq(&ms, kc->dest, seq);
    printKmerConnector(kc, seq);
    printTIdList(kc->idflags);
    tIdList* tlist = kc->idflags;
    while (tlist){
      kcLL* tmp = (kcLL*) tlist->trace.circular;
      if (tmp){
        printf("Circ%u: ", tlist->trace.n);
        while (tmp){
          printf("%.8x\t", tmp->kc->uid);
          tmp = tmp->next;
        }
        printf("\n");
      }
      tlist = tlist->next;
    }
    free(seq);
    kc = kc->next;
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

void writeOut(memstruct* ms, char* fname){
  FILE* fout = fopen(fname, "w");
  fwrite(ms->kmerArray, ms->nBytes, 1, fout);
}

/*void readIn(memstruct* ms, char* fname){
  FILE* file = fopen(fname, "r");
  fseek(file, 0, SEEK_END);
  off_t fileLen = ftell(file);
  fseek(file, 0, SEEK_SET);
  if (fileLen < ms->nPos){
    fprintf(stderr, "File not suitable for %d-mers\n", (int) ms->kmerSize);
    exit(1);
  }
  fread(ms->kmerArray, ms->nPos, 1, file);
  fclose(file);
}*/

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
  uint8_t baseVal = base(k->ms, b);
  if (DEBUG){
    printf("Base %c, val %u\n", (int) b[0], baseVal);
  }
  if (!isBase(k->ms, b)){
    if (DEBUG){
      fprintf(stderr, "Non-standard base %d\n", (int) b[0]);
    }
    resetTrace(kp);
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
  return result;
}

kmerHolder* initKmer(uint8_t kmerSize){
  kmerHolder* result = (kmerHolder*) malloc(sizeof(kmerHolder));
  result->kmerSize = kmerSize;
  if (kmerSize >= KMERBUFFERSIZE){
    fprintf(stderr, "kmer size too large, max is %d", KMERBUFFERSIZE);
    exit(1);
  }
  result->buffer = (unsigned char*) calloc(KMERBUFFERSIZE, sizeof(unsigned char));
  result->ms = initFourBase();
  setKmerArray(&result, kmerSize);
  result->posInMemBuffer = 0;
  result->cVal = 0;
  return result;
}



#endif // KMER_H_INCLUDED
