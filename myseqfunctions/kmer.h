#ifndef KMER_H_INCLUDED
#define KMER_H_INCLUDED

#ifndef DEBUG
#define DEBUG 1
#endif /* DEBUG */

#define ENCLOSE(a) printf("\n---before, line %d\n", __LINE__); a; printf("\nafter, line %d---\n", __LINE__);
#define D_(a, ...) if (DEBUG >= a) printf(__VA_ARGS__);
#define E_(a, b) if (DEBUG >= a){ printf("At %d: ", __LINE__); b; }
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
int _CTR = 1;
uint32_t _UID = 0;
#define SCALAR(a, b) _CTR = 1; \
                  while (a->next){\
                    _CTR++; \
                    a = a->next;\
                  } \
                  b = _CTR;

#define SETKC(a, b) (((a)->flags = (a)->flags | b))
#define UNSETKC(a, b) ((a)->flags = (a)->flags ^ b)
#define ISKC(a, b) ((a)->flags & b)
#define GETUID _UID; _UID++;



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
  uint8_t flags;
  uint32_t pos;
  struct mll* next;
  struct mll* last;
} kcLL;

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
  bool extendMeUp;      // Which tIds can be extended?
  bool extendMeDn;      // Which tIds can be extended?
  tIdList* extendMe;
  tIdList* traceSet; // Keep tab of trace IDs the current read belongs to
  LISTTYPE cId;      // Current id for kmerConnector
  kcLL* trace;      // Keeps track of the related kmers in the current trace
  uint32_t current; // Last kmer read
  bool start;       // No kc has still been filled
  bool isFirst;     // Is it the first in a trace?
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
  uint8_t nBases;
  uint32_t cVal;
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
void destroyTIdList(tIdList**, void (*callback)(void**));
void delConnector(kmerConnector**);

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
  result->extendMe = NULL;
  result->extendMeUp = false;
  result->extendMeDn = false;
  result->cId = 0;
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
  result->flags = 0x00;
  result->pos   = 0x00000000;
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

void extendCircUp(tmpKcVessel** kvp){
  tmpKcVessel* tmp = *kvp;
  while (tmp){
    kcLL* toUnshift = NULL;
    kcpush(&toUnshift, &tmp->loop);
    toUnshift->next = tmp->dest->circular;
    tmp->dest->circular = toUnshift;
    tmp = tmp->next;
  }
}

void checkForLoops(kcLL** tmpp, tIdList* which){
  kcLL* tmp = *tmpp;
  while (tmp){
    traceVessel* t = _getTraces(&tmp->kc->idflags, which);
    traceVessel* pointer = t;
    while (pointer->tidl){
      if (IS(pointer->tidl, RESERVED) && !IS(pointer->tidl, IN_USE)){
        SET(pointer->tidl, IN_USE);
        if (pointer->tidl->trace.circular && !IS(pointer->tidl, CIRCULAR)){ // Comes from oneMore
          resetKcLL((kcLL**) &pointer->tidl->trace.circular);
          pointer->tidl->trace.circular = NULL;
          D_(2, "Used to be oneMore. Not anyMore\n");
        }
      }
      else{
        SET(pointer->tidl, RESERVED);
      }
      pointer->tidl->trace.nReads++;
      D_(2, "UID: %.8x, trace: %lu\n", (unsigned int) tmp->kc->uid, (LUI) pointer->tidl->trace.n);
      pointer = pointer->next;
    }
    destroyTraceVessel(&t);
    tmp = tmp->next;
  }
  tmp = *tmpp;
  while (tmp){
    traceVessel* t = _getTraces(&tmp->kc->idflags, which);
    traceVessel* pointer = t;
    while (pointer->tidl){
      if (IS(pointer->tidl, RESERVED)) UNSET(pointer->tidl, RESERVED);
      pointer = pointer->next;
    }
    destroyTraceVessel(&t);
    tmp = tmp->next;
  }
}

void findTraceSet(memstruct** msp){
  memstruct* ms = *msp;
  kcLL** tmpp = &ms->status->trace;
  kcLL* tmp = *tmpp;
  bool newTrace = (ms->status->addExistingTraceStatus == 0 || ms->status->addExistingTraceStatus == 3);
  while (tmp){
    if (newTrace){
      if (ms->status->traceSet){
        D_(0, "Warning: traceSet not empty at newTrace\n");
        //printTIdList(ms->status->traceSet);
        destroyTIdList(&ms->status->traceSet, destroyCircular);
      }
      insertInTIdList(&tmp->kc->idflags, ms->status->cId, destroyCircular);
      insertInTIdList(&ms->status->traceSet, ms->status->cId, destroyCircular);
      D_(2, "Adding existing trace %lu\n", (LUI) ms->status->cId);
    }
    else{

      mergeTIdLists(&tmp->kc->idflags, ms->status->traceSet, destroyCircular);
      D_(2, "Inserting traceSet in idflags\n");
      E_(2, printTIdList(ms->status->traceSet));
    }
    tmp = tmp->next;
  }
}

void postProcess(memstruct** msp){
  memstruct* ms = *msp;
  kcLL** tmpp = &ms->status->trace;
  kcLL* tmp = *tmpp;
  tIdList* lTraces = ms->status->traceSet;
  while (lTraces){
    tIdList* t = NULL;
    if (tmp) t = _getTrace(&tmp->kc->idflags, lTraces->trace.n);
    while (tmp && tmp->next){
      bool oneMore = false;
      if (ISKC(tmp->kc, RESERVED)){
        oneMore = true;
        if (tmp) t = _getTrace(&tmp->kc->idflags, lTraces->trace.n);
      }
      while (tmp->next && t && IS(t, CIRCULAR)){
        if (!oneMore) oneMore = true;
        tmp = tmp->next;
        if (tmp) t = _getTrace(&tmp->kc->idflags, lTraces->trace.n);
      }
      if (tmp && tmp->next){
        if (t && oneMore && !t->trace.circular){
          D_(2, "Setting oneMore between %.8x and %.8x\n", (unsigned int) tmp->kc->uid, (unsigned int) tmp->next->kc->uid);
          kcpush((kcLL**) &t->trace.circular, &tmp->next->kc);
        }
        tmp = tmp->next;
        if (tmp) t = _getTrace(&tmp->kc->idflags, lTraces->trace.n);
      }
      tmp = tmp->next;
      if (tmp) t = _getTrace(&tmp->kc->idflags, lTraces->trace.n);
    }
    lTraces = lTraces->next;
  }
}

void cleanTraceStatus(kmerHolder** khp){
  memstruct* ms = (*khp)->ms;
  ms->status->cId = 0;
  destroyTIdList(&ms->status->traceSet, destroyCircular);
  free(ms->status->traceSet);
  destroyTIdList(&ms->status->extendMe, destroyCircular);
  free(ms->status->extendMe);
  ms->status->traceSet = NULL;
  ms->status->addExistingTraceStatus = 0;
  ms->status->start   = true;
  ms->status->isFirst = true;
  ms->status->current = NOKMER;
  resetKmer(khp);
  resetKcLL(&ms->status->trace);
  ms->status->trace = NULL;
}

void resetTrace(kmerHolder** kp){
  memstruct* ms = (*kp)->ms;
  kcLL** tmpp = &ms->status->trace;
  kcLL* tmp = *tmpp;
  bool firstInRead = true;
  bool extendingUp = false;
  bool extendingDn = false;
  bool hasLoops    = false;
  bool canCheckLoop = true;
  bool newTrace = (ms->status->addExistingTraceStatus == 0 || ms->status->addExistingTraceStatus == 3);
  D_(1, "Resetting trace with status %d\n", ms->status->addExistingTraceStatus);
  if (newTrace){
    canCheckLoop = true;
    D_(1, "This is a new trace\n");
  }
  // Are we extending from upstream?
  if (ms->status->extendMeUp && !newTrace){
    extendingUp = true; // Until we get to first_in_trace
    canCheckLoop = true;
    D_(1, "Extending up\n");
  }
  if (tmpp && tmp){
    //Find out and add traceSet
    findTraceSet(&ms);
    E_(1, printTIdList(ms->status->traceSet));
    //Pre-check for loops
    checkForLoops(&ms->status->trace, ms->status->traceSet);
    if (!tmp->next){
      D_(2, "Got to the end of the read\n");
    }
    //printRead(ms);
    kcLL* kcCirc = NULL; // Which kcs in this trace loop
    traceLL* tCirc = newTraceLL(); // Which traces loop
    tmpKcVessel* upXt = NULL;
    while(tmp){
      D_(1, "Adding kc with uid: %.8x\n", tmp->kc->uid);
      if (extendingUp && isTraceFirst(&tmp->kc->idflags, ms->status->traceSet, destroyCircular)){
        extendingUp = false; // Now download the contents of upXt
        extendCircUp(&upXt);
        canCheckLoop = false;
        D_(1, "We are getting into known territory: %.8x\n", tmp->kc->uid);
      }
      if (ms->status->extendMeDn && !extendingDn && isTraceLast(&tmp->kc->idflags, ms->status->traceSet, destroyCircular)){
        extendingDn = true;
        canCheckLoop = true;
        D_(1, "Extending down\n");
      }
      // Check for trace loops
      if (canCheckLoop){
        bool overrideCirc = (extendingUp | extendingDn);
        tIdList* cTraces = circTraces(&tmp->kc->idflags, overrideCirc, destroyCircular);
        if (tmp->next && cTraces){
          hasLoops = true;
          D_(2, "This kc loops with override %u\n", overrideCirc);
          pushTrace(&tCirc, cTraces, destroyCircular);
          kcpush(&kcCirc, &tmp->kc);
          nextTraceLL(&tCirc);
          traceVessel* tloop = _getTraces(&tmp->kc->idflags, cTraces);
          traceVessel* ptr = tloop;
          while (ptr->tidl){
            if (extendingUp){
              tmpKcVessel* prev = newTmpKcVessel(&ptr->tidl->trace, &tmp->next->kc, &upXt);
              upXt = prev;
              D_(2, "Unshifting temp loopers\n");
            }
            else{
              kcpush((kcLL**) &ptr->tidl->trace.circular, &tmp->next->kc);
              D_(2, "Adding to circular\n");
              E_(2, printKcLL(ms, ptr->tidl->trace.circular));
            }
            ptr = ptr->next;
          }
          destroyTraceVessel(&tloop);
        }
        destroyTIdList(&cTraces, destroyCircular);
      }
      else{
        D_(2, "Cannot check for loop\n");
      }
      setAsUsed(&tmp->kc->idflags, ms->status->traceSet);
      if (firstInRead){
        firstInRead = false;
      }
      else{
        unsetAsFirst(&tmp->kc->idflags, ms->status->traceSet);
      }
      if (tmp->next) unsetAsLast(&tmp->kc->idflags, ms->status->traceSet);
      tmp = tmp->next;
    }
    destroytmpKcVessel(&upXt);
    if (newTrace || ms->status->extendMeUp){
      setAsFirst(&ms->status->trace->kc->idflags, ms->status->traceSet);
      D_(1, "Marking as first: %.8x\n", ms->status->trace->kc->uid);
    }
    if (newTrace || ms->status->extendMeDn){
      setAsLast(&ms->status->trace->last->kc->idflags, ms->status->traceSet);
      D_(1, "Marking as last: %.8x\n", ms->status->trace->last->kc->uid);
    }
    /* Set circ bits*/
    kcLL* kctmp    = kcCirc;
    traceLL* ttmp  = tCirc;
    while (kctmp){
      setCircular(&kctmp->kc->idflags, ttmp->tidl);
      D_(1, "Setting %d in %.8x as Circ\n", ttmp->tidl->trace.n, kctmp->kc->uid);
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
  postProcess(&ms);
  /* Clean for next trace */
  cleanTraceStatus(kp);
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
  //D_(2, "yo\t"); printTIdList(fst);
  if (kc->idflags){
    LISTTYPE m = maxInList(kc->idflags);
    if (ms->status->cId <= m){
      ms->status->cId = m + 1;
    }
  }
  if (ms->status->addExistingTraceStatus == 0 || ms->status->addExistingTraceStatus == 3){
    if (kc->idflags){
      destroyTIdList(&ms->status->traceSet, destroyCircular);
      tIdList* fst = traceFirst(kc->idflags, destroyCircular);
      if (fst){
        ms->status->traceSet = fst;
        fst = NULL;
        ms->status->extendMeUp = true;
        ms->status->addExistingTraceStatus = 1;
        D_(2, "Found an existing trace: 0 to 1\n");
      }
      else if (ms->status->isFirst){
        ms->status->isFirst = false;
        ms->status->extendMeUp = false;
        ms->status->traceSet = copyTIdList(kc->idflags, destroyCircular);
        ms->status->addExistingTraceStatus = 1;
        D_(2, "Inside, not extending: 0 to 1\n");
      }
      else{
        ms->status->addExistingTraceStatus = 3;
        ms->status->extendMeUp = false;
        ms->status->extendMeDn = false;
        D_(2, "Found existing trace, not first: to 3\n");
      }
      destroyTIdList(&fst, destroyCircular);
    }
  }
  else if (ms->status->addExistingTraceStatus == 1){
    intersectTIdLists(&ms->status->traceSet, kc->idflags, destroyCircular);
    if (!ms->status->traceSet){
      if (ms->status->extendMe){
        ms->status->addExistingTraceStatus = 2;
        destroyTIdList(&ms->status->traceSet, destroyCircular);
        ms->status->traceSet = copyTIdList(ms->status->extendMe, destroyCircular);
        ms->status->extendMeDn = true;
        D_(2, "Extending existing traces: 1 to 2\n");
      }
      else{
        ms->status->addExistingTraceStatus = 3;
        ms->status->extendMeUp = false;
        ms->status->extendMeDn = false;
        destroyTIdList(&ms->status->traceSet, destroyCircular);
        D_(2, "No traces to extend: 1 to 3\n");
      }
    }
  }
  else if (ms->status->addExistingTraceStatus == 2){
    if (kc->idflags){
      /*if (isIncluded(&ms->status->traceSet, &kc->idflags)){
        destroyTIdList(&ms->status->traceSet, destroyCircular);
        ms->status->traceSet = copyTIdList(kc->idflags, destroyCircular);
        D_(2, "No conflict, ids compatible\n");
      }
      else{
        ms->status->addExistingTraceStatus = 3;
        ms->status->extendMeUp = false;
        ms->status->extendMeDn = false;
        D_(2, "Conflict: 2 to 3\n");
      }*/
      ms->status->addExistingTraceStatus = 3;
      ms->status->extendMeUp = false;
      ms->status->extendMeDn = false;
      destroyTIdList(&ms->status->traceSet, destroyCircular);
      D_(2, "Conflict: 2 to 3\n");
    }
  }
  D_(2, "TraceStatus: %d\n", ms->status->addExistingTraceStatus);
}


kmerConnector* getConnector(memstruct** msp, uint32_t from, uint32_t to){
  memstruct *ms = *msp;
  kmerConnector* result = ms->kmerArray[from];
  D_(2, "Getting connector from %lu to %lu\n", (LUI) from, (LUI) to);
  if (!result){ // This connector did not exist
    D_(2, "This is a new one\n");
    ms->kmerArray[from] = newKmerConnector(to);
    result = ms->kmerArray[from];
    _existingTrace(msp, &result);
  }
  else{
    while (result->dest != to && result->next){
      result = result->next;
    }
    if (result->dest == to){
      D_(2, "Oh, you mean this one\n");
      _existingTrace(msp, &result);
    }
    else{ // This connector did not exist
      kmerConnector* nkc = newKmerConnector(to);
      result->next = nkc;
      result = result->next;
      D_(2, "New destination\n");
      _existingTrace(msp, &result);
    }
  }
  //ms->status->cId = maxInList(ms->status->traceSet);
  destroyTIdList(&ms->status->extendMe, destroyCircular);
  ms->status->extendMe = traceLast(result->idflags, destroyCircular);
  if (from == to){
    SETKC(result, RESERVED);
  }
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
    tIdList* tlist = kc->idflags;
    while (tlist){
      kcLL* tmpll = (kcLL*) tlist->trace.circular;
      if (tmpll){
        printf("Circ%u: ", tlist->trace.n);
        while (tmpll){
          printf("%.8x\t", tmpll->kc->uid);
          tmpll = tmpll->next;
        }
        printf("\n");
      }
      tlist = tlist->next;
    }
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
  D_(2, "Base %c, val %u\n", (int) b[0], baseVal);
  if (!isBase(k->ms, b)){
    D_(2, "Non-standard base %d\n", (int) b[0]);
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

kmerHolder* initKmer(uint8_t kmerSize, uint8_t nBases){
  kmerHolder* result = (kmerHolder*) malloc(sizeof(kmerHolder));
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
  return result;
}



#endif // KMER_H_INCLUDED
