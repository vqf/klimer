#ifndef KMER_H_INCLUDED
#define KMER_H_INCLUDED

#ifndef DEBUG
#define DEBUG 0
#endif /* DEBUG */

#define ENCLOSE(a) printf("\n---before, line %d\n", __LINE__); a; printf("\nafter, line %d---\n", __LINE__);
#define SEQ(a, b) char* _seq = (char*) calloc(12, sizeof(char)); pos2seq(a, b->dest, _seq); printf("%s\n", _seq); free(_seq);

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
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
  LISTTYPE ntid;
  LISTTYPE posInTrace;
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

typedef struct tSet{
  tIdList* starting; // Where each trace starts
  tIdList* current;  // Where each trace is
} tSet;

typedef struct st{
  bool extendMeUp;      // Which tIds can be extended?
  bool extendMeDn;      // Which tIds can be extended?
  tIdList* extendMe;
  tSet* traceSet; // Keep tab of trace IDs the current read belongs to
  LISTTYPE cId;      // Current id for kmerConnector
  kcLL* trace;      // Keeps track of the related kmers in the current trace
  uint32_t firsInTrace;  // Can keep tab of current when following trace
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
void destroyTIdList(tIdList**);
kmerConnector* getConnector(memstruct**, uint32_t, uint32_t);
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
  result->traceSet = (tSet*) calloc(1, sizeof(tSet));
    result->traceSet->starting = NULL;
    result->traceSet->current = NULL;
  result->firsInTrace = NOKMER;
  result->current = NOKMER;
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
  destroyTIdList(&ms->status->traceSet->starting);
  destroyTIdList(&ms->status->traceSet->current);
  free(ms->status->traceSet);
  destroyTIdList(&ms->status->extendMe);
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
  cleanTraceStatus(kp);
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


void nextId(memstruct** msp){
  memstruct* ms = *msp;
  if (ms->status->cId == MAXLISTTYPE){
    if (DEBUG){
      fprintf(stderr, "%s\n", "Overflow in cId");
    }
    return;
  }
  ms->status->cId++;
}


kcLL* newKcPointer(kmerConnector** newkcp, LISTTYPE ntid, LISTTYPE pos){
  kmerConnector* newkc = *newkcp;
  kcLL* result = (kcLL*) calloc(1, sizeof(kcLL));
  result->kc = newkc;
  result->next = NULL;
  result->ntid = ntid;
  result->posInTrace = pos;
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

void kcpush (kcLL** llp, kmerConnector** newkcp, LISTTYPE ntid, LISTTYPE pos){
  kcLL* ll = *llp;
  if (ll){
    kcLL* tmp = ll->last;
    if (!ll->last) tmp = ll;
    tmp->next = newKcPointer(newkcp, ntid, pos);
    ll->last = tmp->next;
  }
  else{
    ll = newKcPointer(newkcp, ntid, pos);
    ll->last = ll;
  }
  *llp = ll;
}

kcLL* kcCopy(kcLL** kclp){
  kcLL* kcl = *kclp;
  kcLL* result = NULL;
  if (kcl){
    //uint32_t pos = kcl->posInTrace;
    while (kcl){
      kcpush(&result, &kcl->kc, kcl->ntid, kcl->posInTrace);
      kcl = kcl->next;
    }
    //result->posInTrace = pos;
  }
  return result;
}

void kcConcat(kcLL** llp, kcLL** llp2){
  kcLL* ll1 = *llp;
  kcLL* ll2 = *llp2;
  if (ll1){
    ll1 = ll1->last;
    ll1->next = ll2;
  }
  else{
    ll1 = ll2;
  }
  llp2 = NULL;
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

kmerConnector* getKcWithTId(memstruct** msp, uint32_t pos, LISTTYPE tid, LISTTYPE posInTrace){
  memstruct* ms = *msp;
  kmerConnector* kc = ms->kmerArray[pos];
  kmerConnector* result = NULL;
  while (kc){
    tIdList* tt = _getTrace(&kc->idflags, tid, posInTrace);
    if (tt){
      if (result){
        if (!IS(tt, FIRST_IN_TRACE) && !IS(tt, LAST_IN_TRACE)){
          result = kc;
        }
      }
      else{
        result = kc;
      }
    }
    kc = kc->next;
  }
  /*if (!result && !includesLast){
    char* seq = (char*) calloc(12, sizeof(char));
    pos2seq(&ms, pos, seq);
    D_(0, "Error: Could not find trace %lu at %lu (%s)\n", (LUI) tid, (LUI) pos, seq);
    free(seq);
    //exit(0);
  }*/
  return result;
}
/*
bool isNextInConflict(memstruct** msp, uint32_t pos, LISTTYPE tid){
  bool result = false;
  memstruct* ms = *msp;
  kmerConnector* kc = ms->kmerArray[pos];
  kmerConnector* nxt = NULL;
  while (kc){
    tIdList* tt = _getTrace(&kc->idflags, tid);
    if (tt){
      if (nxt){
        if (!IS(tt, FIRST_IN_TRACE) && !IS(tt, LAST_IN_TRACE)){
          result = true;
        }
      }
      else{
        nxt = kc;
      }
    }
    kc = kc->next;
  }
  return result;
}
*/
kmerConnector* nextKc(memstruct** msp, kmerConnector** kcp, LISTTYPE i, LISTTYPE pos){
  kmerConnector* kc = *kcp;
  D_(2, "From %.8x\n", (unsigned int) kc->uid);
  kmerConnector* result = NULL;
  result = getKcWithTId(msp, kc->dest, i, pos);
  return result;
}


bool isLastInTrace(memstruct* ms, kmerConnector* kc, LISTTYPE i, LISTTYPE pos){
  tIdList* myl = _getTrace(&kc->idflags, i, pos);
  if (IS(myl, LAST_IN_TRACE)){
    return true;
  }
  return false;
}



void findTraceSet(memstruct** msp){
  memstruct* ms = *msp;
  kcLL** tmpp = &ms->status->trace;
  kcLL* tmp = *tmpp;
  ms->status->current = ms->status->firsInTrace;
  bool newTrace = (ms->status->addExistingTraceStatus == 0 || ms->status->addExistingTraceStatus == 3);
  if (newTrace && ms->status->traceSet){
    D_(0, "Warning: traceSet not empty at newTrace\n");
    //printTIdList(ms->status->traceSet);
    destroyTIdList(&ms->status->traceSet->starting);
    destroyTIdList(&ms->status->traceSet->current);
  }
  while (tmp){
    if (newTrace){
      insertInTIdList(&tmp->kc->idflags, ms->status->cId);
      insertInTIdList(&ms->status->traceSet->starting, ms->status->cId);
      D_(2, "Adding existing trace %lu\n", (LUI) ms->status->cId);
    }
    else{
      mergeTIdLists(&tmp->kc->idflags, ms->status->traceSet->starting);
      D_(2, "Inserting traceSet in idflags\n");
      E_(2, printTIdList(ms->status->traceSet->starting));
    }
    ms->status->current = tmp->kc->dest;
    tmp = tmp->next;
  }
}


void cleanTraceStatus(kmerHolder** khp){
  memstruct* ms = (*khp)->ms;
  //ms->status->cId = 0;
  destroyTIdList(&ms->status->traceSet->starting);
  destroyTIdList(&ms->status->traceSet->current);
  destroyTIdList(&ms->status->extendMe);
  free(ms->status->extendMe);
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
  bool newTrace = (ms->status->addExistingTraceStatus == 0 || ms->status->addExistingTraceStatus == 3);
  D_(1, "Resetting trace with status %d\n", ms->status->addExistingTraceStatus);
  tIdList* tracesToAdd = NULL;
  if (newTrace){
    nextId(&ms);
    D_(1, "This is a new trace\n");
    if (ms->status->traceSet){
       D_(0, "Warning: traceSet not empty at newTrace\n");
      destroyTIdList(&ms->status->traceSet->starting);
      destroyTIdList(&ms->status->traceSet->current);
    }
    LISTTYPE cposInTrace = 1;
    while (tmp){
      tIdList* tl = addTrace(&tmp->kc->idflags, ms->status->cId);
      addPosInTrace(&tl, cposInTrace);
      cposInTrace++;
      tmp = tmp->next;
    }
  }
  else{
    //Find out and add traceSet (passed to traceToAdd).
    tracesToAdd = copyTIdList(&ms->status->traceSet->starting);
    if (ms->status->extendMeUp){
      extendingUp = true; // Until we get to first_in_trace
      D_(1, "Extending Up\n");
    }
  }
  if (tmpp && tmp){
    E_(1, printTIdList(ms->status->traceSet->starting));
    while(tmp){
      D_(1, "Adding kc with uid: %.8x\n", tmp->kc->uid);
      if (extendingUp && isTraceFirst(&tmp->kc->idflags, ms->status->traceSet->starting)){
        extendingUp = false; // Now refresh posInTrace in the rest of the trace

        D_(1, "We are getting into known territory: %.8x\n", tmp->kc->uid);
      }
      if (ms->status->extendMeDn && !extendingDn && isTraceLast(&tmp->kc->idflags, ms->status->traceSet->starting)){
        extendingDn = true;
        D_(1, "Extending down\n");
      }
      if (firstInRead){
        firstInRead = false;
      }
      else{
        unsetAsFirst(&tmp->kc->idflags, ms->status->traceSet->current);
      }
      if (tmp->next){
        unsetAsLast(&tmp->kc->idflags, ms->status->traceSet->current);
      }
      tmp = tmp->next;
    }
    if (newTrace || ms->status->extendMeUp){
      setAsFirst(&ms->status->trace->kc->idflags, ms->status->traceSet->starting);
      D_(1, "Marking as first: %.8x\n", ms->status->trace->kc->uid);
    }
    if (newTrace || ms->status->extendMeDn){
      setAsLast(&ms->status->trace->last->kc->idflags, ms->status->traceSet->starting);
      D_(1, "Marking as last: %.8x\n", ms->status->trace->last->kc->uid);
    }
  }
  destroyTIdList(&tracesToAdd);
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
    destroyTIdList(&tmp->idflags);
    free(tmp);
    tmp = nxt;
  }
  *kcp = NULL;
}

void trimUncompatibleTraces(memstruct** msp, tSet** tp, kmerConnector** nxtp){
  memstruct* ms = *msp;
  tSet* t    = *tp;
  kmerConnector* nxt = *nxtp;
  intersectTIdLists(&t->starting, nxt->idflags);
  intersectTIdLists(&t->current, nxt->idflags);
  tIdList* current = t->current;
  while (current){
    tIdList* tPos = current->posInTrace;
    while (tPos){
      LISTTYPE nnext = tPos->trace.n + 1;
      tIdList* todel = tPos;
      tPos = tPos->next;
      if (!isInTIdList(&nxt->idflags, nnext)){
        delTIdFromList(&t->current->posInTrace, todel->trace.n);
      }
    }
    tIdList* todel = current;
    current = current->next;
    if (!todel->posInTrace){
      delTIdFromList(&t->current, todel->trace.n);
    }
  }
  intersectTIdLists(&t->starting, t->current);
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
 */

void _existingTrace(memstruct** msp, kmerConnector** kcp){
  memstruct* ms = *msp;
  kmerConnector* kc = *kcp;
  if (ms->status->addExistingTraceStatus == 0 || ms->status->addExistingTraceStatus == 3){
    if (kc->idflags){
      destroyTIdList(&ms->status->traceSet->starting);
      destroyTIdList(&ms->status->traceSet->current);
      tIdList* fst = traceFirst(&kc->idflags);
      if (fst){
        ms->status->traceSet->starting = copyTIdList(fst);
        ms->status->traceSet->current  = copyTIdList(fst);
        ms->status->extendMeUp = true;
        ms->status->addExistingTraceStatus = 1;
        D_(2, "Found an existing trace: 0 to 1\n");
      }
      else if (ms->status->isFirst){
        ms->status->isFirst = false;
        ms->status->extendMeUp = false;
        ms->status->traceSet->starting = copyTIdList(kc->idflags);
        ms->status->traceSet->current  = copyTIdList(kc->idflags);
        //TODO: Copy all posintrace
        ms->status->addExistingTraceStatus = 1;
        D_(2, "Inside, not extending: 0 to 1\n");
      }
      else{
        ms->status->addExistingTraceStatus = 3;
        ms->status->extendMeUp = false;
        ms->status->extendMeDn = false;
        D_(2, "Found existing trace, not first: to 3\n");
      }
      destroyTIdList(&fst);
    }
  }
  else if (ms->status->addExistingTraceStatus == 1){
    trimUncompatibleTraces(msp, &ms->status->traceSet, kcp);
    if (!ms->status->traceSet->starting){
      if (ms->status->extendMe){
        ms->status->addExistingTraceStatus = 2;
        intersectTIdLists(&ms->status->traceSet->starting, ms->status->extendMe);
        destroyTIdList(&ms->status->traceSet->current);
        ms->status->traceSet->current = copyTIdList(ms->status->extendMe);
        ms->status->extendMeDn = true;
        D_(2, "Extending existing traces: 1 to 2\n");
      }
      else{
        ms->status->addExistingTraceStatus = 3;
        ms->status->extendMeUp = false;
        ms->status->extendMeDn = false;
        destroyTIdList(&ms->status->traceSet->starting);
        destroyTIdList(&ms->status->traceSet->current);
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
      destroyTIdList(&ms->status->traceSet->starting);
      destroyTIdList(&ms->status->traceSet->current);
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
  destroyTIdList(&ms->status->extendMe);
  ms->status->extendMe = traceLast(result->idflags);
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
  //Init rnd generator
  time_t t;
  srand((unsigned) time(&t));
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
  return result;
}



#endif // KMER_H_INCLUDED
