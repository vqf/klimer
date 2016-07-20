#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifndef NOKMER
#define NOKMER 0xFFFFFFFF
#endif /* NOKMER */

#ifndef DEBUG
#define DEBUG 0
#endif /* DEBUG */

#ifndef BOOL
#define BOOL
typedef enum { false, true } bool;
#endif /* BOOL */

#define ENCLOSE(a) printf("before\n"); a; printf("after\n");


typedef struct iarr{
  uint8_t n;
  struct iarr* next;
} int8Array;

typedef struct kC{
  uint32_t dest;
  uint16_t n; // Number of events supporting connection
  uint8_t flags; // Misc info
  int8Array* id; // Only used if we have to follow ambiguous reads
  struct kC* next;
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
  uint8_t cId;      // Current id for kmerConnector
  kcLL* trace;      // Keeps track of the related kmers in the current trace
  uint32_t current; // Last kmer read
  bool isFirst;     // Is it the first in a trace?
} msStatus;

typedef struct ms{ //Used to store all kmer combinations
  baseInfo* bi;
  uint32_t nPos; // Number of positions
  uint32_t nBytes; // Number of bytes in kmerArray
  msStatus* status;
  uint8_t kmerSize;
  kmerConnector* kmerArray;
} memstruct;

void resetTrace(memstruct**);
void summarize (memstruct*);
void resetKcLL(kcLL**);
void printKmerConnector(kmerConnector*, char*);
void printKc(memstruct*);

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
  return result;
}

void destroyMs(memstruct* ms){
  resetKcLL(&ms->status->trace);
  free(ms->status);
  free(ms->bi->valBase);
  free(ms->bi->baseVal);
  free(ms);
  ms = NULL;
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

void setKmerArray(memstruct* ms, uint8_t wlength){
  /*
   *Initializes and zeroes a block of memory with the values for each
   *k-mer.
  */
  ms->kmerSize = wlength;
  ms->nPos = power(ms->bi->nbases, wlength);
  resetTrace(&ms);
  uint32_t nb = sizeof(kmerConnector) * ms->nPos;
  ms->nBytes = nb;
  if (DEBUG){
    fprintf(stderr, "Allocating %lu bytes...\n", (long unsigned int) nb);
  }
  ms->kmerArray = (kmerConnector*) calloc(ms->nPos, sizeof(kmerConnector));
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

void pos2seq (memstruct* ms, uint32_t pos, char* result){
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
  uint16_t result = ms->kmerArray[pos].n;
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
void addCount(memstruct** msp, uint32_t pos){
  memstruct *ms = *msp;
  if (pos < 0 || pos > ms->nPos){
    //fprintf(stderr, "Too much\n");
    return;
  }
  //fprintf(stderr, "Add to %lu\n", (long unsigned int) pos);
  if (ms->kmerArray[pos].n <= 0xFFFF) ms->kmerArray[pos].n++;
}

/*
 *With relationships
 */


/*
 * int8Array functions
 */
void printArr(int8Array* a){
  if (a == NULL){
    printf("Null int8Array\n");
    return;
  }
  printf("%d", a->n);
  while(a->next != NULL){
    a = a->next;
    printf(", %d", a->n);
  }
  printf("\n");
}

int8Array* newInt8Array(uint8_t val){
  int8Array* result = (int8Array*) calloc(1, sizeof(int8Array));
  result->n = val;
  result->next = NULL;
  return result;
}

void int8push(int8Array** Iarr, uint8_t val){
  if (!Iarr){
    return;
  }
  printf("-- %d -- %p -- \n", val, Iarr);
  int8Array* arr = *Iarr;
  if (arr){
    int8Array* nxt = newInt8Array(val);
    int8Array* p = arr;
    while(p->next){
      p = p->next;
    }
    p->next = nxt;
  }
  else{
    arr = newInt8Array(val);
    *Iarr = arr;
  }
}

bool isInInt8Array(int8Array* arr, uint8_t val){
  bool result = false;
  if (arr->n == val){
    result = true;
  }
  while (!result && arr->next != NULL){
    arr = arr->next;
    if (arr->n == val){
      result = true;
    }
  }
  return result;
}

void delInt8Array(int8Array** todel){
  while (*todel != NULL){
    int8Array* tmp = *todel;
    *todel = (*todel)->next;
    free(tmp);
  }
}



/*
 * int8Array functions //
 */

void nextId(memstruct* ms){
  if (ms->status->cId == 0xFF){
    if (DEBUG){
      fprintf(stderr, "%s\n", "Overflow in cId");
    }
    return;
  }
  ms->status->cId++;
}


kcLL* newKcPointer(kmerConnector* newkc){
  kcLL* result = (kcLL*) calloc(1, sizeof(kcLL));
  result->kc = newkc;
  result->next = NULL;
  result->last = NULL;
  return result;
}

void resetKcLL(kcLL** llp){
  kcLL* ll = *llp;
  while (ll != NULL){
    kcLL** todel = &ll;
    ll = ll->next;
    free(*todel);
    *todel = NULL;
  }
}

void kcpush (kcLL** llp, kmerConnector** newkcp){
  kmerConnector* newkc = *newkcp;
  kcLL* ll = *llp;
  if (ll){
    kcLL* tmp = ll->last;
    tmp->next = newKcPointer(newkc);
    ll->last = tmp->next;
  }
  else{
    ll = newKcPointer(newkc);
    ll->last = ll;
    *llp = ll;
  }
}

void printKc (memstruct* ms){
  kcLL* ll = ms->status->trace;
  char* seq = (char*) calloc(12, sizeof(char));
  pos2seq(ms, ll->kc->dest, seq);
  printf("-----KmerConnector------\n");
  printf("Dest: %s\n", seq);
  free(seq);
  seq = NULL;
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

void resetTrace(memstruct** msp){
  memstruct* ms = *msp;
  kcLL** tmpp = &ms->status->trace;
  kcLL* tmp = *tmpp;
  if (tmp){
    int8push(&tmp->kc->id, ms->status->cId);
    while(tmp->next){
      tmp = tmp->next;
      int8push(&tmp->kc->id, ms->status->cId);
    }
  }
  ms->status->cId = 0;
  ms->status->isFirst = true;
  ms->status->current = NOKMER;
  resetKcLL(&ms->status->trace);
}

  


/*
 * kmerConnector functions
 */

kmerConnector* newKmerConnector(uint32_t to){
  kmerConnector* result = (kmerConnector*) calloc(1, sizeof(kmerConnector));
  result->dest  = to;
  result->n     = 0;
  result->flags = 0;
  result->id    = NULL;
  result->next  = NULL;
  return result;
}

void delConnector(kmerConnector* kc){
  while (kc != NULL){
    kmerConnector* tmp = kc;
    kc = kc->next;
    delInt8Array(&tmp->id);
    free(tmp);
  }
  kc = NULL;
}

kmerConnector* getConnector(memstruct** msp, uint32_t from, uint32_t to){
  memstruct *ms = *msp;
  kmerConnector* result = &ms->kmerArray[from];
  if (result->n == 0){ 
    result->dest = to;
    return result;
  }
  else{
    while (result->dest != to && result->next != NULL){
      result = result->next;
    }
    if (result->dest != to){
      kmerConnector* nkc = newKmerConnector(to);
      result->next = nkc;
      result = result->next;
      nextId(ms);
    }
  }
  //kcLL* kc = ms->status->trace;
  kcpush(&ms->status->trace, &result);
  
  return result;
}

void addRelationship(memstruct** msp, uint32_t to){
  memstruct* ms = *msp;
  uint32_t from = ms->status->current;
  if (ms->status->isFirst){
    ms->status->isFirst = false;
    ms->status->current = to;
    return;
  }
  if (from == NOKMER){
    resetTrace(msp);
  }
  kmerConnector* kc = getConnector(msp, from, to);
  kc->n++;
  ms->status->current = to;
}

void printKmerConnector(kmerConnector* kc, char* seq){
  printf("Dest: %s\nN: %d\nFlags: %d\nIds: ", seq, kc->n, kc->flags);
}

void printKmerConnectors(memstruct* ms, uint32_t pos){
  kmerConnector* kc = &ms->kmerArray[pos];
  while (kc != NULL){
    char* seq = (char*) calloc(ms->kmerSize + 1, sizeof(char));
    pos2seq(ms, kc->dest, seq);
    printKmerConnector(kc, seq);
    printArr(kc->id);
    free(seq);
    kc = kc->next;
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

void readIn(memstruct* ms, char* fname){
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
}

void summarize (memstruct* ms){
  uint32_t i = 0;
  for (i = 0; i < ms->nPos; i++){
    if (ms->kmerArray[i].n > 0){
      char* seq = (char*) calloc(ms->kmerSize + 1, sizeof(char));
      pos2seq(ms, i, seq);
      printf("----\n%s\n", seq);
      free(seq);
      seq = NULL;
      printKmerConnectors(ms, i);
    }
  }
}


