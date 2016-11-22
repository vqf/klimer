#ifndef SNUMBERLIST_H_INCLUDED
#define SNUMBERLIST_H_INCLUDED

#ifndef LUI
#define LUI long unsigned int
#endif /* BOOL */

#ifndef BOOL
#define BOOL
typedef enum { false, true } bool;
#endif /* BOOL */


// Trace flags
#define RESERVED       0x80
#define FIRST_IN_TRACE 0x08
#define LAST_IN_TRACE  0x04
#define CIRCULAR       0x02
#define IN_USE         0x01


// Set, unset and check flags in traces
// a -> tId object
// b -> flag to operate

#define SET(a, b) (((a)->trace.flag = (a)->trace.flag | b))
#define UNSET(a, b) ((a)->trace.flag = (a)->trace.flag ^ b)
#define IS(a, b) ((a)->trace.flag & b)


#define LISTTYPE uint32_t
#define LISTTYPESYMBOL %lu
#define MAXLISTTYPE 0xFFFFFFFF

#include "kmer.h"

typedef struct tId{
  LISTTYPE n;
  uint8_t flag;
  uint32_t nReads; // Number of events supporting connection
  void* circular; // Only used when a trace references a kmer more than once
} tId;


typedef struct iarr{
  tId trace;
  struct iarr* next;
} tIdList;

typedef struct traceLL{
  tIdList* tidl;
  struct traceLL* next;
} traceLL;

typedef traceLL traceVessel; //TraceVessel only takes pointers and its tidl is
                             // not freed

void printTraceLL(traceLL*);
traceVessel* _getTraces(tIdList**, tIdList*);

uint32_t tIdListLength(tIdList** tidp){
  tIdList* tidl = *tidp;
  uint32_t result = 0;
  while (tidl){
    result++;
    tidl = tidl->next;
  }
  return result;
}

void printTIdList(tIdList* a){
  if (!a){
    printf("Null tIdList\n");
    return;
  }
  printf("[%d, %.1x] (%u)", a->trace.n, a->trace.flag, (short unsigned) a->trace.nReads);
  while(a->next){
    a = a->next;
    printf(", [%d, %.1x] (%u)", a->trace.n, a->trace.flag, (short unsigned) a->trace.nReads);
  }
  printf("\n");
}


/*
 * New trace list
 * trace.circular is a null void pointer, which will be cast into
 * a kcLL* when necessary.
 */
tIdList* newTIdList(LISTTYPE val){
  tIdList* result = (tIdList*) calloc(1, sizeof(tIdList));
  result->trace.n = val;
  result->trace.flag = 0;
  result->next = NULL;
  result->trace.circular = NULL;
  return result;
}

void destroyTIdList(tIdList** todel, void (*funcDestroyCirc)(void**)){
  while (*todel){
    tIdList* tmp = *todel;
    *todel = (*todel)->next;
    if (tmp->trace.circular){
      funcDestroyCirc(&tmp->trace.circular);
    }
    free(tmp);
  }
}


traceLL* newTraceLL(){
  traceLL* result = (traceLL*) calloc(1, sizeof(traceLL));
  result->tidl = NULL;
  result->next = NULL;
  return result;
}


LISTTYPE maxInList(tIdList* l){
  if (!l){
    return 0;
  }
  LISTTYPE result = l->trace.n;
  while(l->next){
    l = l->next;
    if (l->trace.n > result){
      result = l->trace.n;
    }
  }
  return result;
}

void insertInTIdList(tIdList** Iarr, LISTTYPE val, void (*callback)(void**)){
  if (!*Iarr){
    *Iarr = newTIdList(val);
    return;
  }
  tIdList* arr = *Iarr;
  if (arr){
    tIdList* nxt = newTIdList(val);
    tIdList* p = arr;
    if (val < p->trace.n){ //Unshift
      nxt->next = *Iarr;
      *Iarr = nxt;
    }
    else{
      while(p->next && p->next->trace.n < val){
        p = p->next;
      }
      if (p->next){ //Insert
        if (p->trace.n == val || p->next->trace.n == val){
          destroyTIdList(&nxt, callback);
          return;
        }
        tIdList* cnt = p->next;
        p->next = nxt;
        nxt->next = cnt;
      }
      else{ //Push
        if (p->trace.n == val){
          destroyTIdList(&nxt, callback);
          return;
        }
        p->next = nxt;
      }
    }
  }
  else{
    arr = newTIdList(val);
    *Iarr = arr;
  }
}

tIdList* copyTIdList(tIdList* tocopy, void (*callback)(void**)){
  tIdList* result = NULL;
  if (!tocopy){
    return result;
  }
  tIdList* tmp = tocopy;
  while(tmp){
    insertInTIdList(&result, tmp->trace.n, callback);
    tmp = tmp->next;
  }
  return result;
}

void mergeTIdLists(tIdList** l1p, tIdList* l2, void (*callback)(void**)){
  if (l1p && *l1p){
    while (l2){
      insertInTIdList(l1p, l2->trace.n, callback);
      l2 = l2->next;
    }
  }
  else if (l2){
    *l1p = copyTIdList(l2, callback);
  }
}


void delTIdFromList(tIdList** parr, LISTTYPE val, void (*callback)(void**)){
  tIdList* arr = *parr;
  tIdList* todel = NULL;
  if (arr->trace.n == val){ //Shift
    todel = arr;
    *parr = arr->next;
    todel->next = NULL;
  }
  else if (arr->trace.n < val){
    while (arr->next && arr->next->trace.n < val){
      arr = arr->next;
    }
    if (arr->next && arr->next->trace.n == val){
      todel = arr->next;
      arr->next = arr->next->next;
      todel->next = NULL;
    }
  }
  destroyTIdList(&todel, callback);
}

void intersectTIdLists(tIdList** l1p, tIdList* l2, void (*callback)(void**)){
  /*
  Both lists must be sorted in increasing order
  NULL l2 lists are transparent
  */
  tIdList* l1 = *l1p;
  tIdList* tmp = l1;
  if (!l2){
    destroyTIdList(l1p, callback);
    *l1p = NULL;
  }
  else{
    while (tmp){
      while (l2 && (l2->trace.n < tmp->trace.n)){
        l2 = l2->next;
      }
      if (l2 && l2->trace.n == tmp->trace.n){
        tmp = tmp->next;
      }
      else{
        tIdList* tmp2 = tmp->next;
        delTIdFromList(l1p, tmp->trace.n, callback);
        tmp = tmp2;
      }
    }
  }
}


bool isInTIdList(tIdList* arr, LISTTYPE val){
  bool result = false;
  while (!result && arr){
    if (arr->trace.n == val){
      result = true;
    }
    else if (arr->trace.n > val){
      return false;
    }
    arr = arr->next;
  }
  return result;
}

bool isIncluded(tIdList** universe, tIdList** subset){
  tIdList* a1 = *universe;
  tIdList* a2 = *subset;
  if (!a1 || !a2) return false;
  if (a1->trace.n > a2->trace.n) return false;
  while (a2){
    while (a1 && a1->trace.n < a2->trace.n){
      a1 = a1->next;
    }
    if (a1 && a2 && (a1->trace.n != a2->trace.n)){
      return false;
    }
    a2 = a2->next;
  }
  return true;
}


//Pushes trace copy into the last traceLL
traceLL* lastTraceLL(traceLL** tllp){
  traceLL* result = *tllp;
  while (result->next){
    result = result->next;
  }
  return result;
}

void pushTrace(traceLL** tllp, tIdList* toadd, void (*callback)(void**)){
  traceLL* nxt = lastTraceLL(tllp);
  insertInTIdList(&nxt->tidl, toadd->trace.n, callback);
}
//Adds traceLL at the end
void nextTraceLL(traceLL** tllp){
  traceLL* nxt = lastTraceLL(tllp);
  nxt->next = newTraceLL();
}




void destroyTraceLL(traceLL** tll, void (*funcDestroyCirc)(void**)){
  traceLL* tmp = *tll;
  while (tmp){
    traceLL* nxt = tmp->next;
    destroyTIdList(&tmp->tidl, funcDestroyCirc);
    free(tmp);
    tmp = nxt;
  }
}

void destroyTraceVessel(traceVessel** tvp){
  traceVessel* tmp = *tvp;
  while (tmp){
    traceVessel* nxt = tmp->next;
    free(tmp);
    tmp = nxt;
  }
}


tIdList* traceLast(tIdList* t, void (*callback)(void**)){
  tIdList* result = NULL;
  tIdList* tmp = t;
  while (tmp){
    if (IS(tmp, LAST_IN_TRACE)){
      insertInTIdList(&result, tmp->trace.n, callback);
    }
    tmp = tmp->next;
  }
  return result;
}

bool isTraceFirst(tIdList** t, tIdList* which, void (*callback)(void**)){
  if (!t || !which) return false;
  bool result = true;
  traceVessel* tmpor = _getTraces(t, which);
  traceVessel* tmp = tmpor;
  while (tmp && tmp->tidl){
    if (!IS(tmp->tidl, FIRST_IN_TRACE)){
      result = false;
      destroyTraceVessel(&tmpor);
      return result;
    }
    tmp = tmp->next;
  }
  destroyTraceVessel(&tmpor);
  return result;
}

bool isTraceLast(tIdList** t, tIdList* which, void (*callback)(void**)){
  if (!t || !which) return false;
  bool result = true;
  traceVessel* tmpor = _getTraces(t, which);
  traceVessel* tmp = tmpor;
  while (tmp && tmp->tidl){
    if (!IS(tmp->tidl, LAST_IN_TRACE)){
      result = false;
      destroyTraceVessel(&tmpor);
      return result;
    }
    tmp = tmp->next;
  }
  destroyTraceVessel(&tmpor);
  return result;
}

tIdList* traceFirst(tIdList** tp, void (*callback)(void**)){
  tIdList* t = *tp;
  tIdList* result = NULL;
  tIdList* tmp = t;
  while (tmp){
    if (IS(tmp, FIRST_IN_TRACE)){
      insertInTIdList(&result, tmp->trace.n, callback);
    }
    tmp = tmp->next;
  }
  return result;
}

void pushTraceInVessel(traceVessel** tvp, tIdList** tp){
  traceVessel* ptr = *tvp;
  tIdList* el = *tp;
  while (ptr->next){
    ptr = ptr->next;
  }
  ptr->tidl = el;
  ptr->next = (traceVessel*) newTraceLL();
}

tIdList* _getTrace(tIdList** tp, LISTTYPE i){
  tIdList* tmp = *tp;
  while (tmp){
    if (tmp->trace.n == i){
      return tmp;
    }
    tmp = tmp->next;
  }
  return NULL;
}

traceVessel* _getTraces(tIdList** tp, tIdList* which){
  tIdList* tmp = *tp;
  traceVessel* result = (traceVessel*) newTraceLL();
  traceVessel* pointer = result;
  LISTTYPE i = which->trace.n;
  while (tmp){
    while (tmp->trace.n > i && which){
      which = which->next;
      if (which) i = which->trace.n;
    }
    if (tmp->trace.n == i){
      pointer->tidl = tmp;
      pointer->next = (traceVessel*) newTraceLL();
      pointer = pointer->next;
    }
    tmp = tmp->next;
  }
  return result;
}

void setCircular (tIdList** t, tIdList* which){
  if (!which) return;
  traceVessel* tmp = _getTraces(t, which);
  traceVessel* ptr = tmp;
  while (ptr->tidl){
    SET(ptr->tidl, CIRCULAR);
    ptr = ptr->next;
  }
  destroyTraceVessel(&tmp);
}

void setAsFirst (tIdList** t, tIdList* which){
  if (!which) return;
  traceVessel* tmp = _getTraces(t, which);
  traceVessel* ptr = tmp;
  while (ptr->tidl){
    SET(ptr->tidl, FIRST_IN_TRACE);
    ptr = ptr->next;
  }
  destroyTraceVessel(&tmp);
}
void setAsLast (tIdList** t, tIdList* which){
  if (!which) return;
  traceVessel* tmp = _getTraces(t, which);
  traceVessel* ptr = tmp;
  while (ptr->tidl){
    SET(ptr->tidl, LAST_IN_TRACE);
    ptr = ptr->next;
  }
  destroyTraceVessel(&tmp);
}

void setAsUsed(tIdList** t, tIdList* which){
  if (!which) return;
  traceVessel* tmp = _getTraces(t, which);
  traceVessel* ptr = tmp;
  while (ptr->tidl){
    SET(ptr->tidl, IN_USE);
    ptr = ptr->next;
  }
  destroyTraceVessel(&tmp);
}

void unsetInUse(tIdList** t, tIdList* which){
  if (!which) return;
  traceVessel* tmp = _getTraces(t, which);
  traceVessel* ptr = tmp;
  while (ptr->tidl){
    if (IS(ptr->tidl, IN_USE)) UNSET(ptr->tidl, IN_USE);
    ptr = ptr->next;
  }
  destroyTraceVessel(&tmp);
}

void printTraceLL(traceLL* toprint){
  printf("==\n");
  while (toprint){
    printTIdList(toprint->tidl);
    printf("\n");
    toprint = toprint->next;
  }
  printf("==\n");
}

tIdList* circTraces(tIdList** t, bool extending, void (*callback)(void**)){
  tIdList* result = NULL;
  tIdList* tmp = *t;
  while (tmp){
    if ((!IS(tmp, CIRCULAR) || extending) && IS(tmp, IN_USE)){
      insertInTIdList(&result, tmp->trace.n, callback);
    }
    tmp = tmp->next;
  }
  return result;
}



void unsetCircular(tIdList** t, tIdList* which){
  if (!which) return;
  traceVessel* tmp = _getTraces(t, which);
  traceVessel* ptr = tmp;
  while (ptr->tidl){
    if (IS(ptr->tidl, CIRCULAR)) UNSET(ptr->tidl, CIRCULAR);
    ptr = ptr->next;
  }
  destroyTraceVessel(&tmp);
}

void unsetAsFirst(tIdList** t, tIdList* which){
  if (!which) return;
  traceVessel* tmp = _getTraces(t, which);
  traceVessel* ptr = tmp;
  while (ptr->tidl){
    if (IS(ptr->tidl, FIRST_IN_TRACE)){
      UNSET(ptr->tidl, FIRST_IN_TRACE);
      D_(1, "Unsetting first in %lu\n", (LUI) ptr->tidl->trace.n);
    }
    ptr = ptr->next;
  }
  destroyTraceVessel(&tmp);
}
void unsetAsLast(tIdList** t, tIdList* which){
  if (!which) return;
  traceVessel* tmp = _getTraces(t, which);
  traceVessel* ptr = tmp;
  while (ptr->tidl){
    if (IS(ptr->tidl, LAST_IN_TRACE)){
      UNSET(ptr->tidl, LAST_IN_TRACE);
      D_(1, "Unsetting last in %lu\n", (LUI) ptr->tidl->trace.n);
    }
    ptr = ptr->next;
  }
  destroyTraceVessel(&tmp);
}

#endif // SNUMBERLIST_H_INCLUDED
