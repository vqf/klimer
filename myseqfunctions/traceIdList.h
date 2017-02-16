#ifndef SNUMBERLIST_H_INCLUDED
#define SNUMBERLIST_H_INCLUDED

#include "mydmacros.h"

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

//#include "kmer.h"

typedef struct tId{
  LISTTYPE n;
  uint8_t  flag;
  uint32_t nReads; // Number of events supporting connection
} tId;


typedef struct iarr{
  tId trace;
  struct iarr* posInTrace;
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
    return;
  }
  while(a){
    if (a->posInTrace){
      printf("[");
    }
    else{
      printf("{");
    }
    printf("%d %01x (%u)", a->trace.n, a->trace.flag, (short unsigned) a->trace.nReads);
    if (a->posInTrace) printf(", ");
    printTIdList(a->posInTrace);
    if (a->posInTrace){
      printf("]");
    }
    else{
      printf("}");
    }
    a = a->next;
    if (a) printf(", ");
  }
  //printf("\n");
}


/*
 * New trace list
 * trace.circular is a null void pointer, which will be cast into
 * a kcLL* when necessary.
 */
tIdList* newTIdList(LISTTYPE val){
  tIdList* result = (tIdList*) calloc(1, sizeof(tIdList));
  result->trace.n = val;
  result->trace.nReads = 0;
  result->next = NULL;
  result->posInTrace = NULL;
  return result;
}

void destroyTIdList(tIdList** todel){
  while (*todel){
    tIdList* tmp = *todel;
    *todel = (*todel)->next;
    if (tmp->posInTrace) destroyTIdList(&tmp->posInTrace);
    free(tmp);
  }
  *todel = NULL;
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

tIdList* insertInTIdList(tIdList** Iarr, LISTTYPE val){
  if (!*Iarr){
    *Iarr = newTIdList(val);
    return *Iarr;
  }
  tIdList* arr = *Iarr;
  if (arr){
    tIdList* p = arr;
    if (val < p->trace.n){ //Unshift
      tIdList* nxt = newTIdList(val);
      nxt->next = *Iarr;
      *Iarr = nxt;
      return nxt;
    }
    else{
      while(p->next && p->next->trace.n < val){
        p = p->next;
      }
      if (p->next){ //Insert
        if (p->trace.n == val || p->next->trace.n == val){
          if (p->trace.n == val) return p;
          if (p->next->trace.n == val) return p->next;
        }
        tIdList* nxt = newTIdList(val);
        tIdList* cnt = p->next;
        p->next = nxt;
        nxt->next = cnt;
        return nxt;
      }
      else{ //Push
        if (p->trace.n == val){
          return p;
        }
        tIdList* nxt = newTIdList(val);
        p->next = nxt;
        return nxt;
      }
    }
  }
  else{
    arr = newTIdList(val);
    *Iarr = arr;
    return arr;
  }
}

tIdList* addTrace(tIdList** tlp, LISTTYPE i){
  tIdList* r = insertInTIdList(tlp, i);
  r->trace.nReads++;
  return r;
}

tIdList* addPosInTrace(tIdList** ap, LISTTYPE pos){
  tIdList* r = insertInTIdList((&(*ap)->posInTrace), pos);
  r->trace.nReads++;
  return r;
}

void nextPos(tIdList** tp){
  tIdList* t = *tp;
  while (t){
    tIdList* cPos = t->posInTrace;
    while (cPos){
      cPos->trace.n++;
      cPos = cPos->next;
    }
    t = t->next;
  }
}


tIdList* copyTIdList(tIdList** tocopyp){
  tIdList* tocopy = *tocopyp;
  tIdList* result = NULL;
  if (!tocopy){
    return result;
  }
  tIdList* tmp = tocopy;
  while(tmp){
    tIdList* addtome = insertInTIdList(&result, tmp->trace.n);
    tIdList* pos = tmp->posInTrace;
    while (pos){
      insertInTIdList(&addtome->posInTrace, pos->trace.n);
      pos = pos->next;
    }
    tmp = tmp->next;
  }
  return result;
}

void mergeTIdLists(tIdList** l1p, tIdList* l2){
  if (l1p && *l1p){
    while (l2){
      insertInTIdList(l1p, l2->trace.n);
      l2 = l2->next;
    }
  }
  else if (l2){
    *l1p = copyTIdList(&l2);
  }
}


void delTIdFromList(tIdList** parr, LISTTYPE val){
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
    if (arr->next->trace.n == val && !arr->next){
      destroyTIdList(parr);
      *parr = NULL;
      return;
    }
    if (arr->next && arr->next->trace.n == val){
      todel = arr->next;
      arr->next = arr->next->next;
      todel->next = NULL;
    }
  }
  destroyTIdList(&todel);
}

void intersectTIdLists(tIdList** l1p, tIdList* l2){
  /*
  Both lists must be sorted in increasing order
  */
  tIdList* l1 = *l1p;
  tIdList* tmp = l1;
  if (!l2){
    destroyTIdList(l1p);
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
        delTIdFromList(l1p, tmp->trace.n);
        tmp = tmp2;
      }
    }
  }
}

tIdList* isInTIdList(tIdList** arrp, LISTTYPE val){
  tIdList *arr = *arrp;
  tIdList* result = NULL;
  while (!result && arr){
    if (arr->trace.n == val){
      result = arr;
    }
    arr = arr->next;
  }
  return result;
}

tIdList* isInTIdListNotInUse(tIdList** arrp, LISTTYPE val){
  tIdList *arr = *arrp;
  tIdList* result = NULL;
  while (!result && arr){
    if (!IS(arr, IN_USE) && arr->trace.n == val){
      result = arr;
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
  if (result){
    while (result->next){
      result = result->next;
    }
  }
  return result;
}

void pushTrace(traceLL** tllp, tIdList** toaddp){
  traceLL* tll = *tllp;
  tIdList* toadd = *toaddp;
  if (!tll){
    *tllp = newTraceLL();
    tll = *tllp;
    tll->tidl = toadd;
  }
  else{
    traceLL* lst = lastTraceLL(tllp);
    lst->next = newTraceLL();
    lst->next->tidl = toadd;
  }
}
//Adds traceLL at the end
void nextTraceLL(traceLL** tllp){
  traceLL* nxt = lastTraceLL(tllp);
  nxt->next = newTraceLL();
}




void destroyTraceLL(traceLL** tll){
  traceLL* tmp = *tll;
  while (tmp){
    traceLL* nxt = tmp->next;
    free(tmp);
    tmp = nxt;
  }
  *tll = NULL;
}

void destroyTraceVessel(traceVessel** tvp){
  traceVessel* tmp = *tvp;
  while (tmp){
    traceVessel* nxt = tmp->next;
    free(tmp);
    tmp = nxt;
  }
  *tvp = NULL;
}


tIdList* traceLast(tIdList** tp){
  tIdList* t = *tp;
  tIdList* result = NULL;
  tIdList* tmp = t;
  while (tmp){
    if (IS(tmp, LAST_IN_TRACE)){
      tIdList* lPos = tmp->posInTrace;
      LISTTYPE mpos = 0;
      while (lPos && lPos->next){
        lPos = lPos->next;
      }
      mpos = lPos->trace.n;
      tIdList* ptr = insertInTIdList(&result, tmp->trace.n);
      insertInTIdList(&ptr->posInTrace, mpos);
    }
    tmp = tmp->next;
  }
  return result;
}

bool isTraceFirst(tIdList** t, tIdList* which){
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

bool isTraceLast(tIdList** t, tIdList* which){
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

tIdList* traceFirst(tIdList** tp){
  tIdList* t = *tp;
  tIdList* result = NULL;
  tIdList* tmp = t;
  while (tmp){
    if (IS(tmp, FIRST_IN_TRACE)){
      if (tmp->posInTrace->trace.n > 1){
        D_(0, "First in trace, but posInTrace is %lu\n", (long unsigned) tmp->posInTrace->trace.n);
      }
      tIdList* ptr = insertInTIdList(&result, tmp->trace.n);
      insertInTIdList(&ptr->posInTrace, 1);
    }
    tmp = tmp->next;
  }
  E_(2, printTIdList(result); printf("\n"););
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

tIdList* _getTrace(tIdList** tp, LISTTYPE i, LISTTYPE pos){
  tIdList* tmp = *tp;
  while (tmp){
    if (tmp->trace.n == i){
      if (!pos) return tmp;
      tIdList* intmp = tmp->posInTrace;
      while (intmp){
        if (!pos || intmp->trace.n == pos){
          return tmp;
        }
        intmp = intmp->next;
      }
    }
    tmp = tmp->next;
  }
  return NULL;
}

tIdList* _getTraceNotInUse(tIdList** tp, LISTTYPE i, LISTTYPE pos){
  tIdList* tmp = *tp;
  while (tmp){
    if (tmp->trace.n == i && !IS(tmp, IN_USE)){
      tIdList* intmp = tmp->posInTrace;
      while (intmp){
        if (intmp->trace.n == pos && !IS(intmp, IN_USE)){
          return tmp;
        }
        intmp = intmp->next;
      }
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


bool isFirst(tIdList** t, tIdList* which){
  if (!which) return false;
  bool result = false;
  traceVessel* tmp = _getTraces(t, which);
  traceVessel* ptr = tmp;
  while (ptr->tidl){
    if (IS(ptr->tidl, FIRST_IN_TRACE)){
      result = true;
    }
    else{
      destroyTraceVessel(&tmp);
      return false;
    }
    ptr = ptr->next;
  }
  destroyTraceVessel(&tmp);
  return result;
}

bool isLast(tIdList** t, tIdList* which){
  if (!which) return false;
  bool result = false;
  traceVessel* tmp = _getTraces(t, which);
  traceVessel* ptr = tmp;
  while (ptr->tidl){
    if (IS(ptr->tidl, LAST_IN_TRACE)){
      result = true;
    }
    else{
      destroyTraceVessel(&tmp);
      return false;
    }
    ptr = ptr->next;
  }
  destroyTraceVessel(&tmp);
  return result;
}

void setAsFirst (tIdList** t, LISTTYPE id, LISTTYPE pos){
  tIdList* tmp = _getTrace(t, id, pos);
  if (tmp){
    SET(tmp, FIRST_IN_TRACE);
    tIdList* tmpos = _getTrace(&tmp->posInTrace, pos, 0);
    if (!tmpos){
      printTIdList(tmp->posInTrace);X_;
    }
    SET(tmpos, FIRST_IN_TRACE);
  }
  else{
    D_(0, "Cannot find first in trace, id %lu, pos %lu\n", (LUI) id, (LUI) pos);
  }
}

void setAsLast (tIdList** t, LISTTYPE id, LISTTYPE pos){
  tIdList* tmp = _getTrace(t, id, pos);
  if (tmp){
    SET(tmp, LAST_IN_TRACE);
    tIdList* tmpos = _getTrace(&tmp->posInTrace, pos, 0);
    SET(tmpos, LAST_IN_TRACE);
  }
  else{
    D_(0, "Cannot find last in trace\n");
  }
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


void unsetAsFirst (tIdList** t, LISTTYPE id, LISTTYPE pos){
  tIdList* tmp = _getTrace(t, id, pos);
  if (tmp){
    UNSET(tmp, FIRST_IN_TRACE);
    tIdList* tmpos = _getTrace(&tmp->posInTrace, pos, 0);
    UNSET(tmpos, FIRST_IN_TRACE);
  }
  else{
    D_(0, "Cannot find last in trace\n");
  }
}

void unsetAsLast (tIdList** t, LISTTYPE id, LISTTYPE pos){
  tIdList* tmp = _getTrace(t, id, pos);
  if (tmp){
    UNSET(tmp, LAST_IN_TRACE);
    tIdList* tmpos = _getTrace(&tmp->posInTrace, pos, 0);
    UNSET(tmpos, LAST_IN_TRACE);
  }
  else{
    D_(0, "Cannot find last in trace\n");
  }
}

/* TraceSet */

typedef struct tSet{
  LISTTYPE traceId;
  LISTTYPE startingPos; // Where each trace starts
  LISTTYPE currentPos;
  LISTTYPE upShift; // If type is 1 or 3, how many kcs until we reach the first existing one
  LISTTYPE dnShift; // If type is 2 or 3, how many kcs until we reach the last existing one
  bool conflict;    // More than one upShift is possible
  bool fixed;       // Should not be checked anymore (type 2, reached end of trace)
  uint8_t type;
  struct tSet* prev;
  struct tSet* next;
  // Type - relationship to trace
  //  0 inside, not extending
  //  1 extendUp
  //  2 extendDn)
  //  3 extendUp + extendDn
} tSet;

tSet* newTSet(LISTTYPE traceId, LISTTYPE pos){
  tSet* result = (tSet*) calloc(1, sizeof(tSet));
  result->traceId = traceId;
  result->startingPos = pos;
  result->currentPos  = pos;
  result->conflict = false;
  result->fixed = false;
  result->upShift = 0;
  result->dnShift = 0;
  result->next = NULL;
  result->prev = NULL;
  return result;
}

void destroyTSet(tSet** todelp){
  tSet* todel = *todelp;
  while (todel){
    tSet* nxt = todel->next;
    free(todel);
    todel = nxt;
  }
  *todelp = NULL;
}

void unshiftTSet(tSet** origp, LISTTYPE traceId, LISTTYPE pos){
  tSet* orig = *origp;
  tSet* toadd = newTSet(traceId, pos);
  toadd->next = orig;
  if (orig) orig->prev = toadd;
  *origp = toadd;
}

void exciseTSet (tSet** headp, tSet** todelp){
  tSet* todel = *todelp;
  if (!todel) return;
  tSet* p = todel->prev;
  tSet* n = todel->next;
  if (!p && !n){
    destroyTSet(headp);
    return;
  }
  if (p){
    p->next = n;
  }
  else{ // shift
    *headp = todel->next;
    todel->next = NULL;
    destroyTSet(&todel);
    (*headp)->prev = NULL;
    return;
  }
  if (n){
    n->prev = p;
  }
  else{  //pop
    if (p){
      p->next = NULL;
      destroyTSet(&todel);
    }
    return;
  }
  todel->next = NULL;
  destroyTSet(&todel);
  todelp = &n;
}

bool _isLower(tSet** one, tSet** two){
  if ((*one)->traceId > (*two)->traceId){
    return false;
  }
  else if ((*one)->traceId == (*two)->traceId){
    if ((*one)->startingPos > (*two)->startingPos){
      return false;
    }
  }
  return true;
}

tSet* copyTSet(tSet* t){
  tSet* result = newTSet(t->traceId, t->startingPos);
  result->currentPos = t->currentPos;
  result->type = t->type;
  result->upShift = t->upShift;
  result->dnShift = t->dnShift;
  return result;
}

void printTSet(tSet* t){
  if (!t) return;
  printf("\n=====\n");
  while (t){
    printf("[%lu (%u), {%lu, %lu}]\n", (LUI) t->traceId,
                                       t->type,
                                       (LUI) t->startingPos,
                                       (LUI) t->currentPos);
    t = t->next;
  }
  printf("=====\n");
}

#endif // SNUMBERLIST_H_INCLUDED
