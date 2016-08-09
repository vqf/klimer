#ifndef SNUMBERLIST_H_INCLUDED
#define SNUMBERLIST_H_INCLUDED


#ifndef BOOL
#define BOOL
typedef enum { false, true } bool;
#endif /* BOOL */


// Trace flags

#define FIRST_IN_TRACE 0x80
#define LAST_IN_TRACE  0x40
#define IN_USE         0x01

// Set, unset and check flags in traces
// a -> tId object
// b -> flag to operate

#define SET(a, b) (((a)->trace.flag = (a)->trace.flag | b))
#define UNSET(a, b) ((a)->trace.flag = (a)->trace.flag ^ b)
#define IS(a, b) ((a)->trace.flag & b)
/*
#define SET_AS_FIRST(a) ((a)->trace.flag = (a)->trace.flag | FIRST_IN_TRACE)
#define UNSET_AS_FIRST(a) ((a)->trace.flag = (a)->trace.flag ^ FIRST_IN_TRACE)
#define IS_FIRST(a) ((a)->trace.flag & FIRST_IN_TRACE)

#define SET_AS_LAST(a) ((a)->trace.flag = (a)->trace.flag | LAST_IN_TRACE)
#define UNSET_AS_LAST(a) ((a)->trace.flag = (a)->trace.flag ^ LAST_IN_TRACE)
#define IS_LAST(a) ((a)->trace.flag & LAST_IN_TRACE)

#define SET_IN_USE(a) ((a)->trace.flag = (a)->trace.flag | IN_USE)
#define UNSET_IN_USE(a) ((a)->trace.flag = (a)->trace.flag ^ IN_USE)
#define IS_IN_USE(a) ((a)->trace.flag & IN_USE)
*/
//


#define LISTTYPE uint8_t
#define MAXLISTTYPE 0xFF

#include "kmer.h"


typedef struct tId{
  LISTTYPE n;
  uint8_t flag;
  struct mll* circular; // Only used when a trace references a kmer more than once
} tId;


typedef struct iarr{
  tId trace;
  struct iarr* next;
} tIdList;


void printTIdList(tIdList* a){
  if (!a){
    printf("Null tIdList\n");
    return;
  }
  printf("[%d, %u]", a->trace.n, (short unsigned) a->trace.flag);
  while(a->next){
    a = a->next;
    printf(", [%d, %u]", a->trace.n, (short unsigned) a->trace.flag);
  }
  printf("\n");
}

tIdList* newTIdList(LISTTYPE val){
  tIdList* result = (tIdList*) calloc(1, sizeof(tIdList));
  result->trace.n = val;
  result->trace.flag = 0;
  result->next = NULL;
  result->trace.circular = NULL;
  return result;
}

void destroyTIdList(tIdList** todel){
  while (*todel != NULL){
    tIdList* tmp = *todel;
    *todel = (*todel)->next;
    resetKcLL(&tmp->trace.circular);
    free(tmp);
  }
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

void insertInTIdList(tIdList** Iarr, uint8_t val){
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
          destroyTIdList(&nxt);
          return;
        }
        tIdList* cnt = p->next;
        p->next = nxt;
        nxt->next = cnt;
      }
      else{ //Push
        if (p->trace.n == val){
          destroyTIdList(&nxt);
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

tIdList* copyTIdList(tIdList* tocopy){
  tIdList* result = NULL;
  if (!tocopy){
    return result;
  }
  LISTTYPE n = tocopy->trace.n;
  result = newTIdList(n);
  tIdList* tmp = tocopy;
  while(tmp->next){
    insertInTIdList(&result, tmp->next->trace.n);
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
    *l1p = copyTIdList(l2);
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
  NULL l2 lists are transparent
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

tIdList* traceLast(tIdList* t){
  tIdList* result = NULL;
  tIdList* tmp = t;
  while (tmp){
    if (IS(tmp, LAST_IN_TRACE)){
      insertInTIdList(&result, tmp->trace.n);
    }
    tmp = tmp->next;
  }
  return result;
}

tIdList* traceFirst(tIdList* t){
  tIdList* result = NULL;
  tIdList* tmp = t;
  while (tmp){
    if (IS(tmp, FIRST_IN_TRACE)){
      insertInTIdList(&result, tmp->trace.n);
    }
    tmp = tmp->next;
  }
  return result;
}

void setAsFirst (tIdList** t){
  tIdList* tmp = *t;
  while (tmp){
    SET(tmp, FIRST_IN_TRACE);
    tmp = tmp->next;
  }
}
void setAsLast (tIdList** t){
  tIdList* tmp = *t;
  while (tmp){
    SET(tmp, LAST_IN_TRACE);
    tmp = tmp->next;
  }
}
void unsetAsFirst(tIdList** t){
  tIdList* tmp = *t;
  while (tmp){
    if (IS(tmp, FIRST_IN_TRACE)) UNSET(tmp, FIRST_IN_TRACE);
    tmp = tmp->next;
  }
}
void unsetAsLast(tIdList** t){
  tIdList* tmp = *t;
  while (tmp){
    if (IS(tmp, LAST_IN_TRACE)) UNSET(tmp, LAST_IN_TRACE);
    tmp = tmp->next;
  }
}

#endif // SNUMBERLIST_H_INCLUDED
