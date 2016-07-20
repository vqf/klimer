#ifndef SNUMBERLIST_H_INCLUDED
#define SNUMBERLIST_H_INCLUDED



#endif // SNUMBERLIST_H_INCLUDED

#ifndef BOOL
#define BOOL
typedef enum { false, true } bool;
#endif /* BOOL */


#define LISTTYPE uint8_t
#define MAXLISTTYPE 0xFF

typedef struct iarr{
  LISTTYPE n;
  struct iarr* next;
} snumberList;


void printSnumberList(snumberList* a){
  if (!a){
    printf("Null snumberList\n");
    return;
  }
  printf("%d", a->n);
  while(a->next){
    a = a->next;
    printf(", %d", a->n);
  }
  printf("\n");
}

snumberList* newSnumberList(LISTTYPE val){
  snumberList* result = (snumberList*) calloc(1, sizeof(snumberList));
  result->n = val;
  result->next = NULL;
  return result;
}

void destroySnumberList(snumberList** todel){
  while (*todel != NULL){
    snumberList* tmp = *todel;
    *todel = (*todel)->next;
    free(tmp);
  }
}

LISTTYPE maxInList(snumberList* l){
  if (!l){
    return 0;
  }
  LISTTYPE result = l->n;
  while(l->next){
    l = l->next;
    if (l->n > result){
      result = l->n;
    }
  }
  return result;
}

void insertInSnumberList(snumberList** Iarr, uint8_t val){
  if (!*Iarr){
    *Iarr = newSnumberList(val);
    return;
  }
  snumberList* arr = *Iarr;
  if (arr){
    snumberList* nxt = newSnumberList(val);
    snumberList* p = arr;
    if (val < p->n){ //Unshift
      nxt->next = p;
      *Iarr = nxt;
    }
    else{
      while(p->next && p->next->n < val){
        p = p->next;
      }
      if (p->next){ //Insert
        if (p->n == val || p->next->n == val){
          return;
        }
        snumberList* cnt = p->next;
        p->next = nxt;
        nxt->next = cnt;
      }
      else{ //Push
        if (p->n == val){
          return;
        }
        p->next = nxt;
      }
    }
  }
  else{
    arr = newSnumberList(val);
    *Iarr = arr;
  }
}

snumberList* copySnumberList(snumberList* tocopy){
  snumberList* result = NULL;
  if (!tocopy){
    return result;
  }
  LISTTYPE n = tocopy->n;
  result = newSnumberList(n);
  snumberList* tmp = tocopy;
  while(tmp->next){
    insertInSnumberList(&result, tmp->next->n);
    tmp = tmp->next;
  }
  return result;
}

void mergeSnumberLists(snumberList** l1p, snumberList* l2){
  if (l1p && *l1p){
    while (l2){
      insertInSnumberList(l1p, l2->n);
      l2 = l2->next;
    }
  }
  else if (l2){
    *l1p = copySnumberList(l2);
  }
}


void delSnumberFromList(snumberList** parr, LISTTYPE val){
  snumberList* arr = *parr;
  snumberList* todel = NULL;
  if (arr->n == val){ //Shift
    todel = arr;
    *parr = arr->next;
    todel->next = NULL;
  }
  else if (arr->n < val){
    while (arr->next && arr->next->n < val){
      arr = arr->next;
    }
    if (arr->next && arr->next->n == val){
      todel = arr->next;
      arr->next = arr->next->next;
      todel->next = NULL;
    }
  }
  destroySnumberList(&todel);
}

void intersectSnumberLists(snumberList** l1p, snumberList* l2){
  /*
  Both lists must be sorted in increasing order
  NULL l2 lists are transparent
  */
  snumberList* l1 = *l1p;
  snumberList* tmp = l1;
  if (!l2){
    destroySnumberList(l1p);
    *l1p = NULL;
  }
  else{
    while (tmp){
      while (l2 && (l2->n < tmp->n)){
        l2 = l2->next;
      }
      if (l2 && l2->n == tmp->n){
        tmp = tmp->next;
      }
      else{
        snumberList* tmp2 = tmp->next;
        delSnumberFromList(l1p, tmp->n);
        tmp = tmp2;
      }
    }
  }
}


bool isInSnumberList(snumberList* arr, LISTTYPE val){
  bool result = false;
  while (!result && arr){
    if (arr->n == val){
      result = true;
    }
    else if (arr->n > val){
      return false;
    }
    arr = arr->next;
  }
  return result;
}


