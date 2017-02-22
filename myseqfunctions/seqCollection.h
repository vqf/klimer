#ifndef SEQCOLLECTION_H_INCLUDED
#define SEQCOLLECTION_H_INCLUDED

#include "mydmacros.h"


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


typedef struct seqCollection{
  kcLL* trace;
  LISTTYPE traceId;
  struct seqCollection* next;
  struct seqCollection* down;
} seqCollection;



void resetKcLL(kcLL** llp){
  kcLL* ll = *llp;
  while (ll){
    kcLL* todel = ll;
    ll = ll->next;
    free(todel);
  }
  llp = NULL;
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







seqCollection* newSeqCollection(){
  seqCollection* result = calloc(1, sizeof(seqCollection));
  result->trace = NULL;
  result->next = NULL;
  result->down = NULL;
  return result;
}

seqCollection* pushSeq(seqCollection** scp, kcLL** tracep){
  // returns a pointer to the last element
  seqCollection* sc = *scp;
  while (sc && sc->next){
    sc = sc->next;
  }
  seqCollection* nxt = newSeqCollection();
  nxt->trace = *tracep;
  sc->next = nxt;
  return nxt;
}

void clearTraceUse(kcLL** kclp){
  kcLL* kcl = *kclp;
  while (kcl){
    kmerConnector* kc = kcl->kc;
    if (kc && ISKC(kc, IN_USE)){
      UNSETKC(kc, IN_USE);
    }
    kcl = kcl->next;
  }
}

void destroySeqCollection(seqCollection** scp){
  seqCollection* sc = *scp;
  while(sc){
    if (sc->down){
      destroySeqCollection(&sc->down);
    }
    seqCollection* nxt = sc->next;
    free(sc);
    sc = nxt;
  }
  scp = NULL;
}



#endif // SEQCOLLECTION_H_INCLUDED
