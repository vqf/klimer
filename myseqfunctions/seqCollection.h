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
  uint32_t cPos;
  struct mll* next;
  struct mll* last;
} kcLL;


typedef struct seqCollection{
  kcLL* trace;
  LISTTYPE traceId;
  LISTTYPE posInTrace;
  bool delme;
  struct seqCollection* next;
  struct seqCollection* down;
  char* seq;
} seqCollection;

void destroySeqCollection(seqCollection**);

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
  result->cPos = pos;
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
      kcpush(&result, &kcl->kc, kcl->ntid, kcl->cPos);
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
  result->posInTrace = 0;
  result->next = NULL;
  result->down = NULL;
  result->delme = false;
  result->seq = NULL;
  return result;
}

seqCollection* pushSeq(seqCollection** scp, kcLL** tracep, LISTTYPE posInTrace){
  // returns a pointer to the last element
  seqCollection* sc = *scp;
  if (!sc){
    sc = newSeqCollection();
    sc->trace = *tracep;
    sc->posInTrace = posInTrace;
    *scp = sc;
    return sc;
  }
  while (sc && sc->down){
    sc = sc->down;
  }
  seqCollection* dn = newSeqCollection();
  dn->trace = *tracep;
  dn->posInTrace = posInTrace;
  sc->down = dn;
  return dn;
}

seqCollection* addSeq(seqCollection** scp, kcLL** tracep){
  // returns a pointer to the last element
  seqCollection* sc = *scp;
  if (!sc){
    sc = newSeqCollection();
    sc->trace = *tracep;
    *scp = sc;
    return sc;
  }
  while (sc && sc->next){
    sc = sc->next;
  }
  seqCollection* nxt = newSeqCollection();
  nxt->trace = *tracep;
  sc->next = nxt;
  return nxt;
}

void pruneSeqCollection(seqCollection** scp){
  seqCollection* sc = *scp;
  seqCollection* psc = NULL;
  int counter = 0;
  while (sc){
    counter++;
    if (sc->delme){
      seqCollection* todel = sc;
      D_(2, "Deleting %d\n", counter);
      sc = sc->down;
      if (psc){
        psc->down = sc;
      }
      else{
        *scp = sc;
      }
      todel->down = NULL;
      todel->next = NULL;
      destroySeqCollection(&todel);
    }
    else{
      psc = sc;
      sc = sc->down;
    }
  }
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
    resetKcLL(&sc->trace);
    free(sc->seq);
    free(sc);
    sc = nxt;
  }
  scp = NULL;
}




#endif // SEQCOLLECTION_H_INCLUDED
