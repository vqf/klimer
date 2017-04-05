#include <sys/types.h>
#include <stdio.h>
#include "IOLayer.h"


#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif /* _LARGEFILE64_SOURCE */

#define NAMESIZE 0xFFFF
#define MAXREAD 0xFFFF


#ifndef BOOL
#define BOOL
typedef enum { false, true } bool;
#endif /* BOOL */


typedef struct fqr{
  fpointer* fp;
  char* fileName;
  char* seqName;
  char* fastaSeq;
  char* comments;
  char* quals;
  uint32_t pos; //Position in buffer
  uint32_t seqlen;
} fastqReader;

/*uint64_t fSize(fastqReader** fqp){
  fastqReader* fqr = *fqp;
  uint64_t result = 0;
  FILE * f = fopen(fqr->fileName, "rb");
  fseek(f, 0, SEEK_END);
  result = (uint64_t) ftell(f);
  fclose(f);
  return result;
}

uint8_t doTell(fastqReader** fqp){ //Percentage read in file
  fastqReader* fqr = *fqp;
  uint8_t result = 0;
  uint64_t fpos = 0;
  if (fqr->type == gzipped){
    fpos = (uint64_t) gzoffset(*(fqr->zfp));
  }
  else{
    fpos = (uint64_t) ftell(fqr->fp);
  }
  uint64_t fsz = fSize(&fqr);
  result = (uint8_t) (fpos / ((uint64_t) (fsz / 100)));
  return result;
}*/

fastqReader* initFastqReader(char* fName){
  fastqReader* result = calloc(1, sizeof(fastqReader));
  result->fileName = fName;
  result->seqName  = calloc(NAMESIZE, sizeof(char));
  result->fastaSeq = calloc(MAXREAD, sizeof(char));
  result->comments = calloc(NAMESIZE, sizeof(char));
  result->quals    = calloc(MAXREAD, sizeof(char));
  result->pos = 0;
  result->seqlen = 0;
  result->fp = newFPointer(fName);
  return result;
}

void _resetFq(fastqReader** fqp){
  fastqReader* fq = *fqp;
  fq->pos = 0;
  fq->seqlen = strlen(fq->fastaSeq);
}

bool endTrace(fastqReader* fq){
  if (!fq->seqlen || fq->pos >= fq->seqlen ){
    return true;
  }
  return false;
}

bool looksLikeFastq(char* fname){
  bool result = false;
  bool goon = true;
  fpointer* tmpfp = newFPointer(fname);
  char* line = (char*) calloc(60, sizeof(char));
  while (goon){
    char* r = myfgets(&tmpfp, &line, 60);
    if (!r) goon = false;
    char* ec = &line[0];
    if (*ec > 0x20){
      if (*ec == 0x40) result = true;
      goon = false;
    }
  }
  free(line);
  destroyFPointer(&tmpfp);
  return result;
}

bool getNextFqRead(fastqReader** fqp){
  fastqReader* fq = *fqp;
  bool result = false;
  if (fq == NULL){
    fprintf(stderr, "Non-existent fastqReader\n");
    exit(0);
  }
  if (myfgets(&fq->fp, &fq->seqName,  NAMESIZE) &&
      myfgets(&fq->fp, &fq->fastaSeq,  MAXREAD) &&
      myfgets(&fq->fp, &fq->comments, NAMESIZE) &&
      myfgets(&fq->fp, &fq->quals,     MAXREAD)){
    result = true;
  }
  _resetFq(&fq);
  return result;
}

char getNextFqBase(fastqReader** fqp){
  fastqReader* fq = *fqp;
  char result = '\0';
  if (!fq->seqlen || fq->pos >= fq->seqlen ){
    return result;
  }
  result = fq->fastaSeq[fq->pos];
  fq->pos++;
  return result;
}

char getFqQual(fastqReader** fqp){
  fastqReader* fq = *fqp;
  char result = '\0';
  if (!fq->seqlen || fq->pos >= fq->seqlen ){
    return result;
  }
  result = fq->quals[fq->pos];
  return result;
}

void destroyFastqReader(fastqReader** fqrp){
  if (fqrp && *fqrp){
    fastqReader* fqr = *fqrp;
    free(fqr->comments);
    free(fqr->fastaSeq);
    free(fqr->quals);
    free(fqr->seqName);
    destroyFPointer(&fqr->fp);
    free(fqr);
    fqrp = NULL;
  }
}
