#include <zlib.h>
#define NAMESIZE 0xFFFF
#define MAXREAD 0xFFFF

#ifndef BOOL
#define BOOL
typedef enum { false, true } bool;
#endif /* BOOL */
typedef enum { normal, gzipped } instream;

typedef struct fqr{
  FILE* fp;
  gzFile* zfp;
  char* seqName;
  char* fastaSeq;
  char* comments;
  char* quals;
  instream type;
  uint32_t pos; //Position in buffer
  uint32_t seqlen;
} fastqReader;

fastqReader* initFastqReader(char* fName){
  fastqReader* result = calloc(1, sizeof(fastqReader));
  result->seqName  = calloc(NAMESIZE, sizeof(char));
  result->fastaSeq = calloc(MAXREAD, sizeof(char));
  result->comments = calloc(NAMESIZE, sizeof(char));
  result->quals    = calloc(MAXREAD, sizeof(char));
  result->pos = 0;
  result->seqlen = 0;
  result->type = normal;
  if (strcmp(fName, "-") == 0){
    result->fp = stdin;
  }
  else{
    FILE* fip = fopen(fName, "r");
    if (!fip){
      fprintf(stderr, "Could not open %s\n", fName);
      return NULL;
    }
    uint8_t magic[2];
    fread(magic, 2, 1, fip);
    //fprintf(stderr, "Magic: %x, %x\n", magic[0], (uint8_t) magic[1]);
    if (magic[0] == 0x1f && magic[1] == 0x8b){  //gzipped
      result->type = gzipped;
      fclose(fip);
      fip = NULL;
    }
    if (fip){
      result->fp = fip;
    }
    else if (result->type == gzipped){
      gzFile zfip = gzopen(fName, "r");
      result->zfp = &zfip;
      //fprintf(stderr, "Gzipped\n");
    }
    else{
      fprintf(stderr, "Could not open %s\n", fName);
    }
  }
  return result;
}

void _resetFq(fastqReader* fq){
  fq->pos = 0;
  fq->seqlen = strlen(fq->fastaSeq);
}

bool endTrace(fastqReader* fq){
  if (!fq->seqlen || fq->pos >= fq->seqlen ){
    return true;
  }
  return false;
}

instream ftype(fastqReader* fq){
  return fq->type;
}

bool getNextFqRead(fastqReader* fq){
  bool result = false;
  if (fq == NULL){
    fprintf(stderr, "Non-existent fastqReader\n");
    exit(0);
  }
  if (fq->type == gzipped){
    if (gzgets(*fq->zfp, fq->seqName,  NAMESIZE) &&
        gzgets(*fq->zfp, fq->fastaSeq,  MAXREAD) &&
        gzgets(*fq->zfp, fq->comments, NAMESIZE) &&
        gzgets(*fq->zfp, fq->quals,     MAXREAD)
        ){
      _resetFq(fq);
      result = true;
    }
  }
  else if ( fgets(fq->seqName,  NAMESIZE, fq->fp) &&
            fgets(fq->fastaSeq,  MAXREAD, fq->fp) &&
            fgets(fq->comments, NAMESIZE, fq->fp) &&
            fgets(fq->quals,     MAXREAD, fq->fp)
            ){
    _resetFq(fq);
    result = true;
  }
  return result;
}

char getNextFqBase(fastqReader* fq){
  char result = '\0';
  if (!fq->seqlen || fq->pos >= fq->seqlen ){
    return result;
  }
  result = fq->fastaSeq[fq->pos];
  fq->pos++;
  return result;
}

char getNextFqQual(fastqReader* fq){
  char result = '\0';
  if (!fq->seqlen || fq->pos >= fq->seqlen ){
    return result;
  }
  result = fq->quals[fq->pos];
  fq->pos++;
  return result;
}

