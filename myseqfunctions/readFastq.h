#include <sys/types.h>

#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif /* _LARGEFILE64_SOURCE */

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
  char* fileName;
  char* seqName;
  char* fastaSeq;
  char* comments;
  char* quals;
  instream type;
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
      fclose(fip);
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

instream ftype(fastqReader* fq){
  return fq->type;
}

bool getNextFqRead(fastqReader** fqp){
  fastqReader* fq = *fqp;
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
      _resetFq(&fq);
      result = true;
    }
  }
  else if ( fgets(fq->seqName,  NAMESIZE, fq->fp) &&
            fgets(fq->fastaSeq,  MAXREAD, fq->fp) &&
            fgets(fq->comments, NAMESIZE, fq->fp) &&
            fgets(fq->quals,     MAXREAD, fq->fp)
            ){
    _resetFq(&fq);
    result = true;
  }
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

char getNextFqQual(fastqReader** fqp){
  fastqReader* fq = *fqp;
  char result = '\0';
  if (!fq->seqlen || fq->pos >= fq->seqlen ){
    return result;
  }
  result = fq->quals[fq->pos];
  fq->pos++;
  return result;
}

void destroyFastqReader(fastqReader** fqrp){
  fastqReader* fqr = *fqrp;
  free(fqr->comments);
  free(fqr->fastaSeq);
  free(fqr->quals);
  free(fqr->seqName);
  if (fqr->type == normal){
    fclose(fqr->fp);
  }
  else if (fqr->type == gzipped){
    gzclose(*(fqr->zfp));
  }
  free(fqr);
}
