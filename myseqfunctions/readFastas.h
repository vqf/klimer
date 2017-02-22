#include <sys/types.h>
#include "mydmacros.h"

#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif /* _LARGEFILE64_SOURCE */

#define BSIZE 0xFFFF




typedef struct fr{
  FILE* fp;
  char* buffer;
  char* cname;
  uint32_t bufferSize;
  uint32_t bpos; //Position in buffer
  uint32_t gpos; //Position in genome
  off_t cchrpos;
  bool goon;
  bool newchr;
} fastaReader;

bool _readLine(fastaReader*);
char _getNextBase(fastaReader*);

void noAction(memstruct* ms, uint32_t a){
  return;
}

bool looksLikeFasta(char* fname){
  bool result = false;
  FILE* fp = fopen(fname, "r");

  fclose(fp);
}

fastaReader* newFastaReader(FILE* fp, uint32_t bSize){
  if (bSize == 0) bSize = BSIZE;
  fastaReader* result = (fastaReader*) calloc(1, sizeof(fastaReader));
  result->fp = fp;
  result->bufferSize = bSize;
  result->buffer = NULL;
  result->cname = NULL;
  result->bpos = 0;
  result->gpos = 0;
  result->goon = true;
  result->newchr = true;
  return result;
}

void destroyFastaReader(fastaReader** frp){
  D_(2, "Destroying fasta reader\n");
  fastaReader* fr = *frp;
  if (fr->buffer) free(fr->buffer);
  if (fr->cname) free(fr->cname);
  free(fr);
  frp = NULL;
}

void printchars(char *seq){
  int i = 0;
  char t = seq[i];
  while (t != 0x00){
    printf("%02hhx ", t);
    ++i;
    t = seq[i];
  }
  printf("\n");
}


void _getCname(fastaReader* fr){
  /*
  Trims the first character and any symbol with chr lower than 0x10
  */
  char* start = fr->buffer;
  ++start;
  int l = strlen(start);
  int i = 0; int j = 0;

  char* dest = calloc(l+1, sizeof(char)); //I am trimming the first character
  char b = start[i];
  while (i < l && (
         (b >= 0x41 && b <= 0x7a) || //\w
         (b >= 0x30 && b <= 0x39) || //\d
         (b == 0x2d || b == 0x5f)    //_-
         )){
    dest[j] = b;
    ++j;
    ++i;
    b = start[i];
  }
  if (fr->cname) free(fr->cname);
  fr->cname = dest;
  D_(2, "New chromosome %s\n", fr->cname);
}

bool _readLine(fastaReader* fr){
  D_(2, "Reading new file line\n");
  fr->bpos = 0;
  bool result = true;
  if (!fr->buffer && fr->goon){
    D_(2, "Creating %lu - byte buffer\n", (LUI) fr->bufferSize);
    fr->buffer = (char*) calloc(fr->bufferSize, sizeof(char));
  }
  if (fgets(fr->buffer, fr->bufferSize, fr->fp) == NULL){ //Finished file
    free(fr->buffer);
    fr->buffer = NULL;
    fr->goon = false;
    fr->newchr = false;
    fr->gpos = 0;
    result = false;
  }
  else{
    char fst = fr->buffer[fr->bpos];
    D_(2, "%c, pos %d\n", fr->buffer, fr->bpos);
    if (fst == 0x3e){ //>
      //if (fr->cname) free(fr->cname);
      _getCname(fr);
      D_(2, "Chr %s at %li\n", fr->cname, (LUI) ftell (fr->fp));
      fr->newchr = true;
      fr->cchrpos = ftell(fr->fp);
      //fr->goon = false;
      _readLine(fr);
    }
  }
  return result;
}


char _getNextBase(fastaReader* fr){
  char result = '\0';
  if (fr->goon){
    if (!fr->buffer){
      if (!_readLine(fr)){
        fr->goon = false;
        D_(2, "No more to read\n");
        return result;
      }
    }
    bool cont = true;
    if (fr->bpos < strlen(fr->buffer)){
      result = fr->buffer[fr->bpos];
      D_(3, "Giving pos %lu from buffer\n", fr->bpos);
      fr->bpos++;
    }
    else if (_readLine(fr)){
      D_(2, "Had to read another line\n");
      result = _getNextBase(fr);
      cont = false;
    }
    else{
      D_(2, "No more to read\n");
      return result;
    }
    if (result >= 0x41 && result <= 0x7a){ //\w
      if (cont) fr->gpos++;
      //printf("--%d\t%d\t%c\n", (int) fr->gpos, (int) fr->bpos, result);
    }
    else{
      result = _getNextBase(fr);
    }
  }
  return result;
}

bool newChr(fastaReader* fr){
  return fr->newchr;
}

char* chromName(fastaReader* fr){
  return fr->cname;
}

char getNextBase(fastaReader* fr){
  if (fr->newchr){
    fr->newchr = false;
  }
  char r = _getNextBase(fr);
  return r;
}

void resetTo(fastaReader* fr, off_t fpos, char* chrname, uint32_t pos){
  fseek(fr->fp, fpos, SEEK_SET);
  fr->gpos = pos;
  fr->goon = true;
}

off_t gotoChr(fastaReader* fr, char* chrname){
  rewind(fr->fp);
  bool success = false;
  fr->goon = true;
  while (!success && _readLine(fr)){
    if (fr->newchr && *fr->cname == *chrname){
      success = true;
    }
  }
  if (!success) fprintf(stderr, "Could not find chr %s\n", chrname);
  off_t result = ftell(fr->fp);
  resetTo(fr, result, chrname, 0);
  return result;
}




