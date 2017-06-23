#ifndef IOLAYER_H_INCLUDED
#define IOLAYER_H_INCLUDED
#include <sys/types.h>
#include <stdio.h>
#include "mydmacros.h"
#include <zlib.h>
#include <errno.h>

#ifndef OS
#define OS 0
#ifdef _WIN32
#undef OS
#define OS 1
#endif // _WIN32
#endif // OS

#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif /* _LARGEFILE64_SOURCE */

#define BSIZE 0xFFFF

typedef enum { normal, gzipped, zipped } instream;

typedef struct fpointer{
  FILE* fp;
  gzFile zfp;
  instream type;
} fpointer;

fpointer* newFPointer(char* fName){
  fpointer* result = (fpointer*) calloc(1, sizeof(fpointer));
  if (strcmp(fName, "-") == 0){
    result->fp = stdin;
    D_(1, "Reading from stdin\n");
  }
  else{
    FILE* fip = fopen(fName, "r");
    if (!fip){
      D_(0, "Could not open %s\n", fName);
      return NULL;
    }
    uint8_t magic[2];
    fread(magic, 2, 1, fip);
    //fprintf(stderr, "Magic: %x, %x\n", magic[0], (uint8_t) magic[1]);
    if (magic[0] == 0x1f && magic[1] == 0x8b){  //gzipped
      result->type = gzipped;
      D_(1, "Gzipped\n");
      fclose(fip);
      fip = NULL;
    }
    if (OS == 1){ //Windows
      if (magic[0] == 0x50 && magic[1] == 0x4b){  //zipped
        result->type = zipped;
        fclose(fip);
        fip = NULL;
        D_(0, "zipped, not supported yet\n");
      }
    }
    if (fip){
      fseek(fip, 0, 0);
      result->fp = fip;
      D_(1, "Normal text\n");
    }
    else if (result->type == gzipped){
      gzFile zfip = gzopen(fName, "r");
      if (! zfip) {
        fprintf (stderr, "gzopen of '%s' failed: %s.\n", fName,
                 strerror (errno));
            exit (EXIT_FAILURE);
      }
      result->zfp = zfip;
      //fprintf(stderr, "Gzipped\n");
    }
    else{
      D_(0, "Could not open %s\n", fName);
    }
  }
  return result;
}
//#define newFPointer(a) newFPointer(a); printf("File %s, line %d\n", __FILE__, __LINE__);

void destroyFPointer(fpointer** frp){
  if (frp && *frp){
    fpointer* fr = *frp;
    if (fr->fp) fclose(fr->fp);
    if (fr->zfp) gzclose(fr->zfp);
    free(fr);
    frp = NULL;
  }
}


char* myfgets(fpointer** frp, char** dest, int bufferSize){
  fpointer* fr = *frp;
  char* d = *dest;
  char* result = NULL;
  if (bufferSize == 0){
    bufferSize = BSIZE;
  }
  if (fr->type == gzipped){
    result = gzgets(fr->zfp, d, bufferSize);
  }
  else{
    result = fgets(d, bufferSize, fr->fp);
    D_(2, "%s\n", d);
  }
  return result;
}

instream ftype(fpointer* fq){
  return fq->type;
}




#endif // IOLAYER_H_INCLUDED
