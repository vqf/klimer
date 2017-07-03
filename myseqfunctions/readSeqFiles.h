#ifndef READSEQFILES_H_INCLUDED
#define READSEQFILES_H_INCLUDED

#include "stdio.h"
#include "readFastq.h"
#include "readFastas.h"


typedef enum { fasta, fastq } seqFileType;

typedef struct seqReader{
  fastaReader* far;
  fastqReader* fqr;
  seqFileType ftype;
} seqReader;

seqReader* newSeqReader(char* infile){
  seqReader* result = (seqReader*) calloc(1, sizeof(seqReader));
  result->far = NULL;
  result->fqr = NULL;
  if (looksLikeFasta(infile)){
    result->far = initFastaReader(infile, 0);
    result->ftype = fasta;
    D_(1, "Looks like Fasta\n");
  }
  else if (looksLikeFastq(infile)){
    result->fqr = initFastqReader(infile);
    result->ftype = fastq;
    D_(1, "Looks like Fastq\n");
  }
  else{
    DIE("Cannot recognize format\n");
  }
  return result;
}

//#define newSeqReader(a) newSeqReader(a); printf("File %s, line %d\n", __FILE__, __LINE__);

void destroySeqReader(seqReader** sfp){
  seqReader* sf = *sfp;
  destroyFastaReader(&sf->far);
  destroyFastqReader(&sf->fqr);
  free(sf);
  sfp = NULL;
}

char getNextBase(seqReader** fqp){
  seqReader* fq = *fqp;
  char result = '\0';
  if (fq->ftype == fasta){
    result = _getNextBase(fq->far);
  }
  else if (fq->ftype == fastq){
    if (!fq->fqr->seqlen || fq->fqr->pos >= fq->fqr->seqlen ){
      return result;
    }
    result = fq->fqr->fastaSeq[fq->fqr->pos];
    fq->fqr->pos++;
  }
  return result;
}

bool getNextRead(seqReader** fqp){
  seqReader* fq = *fqp;
  bool result = false;
  D_(1, "Ftype: %d\n", fq->ftype);
  if (fq->ftype == fasta){
    if (fq->far->newchr){
      fq->far->newchr = false;
      result = true;
    }
  }
  else if (fq->ftype == fastq){
    result = getNextFqRead(&fq->fqr);
  }
  return result;
}

char* chromName(seqReader** srp){
  seqReader* sr = *srp;
  char* result = NULL;
  if (sr->ftype == fasta){
    result = sr->far->cname;
  }
  else if (sr->ftype == fastq){
    result = sr->fqr->seqName;
  }
  return result;
}


#endif // READSEQFILES_H_INCLUDED
