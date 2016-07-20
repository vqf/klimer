#define BSIZE 0xFFFF

char* cname = NULL; // Current chr name
uint32_t chrpos = 0;



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

bool readLine(fastaReader*);

void noAction(memstruct* ms, uint32_t a){
  return;
}

fastaReader* initFastaReader(FILE* fp, uint32_t bSize){
  if (bSize == 0) bSize = BSIZE;
  fastaReader* result = malloc(sizeof(fastaReader));
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


void getCname(fastaReader* fr){
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
}

bool readLine(fastaReader* fr){
  fr->bpos = 0;
  bool result = true;
  if (!fr->buffer && fr->goon){
    fr->buffer = (char*) calloc(fr->bufferSize, sizeof(char));
  }
  if (fgets(fr->buffer, fr->bufferSize, fr->fp) == NULL){ //Finished file
    free(fr->buffer);
    fr->buffer = NULL;
    fr->goon = false;
    fr->gpos = 0;
    result = false;
  }
  else{
    char fst = fr->buffer[fr->bpos];
    if (DEBUG) printf("%s, pos %d\n", fr->buffer, fr->bpos);
    if (fst == 0x3e){ //>
      //if (fr->cname) free(fr->cname);
      getCname(fr);
      if (DEBUG) fprintf(stderr, "Chr %s at %li\n", fr->cname, (long int) ftello (fr->fp));
      fr->newchr = true;
      fr->cchrpos = ftello(fr->fp);
      readLine(fr);
    }
  }
  return result;
}


char _getNextBase(fastaReader* fr){
  char result = 0x00;
  if (fr->goon){
    if (!fr->buffer){
      if (!readLine(fr)){
        fr->goon = false;
        return result;
      }
    }
    bool count = true;
    if (fr->bpos < strlen(fr->buffer)){
      result = fr->buffer[fr->bpos];
      fr->bpos++;
    }
    else if (readLine(fr)){
      result = _getNextBase(fr);
      count = false;
    }
    else{
      return result;
    }
    if (result >= 0x41 && result <= 0x7a){ //\w
      if (count) fr->gpos++;
      //printf("--%d\t%d\t%c\n", (int) fr->gpos, (int) fr->bpos, result);
    }
    else{
      result = _getNextBase(fr);
    }
  }
  return result;
}

char getNextBase(fastaReader* fr){
  if (fr->newchr == true) fr->newchr = false;
  char r = _getNextBase(fr);
  return r;
}

void resetTo(fastaReader* fr, off_t fpos, char* chrname, uint32_t pos){
  fseeko(fr->fp, fpos, SEEK_SET);
  fr->gpos = pos;
  fr->goon = true;
}

off_t gotoChr(fastaReader* fr, char* chrname){
  rewind(fr->fp);
  bool success = false;
  fr->goon = true;
  while (!success && readLine(fr)){
    if (fr->newchr && *fr->cname == *chrname){
      success = true;
    }
  }
  if (!success) fprintf(stderr, "Could not find chr %s\n", chrname);
  off_t result = ftello(fr->fp);
  resetTo(fr, result, chrname, 0);
  return result;
}




