#ifndef KHDRAW_H_INCLUDED
#define KHDRAW_H_INCLUDED

#include <inttypes.h>
#include "kmer.h"

typedef struct lines{
  char* line;
  struct lines* next;
} lines;


typedef struct svgKh{
  uint32_t R0;
  kmerHolder* kh;
  char* header;
  lines* code;
  char* closer;
} svgKh;


lines* newLines(char* txt){
  lines* result = (lines*) calloc(1, sizeof(lines));
  result->line = txt;
  result->next = NULL;
  return result;
}

void pushString(lines** lp, char* txt){
  lines* l = *lp;
  if (!l){
    lines* n = newLines(txt);
    *lp = n;
    return;
  }
  while (l->next){
    l = l->next;
  }
  l->next = newLines(txt);
}

char* linesToStr(lines** lp){
  lines* l = *lp;
  uint32_t sl = 0;
  while (l){
    sl += strlen(l->line);
    l = l->next;
  }
  char* result = (char*) calloc(sl + 1, sizeof(char));
  l = *lp;
  uint32_t c = 0;
  while (l){
    char* ptr = &result[c];
    memcpy(ptr, l->line, strlen(l->line));
    c += strlen(l->line);
    l = l->next;
  }
  return result;
}

void destroyLines(lines** lp){
  if (!lp) return;
  lines* l = *lp;
  if (!l) return;
  while (l){
    lines* nxt = l->next;
    free(l);
    l = nxt;
  }
}

svgKh* newSvgKh(kmerHolder** khp, uint32_t R0,
                uint32_t width, uint32_t height){
  svgKh* result = (svgKh*) calloc(1, sizeof(svgKh));
  result->kh = *khp;
  result->R0 = R0;
  char* header = (char*) calloc(70, sizeof(char));
  sprintf(header, "<svg width=\"%"PRIu32"px\" height=\"%"PRIu32"px\">\n", width, height);
  result->header = header;
  result->code = newLines(header);
  pushString(&result->code, "<desc>kmerHold depiction</desc>\n");
  pushString(&result->code, "</svg>\n");
  result->closer = "</svg>\n";
  return result;
}

char* getSvg(svgKh** sp){
  svgKh* s = *sp;

  char* result = linesToStr(&s->code);
  return result;
}

void destroySvgKh(svgKh** sp){
  svgKh* s = *sp;
  free(s->header);
  destroyLines(&s->code);
  free(s);
  *sp = NULL;
}

#endif // KHDRAW_H_INCLUDED
