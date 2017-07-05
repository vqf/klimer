#ifndef KHDRAW_H_INCLUDED
#define KHDRAW_H_INCLUDED

#include <inttypes.h>
#include "kmer.h"
#include <math.h>

#define PI 3.14159265

typedef struct{
  float x;
  float y;
} point;

typedef struct lines{
  char* line;
  struct lines* next;
} lines;


typedef struct svgKh{
  uint32_t R0;
  uint32_t w;
  uint32_t h;
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

void pushLine(lines** lp, char* txt){
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

void adjust(char** s){
  char* tmp = *s;
  char* result = (char*) calloc(strlen(tmp) + 1, sizeof(char));
  memcpy(result, tmp, strlen(tmp));
  free(tmp);
  *s = result;
}

char* buildCirc(uint32_t cx, uint32_t cy, uint32_t R){
  char* result = (char*) calloc(100, sizeof(char));
  sprintf(result, "<circle cx=\"%"PRIu32"\" cy=\"%"PRIu32"\" r=\"%"PRIu32"\" class=\"tdisk\" stroke=\"black\" fill=\"none\"/>\n",
          cx, cy, R);
  adjust(&result);
  return result;
}

char* buildSegment(svgKh** sp, point p1, point p2, point p3, point p4, char* color){
  svgKh* s = *sp;
  char* result = (char*) calloc(500, sizeof(char));
  char* tcolor = "black";
  float cx = (float) (s->w/2);
  float cy = (float) (s->h/2);
  float strw = 0.1;
  point q = {(p1.x + p2.x) / 2, (p1.y + p2.y) / 2};
  point rf = {(q.x + cx) / 2, (q.y + cy) / 2};
  sprintf(result, "<path class=\"tick\" d=\"M %.2f %.2f L %.2f %.2f\" \
                   stroke-width=\"%.2f\" stroke=\"%s\" />\n\
                   <path class=\"myarc\" d=\"M %.2f %.2f Q %.2f %.2f %.2f %.2f\" \
                   stroke-width=\"0.2\" fill=\"none\" stroke=\"%s\" />\n",
          p1.x, p1.y, p2.x, p2.y, strw, tcolor,
          p2.x, p2.y, rf.x, rf.y, p3.x, p3.y, color
          );
  adjust(&result);
  return result;
}

char* setColor(uint32_t v){
  uint32_t top = 2;
  if (v > top) v = top;
  uint8_t red = 255 * v / top;
  uint8_t blue = 255 - red;
  char* result = (char*) calloc(500, sizeof(char));
  sprintf(result, "rgb(%"PRIu8",0,%"PRIu8")", red, blue);
  adjust(&result);
  return result;
}

void populateSvg(svgKh** sp){
  svgKh* s = *sp;
  lines* l = s->code;
  float r = (float) s->R0;
  float cx = (float) s->w/2;
  float cy = (float) s->h/2;
  pushLine(&l, buildCirc(cx, cy, s->R0));
  memstruct* ms = s->kh->ms;
  uint32_t np = ms->nPos;
  float cang = (float) (2 * PI / np);
  float tick = 5;
  for (uint32_t i = 0; i < np; i++){
    kmerConnector* kc = ms->kmerArray[i];
    while (kc){
      point s1 = {(float) (cx + (r + tick/2) * cos(i * cang)),
                 (float) (cy + (r + tick/2) * sin(i * cang))};
      point t1 = {(float) (cx + (r - tick/2) * cos(i * cang)),
                 (float) (cy + (r - tick/2) * sin(i * cang))};
      uint32_t j = kc->dest;
      uint32_t nf = 0;
      SCALAR(kc->idflags, nf);
      char* scolor = setColor(nf);
      point s2 = {(float) (cx + (r + tick/2) * cos(j * cang)),
                 (float) (cy + (r + tick/2) * sin(j * cang))};
      point t2 = {(float) (cx + (r - tick/2) * cos(j * cang)),
                 (float) (cy + (r - tick/2) * sin(j * cang))};
      pushLine(&l, buildSegment(sp, s1, t1, t2, s2, scolor));
      free(scolor);
      kc = kc->next;
    }
  }
}

svgKh* newSvgKh(kmerHolder** khp, uint32_t R0,
                uint32_t width, uint32_t height){
  svgKh* result = (svgKh*) calloc(1, sizeof(svgKh));
  result->kh = *khp;
  result->R0 = R0;
  result->w = width;
  result->h = height;
  char* header = (char*) calloc(70, sizeof(char));
  sprintf(header, "<svg width=\"%"PRIu32"px\" height=\"%"PRIu32"px\">\n", width, height);
  result->header = header;
  result->code = newLines(header);
  pushLine(&result->code, "<desc>kmerHold depiction</desc>\n");
  populateSvg(&result);
  pushLine(&result->code, "</svg>\n");
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
