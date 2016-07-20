#ifndef NEATHTML_H_INCLUDED
#define NEATHTML_H_INCLUDED



#endif // NEATHTML_H_INCLUDED

typedef struct tg{
  char* tagName;
  char* id;
  char* clss;
  char* cx;
  char* cy;
  char* text;
  struct tg* content;
  struct tg* next;
} tag;



tag* newTag(char* tname){
  tag* result = (tag*) calloc(1, sizeof(tag));
  result->tagName = tname;
  result->next    = NULL;
  result->id      = NULL;
  result->clss    = NULL;
  result->content = NULL;
  result->text    = NULL;
  return result;
}
void addContent(tag** tp, tag* s){
  tag* t = *tp;
  t->content = s;
}
void addTag(tag** tp, tag* s){
  tag* t = *tp;
  t->next = s;
}

void addText(tag** tp, char* txt){
  tag* t = *tp;
  if (strcmp(t->tagName, "text") != 0){
    fprintf(stderr, "Warning: Tried to add text to %s\n", t->tagName);
  }
  size_t l = strlen(txt);
  t->text = (char*) calloc(l + 1, sizeof(char));
  strcpy(t->text, txt);
}

void destroyTag(tag** tp){
  tag* t = *tp;
  while (t->next){
    destroyTag(&t->next);
  }
  while (t->content){
    destroyTag(&t->content);
  }
  free(t->text);
  free(*tp);
  *tp = NULL;
}

void concat(char** c1p, char* c2){
  char* result = NULL;
  char* c1 = *c1p;
  if (c1 && c2){
    size_t l1 = strlen(c1);
    size_t l2 = strlen(c2);
    result = (char*) calloc(l1 + l2 + 1, sizeof(char));
    if (result){
      memcpy(result, c1, l1);
      memcpy(result + l1, c2, l2);
      free(c1);
      *c1p = result;
    }
  }
}

char* sprintTag(tag* t){
  if (!t){
    return "";
  }
  char* result = (char*) calloc(1, sizeof(char));
  do{
    if (strcmp(t->tagName, "text") == 0){
      concat(&result, t->text);
      return result;
    }
    else{
      char* ttag = (char*) calloc(1, sizeof(char));
      concat(&ttag, "<");
      concat(&ttag, t->tagName);
      concat(&ttag, ">\n");
      tag* ctag = t->content;
      while(ctag){
        char* dcont = sprintTag(ctag);
        concat(&ttag, dcont);
        ctag = ctag->content;
      }
      concat(&ttag, "\n</");
      concat(&ttag, t->tagName);
      concat(&ttag, ">\n");
      concat(&result, ttag);
      free(ttag);
    }
    t = t->next;
  } while (t);
  return result;
}

