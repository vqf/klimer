#ifndef MYDMACROS_H_INCLUDED
#define MYDMACROS_H_INCLUDED

#ifndef LUI
#define LUI long unsigned int
#endif /* BOOL */

#ifndef DEBUG
#define DEBUG 0
#endif /* DEBUG */


#ifndef BOOL
#define BOOL
typedef enum { false, true } bool;
#endif /* BOOL */


#define DOUT stdout
#define D_(a, ...) if (DEBUG >= a) { fprintf(DOUT, "=D_%d: ", __LINE__); fprintf(DOUT, __VA_ARGS__);}
#define E_(a, b) if (DEBUG >= a){ printf("At %d: ", __LINE__); b; }
#define P_(a) printf("%s at %d: %p\n", #a, __LINE__, a);
#define X_ printf("\n=====\nExit at %d\n====\n", __LINE__); exit(0);
#define DIE(...) printf(__VA_ARGS__); printf("File %s, line %d\n", __FILE__, __LINE__); exit(0);

#define SCALAR(a, b) _CTR = 1; \
                  while (a->next){\
                    _CTR++; \
                    a = a->next;\
                  } \
                  b = _CTR;

#define SETKC(a, b) (((a)->flags = (a)->flags | b))
#define UNSETKC(a, b) ((a)->flags = (a)->flags ^ b)
#define ISKC(a, b) ((a)->flags & b)
#define GETUID _UID; _UID++;
#define STOP(...) printf(__VA_ARGS__); getchar();
#endif // MYDMACROS_H_INCLUDED
