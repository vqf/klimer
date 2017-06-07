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

int _CTR = 1; // Counter for SCALAR
int _T = 0;   // Used to trace how many times a line is passed
int _L = 0;

#define DOUT stdout
#define D_(a, ...) if (DEBUG >= a) { fprintf(DOUT, "=D_%d: ", __LINE__); fprintf(DOUT, __VA_ARGS__);}
#define E_(a, b) if (DEBUG >= a){ printf("At %d: ", __LINE__); b; }
#define P_(a) printf("%s at %d: %p\n", #a, __LINE__, a);
#define X_ printf("\n=====\nExit at %d\n====\n", __LINE__); exit(0);
#define DIE(...) printf(__VA_ARGS__); printf("File %s, line %d\n", __FILE__, __LINE__); exit(0);
#define T_ _T++; if (_L && _L != __LINE__) printf("The _T macro should be used only once\n"); _L = __LINE__;
#define TRACE printf("Line %d used %d times\n", _L, _T);

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
