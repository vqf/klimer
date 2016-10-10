#ifndef DEBUGFILEIO_H_INCLUDED
#define DEBUGFILEIO_H_INCLUDED

uint32_t DBGMEREAD = 0;
uint32_t DBGMEWRIT = 0;

#define fread(a, b, c, d) \
  debfread(a, b, c, d) && ((DBGMEREAD = DBGMEREAD + b * c) && printf("Read %d bytes, %s and %s (%.8x)\n", DBGMEREAD, #b, #c, *(a)))
#define debfread fread

#define fwrite(a, b, c, d) \
  debfwrite(a, b, c, d) && ((DBGMEWRIT = DBGMEWRIT + b * c) && printf("Written %d bytes, %s and %s (%.8x)\n", DBGMEWRIT, #b, #c, *(a)))
#define debfwrite fwrite


#endif // DEBUGFILEIO_H_INCLUDED
