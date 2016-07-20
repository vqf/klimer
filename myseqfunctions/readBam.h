#include <bam.h>

struct GS{
  char* chr;
  uint32_t start;
  uint32_t end;
};

struct BR{
  bamFile in;
  bam_index_t *ind;
  bam1_t *b;
  bam_header_t *header
};

typedef GS gSegment;
typedef BR bamReader;

bamReader* openBAM(char* bamName, char* bamIndex){
  bamReader *result;
  result->in  = bam_open(bamName, "r");
  result->b   = bam_init1();
  result->ind = bam_index_load(bamIndex);
  result->header = bam_header_read(result->in);
  return result;
}