#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerRead.h"
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/khDraw.h"


int main(int argc,char** argv){
  if (argc < 2){
    fprintf(stderr, "Use: khSummary k11_file\n");
    return 1;
  }
  kmerHolder* kh = readIn(argv[1]);
  svgKh* fig = newSvgKh(&kh, 100, 500, 500);
  char* code = getSvg(&fig);
  printf("%s\n", code);
  free(code);
  destroySvgKh(&fig);
  destroyKh(&kh);
  return 0;
}
