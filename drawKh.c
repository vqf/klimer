#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerRead.h"
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/khDraw.h"


int main(int argc,char** argv){
  char* outfile = "/users/vqf/desktop/delme.svg";
  if (argc < 2){
    fprintf(stderr, "Use: khSummary k11_file\n");
    return 1;
  }
  kmerHolder* kh = readIn(argv[1]);
  for (int i = 2; i < argc; i++){
    char* rg = argv[i];
    D_(0, "%s\n", rg);
    if (rg[0] == 0x2d){
      if (EQ(rg, "-v")) DEBUG = 1;
      if (EQ(rg, "-vv")) DEBUG = 2;
    }
    else{
      outfile = rg;
    }
  }
  svgKh* fig = newSvgKh(&kh, 100, 500, 500);
  char* code = getSvg(&fig);
  printf("Writing...\n");
  FILE* fout = fopen(outfile, "w");
  fwrite(code, sizeof(char), strlen(code), fout);
  fclose(fout);
  printf("Done\n");
  free(code);
  destroySvgKh(&fig);
  destroyKh(&kh);
  return 0;
}
