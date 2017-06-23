#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerRead.h"
#include "myseqfunctions/kmerIO.h"
#include <inttypes.h>
#include <unistd.h>


int main(int argc,char** argv){
  if (argc < 2){
    fprintf(stderr, "Use: khSummary k11_file\n");
    return 1;
  }
  D_(0, "Reading %s\n", argv[1]);
  kmerHolder* kh = readIn(argv[1]);
  D_(0, "Done reading kmerHolder\n");
  D_(0, "kmerSize: %" PRIu8 "\n", kh->kmerSize);
  D_(0, "Number of bases: %" PRIu8 "\n", kh->nBases);
  D_(0, "Traces per kmer:\n");
  uint32_t nT = nTraces(&kh, true);
  D_(0, "Number of traces: %" PRIu32 "\n", nT);
  destroyKh(&kh);
}
