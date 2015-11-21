#ifndef MOCK_SEQREADER_HPP
#define MOCK_SEQREADER_HPP

#include "seqreader.hpp"
#include "rawseqs.h"
#include "kraken_headers.hpp"

namespace kraken {

  class MockSequenceReader
  {
     public:
     MockSequenceReader(int num_threads);
     DNASequence next_sequence(int thread_num);
     bool is_valid(int thread_num);
    
     private:
     bool valid;
     const static int num_reads=2000;
     const static int reads_per_thread = 200010;
     int total_threads;
     int total_reads;
     int* thread_index;
     int* read_index;
  };
}

#endif
