#include "mockseqreader.hpp"
#include "kraken_headers.hpp"
#include "seqreader.hpp"

using namespace std;
namespace kraken {
  MockSequenceReader::MockSequenceReader(int num_threads)
  {
    total_threads = num_threads;
    total_reads = total_threads * reads_per_thread; 
    thread_index = (int*) calloc(num_threads,sizeof(int));
    read_index = (int*) calloc(num_reads,sizeof(int));
    for(int i=0;i<num_threads;i++)
    { 
      thread_index[i]=0;
      if(i<num_reads)
        read_index[i]=0;
    }   
  }
  
  bool MockSequenceReader::is_valid(int thread_num) {
      if(thread_index[thread_num] < reads_per_thread)
         return true;
      return false;
  }
  
  DNASequence MockSequenceReader::next_sequence(int thread_num)
  {
    DNASequence dna;
    if(thread_index[thread_num] < reads_per_thread)
    {
      if(read_index[thread_num] >= num_reads)
          read_index[thread_num]=0;
      dna.header_line = raw_list[read_index[thread_num]].id;
      istringstream seq_id(dna.header_line);
      seq_id >> dna.id;
      dna.seq = raw_list[read_index[thread_num]].seq;
      dna.quals = raw_list[read_index[thread_num]].qual;
      read_index[thread_num]++;
      thread_index[thread_num]++;
    }
    return dna;
  }
}
