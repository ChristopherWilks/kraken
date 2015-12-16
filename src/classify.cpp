/*
 * Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken taxonomic sequence classification system.
 *
 * Kraken is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kraken is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <tbb/tbb.h>
#include <tbb/task_group.h>
#include <tbb/mutex.h>
#include <tbb/spin_mutex.h>
#include <tbb/task_group.h>
#include <tbb/atomic.h>
#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "quickfile.hpp"
#include "seqreader.hpp"
#include "ThreadProfile.h"
#include "tbb_classify.h"

#define MUTEX_T tbb::mutex
//#define MUTEX_T tbb::spin_mutex

const size_t DEF_WORK_UNIT_SIZE = 500000;

using namespace std;
using namespace kraken;

void parse_command_line(int argc, char **argv);
void usage(int exit_code=EX_USAGE);
void process_file(char *filename,int numOfThreads);
//void* pclassify(void* args_);
void classify_sequence(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss);
string hitlist_string(vector<uint32_t> &taxa, vector<uint8_t> &ambig);
set<uint32_t> get_ancestry(uint32_t taxon);
void report_stats(struct timeval time1, struct timeval time2);

int Num_threads = 1;
string DB_filename, Index_filename, Nodes_filename;
bool Quick_mode = false;
bool Fastq_input = false;
bool Print_classified = false;
bool Print_unclassified = false;
bool Print_kraken = true;
bool Populate_memory = false;
bool Only_classified_kraken_output = false;
uint32_t Minimum_hit_count = 1;
map<uint32_t, uint32_t> Parent_map;
KrakenDB Database;
string Classified_output_file, Unclassified_output_file, Kraken_output_file;
ostream *Classified_output;
ostream *Unclassified_output;
ostream *Kraken_output;
size_t Work_unit_size = DEF_WORK_UNIT_SIZE;

uint64_t total_classified = 0;
tbb::atomic<uint64_t> total_classified_atomic;

uint64_t total_sequences = 0;
uint64_t total_bases = 0;
	
int main(int argc, char **argv) {
  /*#ifdef _OPENMP
  omp_set_num_threads(1);
  #endif*/
  total_classified_atomic=0;

  parse_command_line(argc, argv);
  if (! Nodes_filename.empty())
    Parent_map = build_parent_map(Nodes_filename);

  if (Populate_memory)
    cerr << "Loading database... ";

  QuickFile db_file;
  db_file.open_file(DB_filename);
  if (Populate_memory)
    db_file.load_file();
  Database = KrakenDB(db_file.ptr());
  KmerScanner::set_k(Database.get_k());

  QuickFile idx_file;
  idx_file.open_file(Index_filename);
  if (Populate_memory)
    idx_file.load_file();
  KrakenDBIndex db_index(idx_file.ptr());
  Database.set_index(&db_index);

  if (Populate_memory)
    cerr << "complete." << endl;

  if (Print_classified) {
    if (Classified_output_file == "-")
      Classified_output = &cout;
    else
      Classified_output = new ofstream(Classified_output_file.c_str());
  }

  if (Print_unclassified) {
    if (Unclassified_output_file == "-")
      Unclassified_output = &cout;
    else
      Unclassified_output = new ofstream(Unclassified_output_file.c_str());
  }

  if (! Kraken_output_file.empty()) {
    if (Kraken_output_file == "-")
      Print_kraken = false;
    else
      Kraken_output = new ofstream(Kraken_output_file.c_str());
  }
  else
    Kraken_output = &cout;

  struct timeval tv1, tv2;
  gettimeofday(&tv1, NULL);
  for (int i = optind; i < argc; i++)
    process_file(argv[i],Num_threads);
  gettimeofday(&tv2, NULL);

  total_classified = total_classified_atomic; 
  report_stats(tv1, tv2);

  return 0;
}

void report_stats(struct timeval time1, struct timeval time2) {
  time2.tv_usec -= time1.tv_usec;
  time2.tv_sec -= time1.tv_sec;
  if (time2.tv_usec < 0) {
    time2.tv_sec--;
    time2.tv_usec += 1000000;
  }
  double seconds = time2.tv_usec;
  seconds /= 1e6;
  seconds += time2.tv_sec;

  cerr << "\r";
  fprintf(stderr, 
          "%llu sequences (%.2f Mbp) processed in %.3fs (%.1f Kseq/m, %.2f Mbp/m).\n",
          (unsigned long long) total_sequences, total_bases / 1.0e6, seconds,
          total_sequences / 1.0e3 / (seconds / 60),
          total_bases / 1.0e6 / (seconds / 60) );
  fprintf(stderr, "  %llu sequences classified (%.2f%%)\n",
          (unsigned long long) total_classified, total_classified * 100.0 / total_sequences);
  fprintf(stderr, "  %llu sequences unclassified (%.2f%%)\n",
          (unsigned long long) (total_sequences - total_classified),
          (total_sequences - total_classified) * 100.0 / total_sequences);
}

struct ReadKmersArg
{
	//pthread_t thread_id;
        int thread_num;
        //DNASequenceReader* reader;
        FastqReader* reader;
	MUTEX_T* readerLock;
	MUTEX_T* writerLock;
	MUTEX_T* classifiedLock;
        int read_length;
	DNASequence** dnaA;
	int num_reads;
	//pthread_mutex_t* readerLock;
};


void clear_dnaseqs(DNASequence** dnaA,int num_reads)
{
     int z = 0;
     for ( z = 0 ; z < num_reads; z++)
     {
	dnaA[z]->seq.clear();	
	dnaA[z]->id.clear();	
	dnaA[z]->header_line.clear();	
	dnaA[z]->quals.clear();
     }
}

void process_file(char *filename,int numOfThreads) {
  string file_str(filename);
  //DNASequenceReader *reader;
  FastqReader* reader;

  if (Fastq_input)
    reader = new FastqReader(file_str);
  //else
  //  reader = new FastaReader(file_str);

  ReadKmersArg* argss = ( ReadKmersArg* ) calloc(numOfThreads, sizeof(struct ReadKmersArg)); 
  
  tbb::task_group tbb_grp;
  MUTEX_T mutexSampleKmers;
  MUTEX_T mutexWriteKmers;
  int i=0;
  int z=0;
  int read_length=100;
  int num_reads = DEF_WORK_UNIT_SIZE/read_length;
  for ( i = 0 ; i < numOfThreads; ++i )
  {
     DNASequence **dnaA=(DNASequence**) malloc(num_reads*sizeof(DNASequence*));
     /*for ( z = 0 ; z < num_reads; z++)
     {
	//printf("%d %d %d\n",z,read_length,num_reads);
        DNASequence d;
        dnaA[z]=&d;
        dnaA[z]->seq.reserve(read_length+1);
        dnaA[z]->quals.reserve(read_length+1);
        dnaA[z]->header_line.reserve(read_length+1);
        dnaA[z]->id.reserve(read_length+1);
     }*/
     //printf("after z\n");
     argss[i].reader=reader;
     argss[i].readerLock=&mutexSampleKmers;
     argss[i].writerLock=&mutexWriteKmers;
     argss[i].thread_num=i+1;
     argss[i].read_length=read_length; //TODO: need to "peak" at the file and determine what the read length should be
     argss[i].num_reads=num_reads; //TODO: need to "peak" at the file and determine what the read length should be
     argss[i].dnaA=dnaA;
     //pthread_create( &(argss[i].thread_id), &pthreadAttr, pclassify, (void *) &argss[i] );
     tbb_grp.run(pclassify(&argss[i]));
  }
  tbb_grp.wait();
  
  delete reader;
}

  //#pragma omp parallel
//#static void* pclassify(void* args_) //DNASequenceReader *reader, void *arg)
void pclassify::operator()() //DNASequenceReader *reader, void *arg)
{
    struct ReadKmersArg* args=(struct ReadKmersArg*) args_;
    char* msg=(char*) calloc(1024,sizeof(char));
    sprintf(msg,"PClassify thread %d",args->thread_num);
    ThreadProfile *tp = new ThreadProfile(msg);


    //DNASequenceReader* reader = args->reader;
    FastqReader* reader = args->reader;
    vector<DNASequence> work_unit;
    ostringstream kraken_output_ss, classified_output_ss, unclassified_output_ss;
    int read_length=101;
    DNASequence dna;
    /*dna.seq.reserve(read_length);
    dna.id.reserve(read_length);
    dna.header_line.reserve(read_length);
    dna.quals.reserve(read_length);*/

    char* seq=(char*) calloc(101,sizeof(char));
    char* id=(char*) calloc(101,sizeof(char));
    char* header_line=(char*) calloc(101,sizeof(char));
    char* quals=(char*) calloc(101,sizeof(char));

    /*dna.seq="AGAAATGGCTTGATGACTAGTAGGAATAAGGGGGAGAAAGTAAGTGAAAATTAAATTGAAGTAAAGAAAAAATGAAAAATAAAATAAAAAAGGAAGGAAG";
    dna.quals="FFDADFG?FGFA5AA8:>>@25555@5=A@DDDD98.<7@ADDDD?9<1A=CCC<GG?=FGDDDI=6=0)08<C6DGF6F9=?@?>?>BB>BE?GGG>GG";
    dna.id="SRR034966.253";
    dna.header_line="SRR034966.253 090421_HWI-EAS255_9111_FC400PT_PE_7_1_10_2018 length=100";*/
   
    //for each thread (here) reserve it's own block of DNA strings to fit up to the total DNA chars
    //"peak" at the length of a DNA string in the single thread portion to set this up if need be
    //we could also reserve one large 500000 sized string, keep 2 extra arrays == num sequences
    //array 1 stores all the headers, array 2 stores offset to start/length of each string in the block
    int rounds=0;
    while (reader->is_valid()) { // and rounds  < 40) {
      rounds++;
      tp->update();
      work_unit.clear();
      
      /*dna.seq.clear();
      dna.header_line.clear();
      dna.id.clear();
      dna.quals.clear();*/
 
      //clear_dnaseqs(args->dnaA,args->num_reads);
      size_t total_nt = 0;
      int a = 0;
      //#pragma omp critical(get_input)
      //pthread_mutex_lock(args->readerLock);
      args->readerLock->lock();
      {
        //while (total_nt < Work_unit_size and a < args->num_reads) {
        while (total_nt < Work_unit_size) {
          reader->next_sequence(&seq,&id,&header_line,&quals);
          //DNASequence dna;
          dna.seq=string(seq); 
          dna.id=string(id); 
          dna.header_line=string(header_line); 
          dna.quals=string(quals); 
          //*(args->dnaA[a])=reader->next_sequence(*(args->dnaA[a]));
          if (! reader->is_valid())
            break;
          work_unit.push_back(dna);
          total_nt += dna.seq.size();
          /*work_unit.push_back(*(args->dnaA[a]));
          total_nt += args->dnaA[a]->seq.size();
          a++;*/
        }
      }
      args->readerLock->unlock();
      //pthread_mutex_unlock(args->readerLock);
      if (total_nt == 0)
        break;
      
      kraken_output_ss.str("");
      classified_output_ss.str("");
      unclassified_output_ss.str("");
      for (size_t j = 0; j < work_unit.size(); j++)
        classify_sequence( work_unit[j], kraken_output_ss,
                           classified_output_ss, unclassified_output_ss );

      //#pragma omp critical(write_output)
      //pthread_mutex_lock(args->readerLock);
      args->writerLock->lock();
      {
        if (Print_kraken)
          (*Kraken_output) << kraken_output_ss.str();
        if (Print_classified)
          (*Classified_output) << classified_output_ss.str();
        if (Print_unclassified)
          (*Unclassified_output) << unclassified_output_ss.str();
        total_sequences += work_unit.size();
        total_bases += total_nt;
        cerr << "\r" << args->thread_num << " Processed " << total_sequences << " sequences (" << total_bases << " bp) ...";
      }
      args->writerLock->unlock();
      //pthread_mutex_unlock(args->readerLock);
    }
    tp->finish(); 
    delete tp;
    free(msg);
  //return NULL;
}  // end parallel section


void classify_sequence(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss) {
  vector<uint32_t> taxa;
  vector<uint8_t> ambig_list;
  map<uint32_t, uint32_t> hit_counts;
  uint64_t *kmer_ptr;
  uint32_t taxon = 0;
  uint32_t hits = 0;  // only maintained if in quick mode

  uint64_t current_bin_key;
  int64_t current_min_pos = 1;
  int64_t current_max_pos = 0;

  if (dna.seq.size() >= Database.get_k()) {
    KmerScanner scanner(dna.seq);
    while ((kmer_ptr = scanner.next_kmer()) != NULL) {
      taxon = 0;
      if (scanner.ambig_kmer()) {
        ambig_list.push_back(1);
      }
      else {
        ambig_list.push_back(0);
        uint32_t *val_ptr = Database.kmer_query(
                              Database.canonical_representation(*kmer_ptr),
                              &current_bin_key,
                              &current_min_pos, &current_max_pos
                            );
        taxon = val_ptr ? *val_ptr : 0;
        if (taxon) {
          hit_counts[taxon]++;
          if (Quick_mode && ++hits >= Minimum_hit_count)
            break;
        }
      }
      taxa.push_back(taxon);
    }
  }

  uint32_t call = 0;
  if (Quick_mode)
    call = hits >= Minimum_hit_count ? taxon : 0;
  else
    call = resolve_tree(hit_counts, Parent_map);

  if (call)
    total_classified_atomic.fetch_and_increment();
    //#pragma omp atomic
    //total_classified++;

  if (Print_unclassified || Print_classified) {
    ostringstream *oss_ptr = call ? &coss : &uoss;
    bool print = call ? Print_classified : Print_unclassified;
    if (print) {
      if (Fastq_input) {
        (*oss_ptr) << "@" << dna.header_line << endl
            << dna.seq << endl
            << "+" << endl
            << dna.quals << endl;
      }
      else {
        (*oss_ptr) << ">" << dna.header_line << endl
            << dna.seq << endl;
      }
    }
  }

  if (! Print_kraken)
    return;

  if (call) {
    koss << "C\t";
  }
  else {
    if (Only_classified_kraken_output)
      return;
    koss << "U\t";
  }
  koss << dna.id << "\t" << call << "\t" << dna.seq.size() << "\t";

  if (Quick_mode) {
    koss << "Q:" << hits;
  }
  else {
    if (taxa.empty())
      koss << "0:0";
    else
      koss << hitlist_string(taxa, ambig_list);
  }

  koss << endl;
}

string hitlist_string(vector<uint32_t> &taxa, vector<uint8_t> &ambig)
{
  int64_t last_code;
  int code_count = 1;
  ostringstream hitlist;

  if (ambig[0])   { last_code = -1; }
  else            { last_code = taxa[0]; }

  for (size_t i = 1; i < taxa.size(); i++) {
    int64_t code;
    if (ambig[i]) { code = -1; }
    else          { code = taxa[i]; }

    if (code == last_code) {
      code_count++;
    }
    else {
      if (last_code >= 0) {
        hitlist << last_code << ":" << code_count << " ";
      }
      else {
        hitlist << "A:" << code_count << " ";
      }
      code_count = 1;
      last_code = code;
    }
  }
  if (last_code >= 0) {
    hitlist << last_code << ":" << code_count;
  }
  else {
    hitlist << "A:" << code_count;
  }
  return hitlist.str();
}

set<uint32_t> get_ancestry(uint32_t taxon) {
  set<uint32_t> path;

  while (taxon > 0) {
    path.insert(taxon);
    taxon = Parent_map[taxon];
  }
  return path;
}

void parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "d:i:t:u:n:m:o:qfcC:U:M")) != -1) {
    switch (opt) {
      case 'd' :
        DB_filename = optarg;
        break;
      case 'i' :
        Index_filename = optarg;
        break;
      case 't' :
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive thread count");
        /*#ifdef _OPENMP
        if (sig > omp_get_num_procs())
          errx(EX_USAGE, "thread count exceeds number of processors");*/
        Num_threads = sig;
        /*omp_set_num_threads(Num_threads);
        #endif*/
        break;
      case 'n' :
        Nodes_filename = optarg;
        break;
      case 'q' :
        Quick_mode = true;
        break;
      case 'm' :
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive minimum hit count");
        Minimum_hit_count = sig;
        break;
      case 'f' :
        Fastq_input = true;
        break;
      case 'c' :
        Only_classified_kraken_output = true;
        break;
      case 'C' :
        Print_classified = true;
        Classified_output_file = optarg;
        break;
      case 'U' :
        Print_unclassified = true;
        Unclassified_output_file = optarg;
        break;
      case 'o' :
        Kraken_output_file = optarg;
        break;
      case 'u' :
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive work unit size");
        Work_unit_size = sig;
        break;
      case 'M' :
        Populate_memory = true;
        break;
      default:
        usage();
        break;
    }
  }

  if (DB_filename.empty()) {
    cerr << "Missing mandatory option -d" << endl;
    usage();
  }
  if (Index_filename.empty()) {
    cerr << "Missing mandatory option -i" << endl;
    usage();
  }
  if (Nodes_filename.empty() && ! Quick_mode) {
    cerr << "Must specify one of -q or -n" << endl;
    usage();
  }
  if (optind == argc) {
    cerr << "No sequence data files specified" << endl;
  }
}

void usage(int exit_code) {
  cerr << "Usage: classify [options] <fasta/fastq file(s)>" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "* -d filename      Kraken DB filename" << endl
       << "* -i filename      Kraken DB index filename" << endl
       << "  -n filename      NCBI Taxonomy nodes file" << endl
       << "  -o filename      Output file for Kraken output" << endl
       << "  -t #             Number of threads" << endl
       << "  -u #             Thread work unit size (in bp)" << endl
       << "  -q               Quick operation" << endl
       << "  -m #             Minimum hit count (ignored w/o -q)" << endl
       << "  -C filename      Print classified sequences" << endl
       << "  -U filename      Print unclassified sequences" << endl
       << "  -f               Input is in FASTQ format" << endl
       << "  -c               Only include classified reads in output" << endl
       << "  -M               Preload database files" << endl
       << "  -h               Print this message" << endl
       << endl
       << "At least one FASTA or FASTQ file must be specified." << endl
       << "Kraken output is to standard output by default." << endl;
  exit(exit_code);
}
