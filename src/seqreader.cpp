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

#include "kraken_headers.hpp"
#include "seqreader.hpp"

using namespace std;

namespace kraken {
  FastaReader::FastaReader(string filename) {
    int read_length=201;
    dna.seq.reserve(read_length);
    dna.id.reserve(read_length);
    dna.header_line.reserve(read_length);
    dna.quals.reserve(read_length);
    file.open(filename.c_str());
    if (file.rdstate() & ifstream::failbit) {
      err(EX_NOINPUT, "can't open %s", filename.c_str());
    }
    valid = true;
  }
  
  DNASequence FastaReader::next_sequence() {
    DNASequence dna_;
    return this->next_sequence(dna_);
  }

  DNASequence FastaReader::next_sequence(DNASequence &dna_) {
    //DNASequence dna;
    dna.seq.clear();
    dna.header_line.clear();
    dna.id.clear();
    dna.quals.clear();

    if (! valid || ! file.good()) {
      valid = false;
      return dna;
    }
    string line;

    if (linebuffer.empty()) {
      getline(file, line);
    }
    else {
      line = linebuffer;
      linebuffer.clear();
    }

    if (line[0] != '>') {
      warnx("malformed fasta file - expected header char > not found");
      valid = false;
      return dna;
    }
    dna.header_line = line.substr(1);
    dna.id = dna.header_line;
    /*istringstream seq_id(dna.header_line);
    seq_id >> dna.id;*/
    
    //ostringstream seq_ss;
    //cwilks: instead of getline, use fcntl to UNLOCK file and then do a read one byte at a time
    //looking for EOL
    while (file.good()) {
      getline(file, line);
      if (line[0] == '>') {
        linebuffer = line;
        break;
      }
      else {
        //seq_ss << line;
        dna.seq += line;
      }
    }
    //dna.seq = seq_ss.str();

    if (dna.seq.empty()) {
      warnx("malformed fasta file - zero-length record (%s)", dna.id.c_str());
      valid = false;
      return dna;
    }

    return dna;
  }

  bool FastaReader::is_valid() {
    return valid;
  }

  FastqReader::FastqReader(string filename) {
    int read_length=101;
    char* seq=(char*) calloc(read_length,sizeof(char));
    dna.seq.reserve(read_length);
    dna.id.reserve(read_length);
    dna.header_line.reserve(read_length);
    dna.quals.reserve(read_length);
    file.open(filename.c_str());
    if (file.rdstate() & ifstream::failbit) {
      err(EX_NOINPUT, "can't open %s", filename.c_str());
    }
    valid = true;
  }
  
  DNASequence FastqReader::next_sequence() {
    DNASequence dna_;
    return this->next_sequence(dna_);
  }
  
  DNASequence FastqReader::next_sequence(DNASequence &dna_) {
    char** seq=(char**) calloc(101,sizeof(char*));
    char** id=(char**) calloc(101,sizeof(char*));
    char** header_line=(char**) calloc(101,sizeof(char*));
    char** quals=(char**) calloc(101,sizeof(char*));
    return this->next_sequence(seq,id,header_line,quals);
  } 

  DNASequence FastqReader::next_sequence(char** seq,char** id,char** header_line,char** quals) {
    //DNASequence dna;
    dna.seq.clear();
    dna.header_line.clear();
    dna.id.clear();
    dna.quals.clear();
    if (! valid || ! file.good()) {
      valid = false;
      return dna;
    }

    string line;
    getline(file, line);
    if (line.empty()) {
      valid = false;  // Sometimes FASTQ files have empty last lines
      return dna;
    }
    if (line[0] != '@') {
      if (line[0] != '\r')
        warnx("malformed fastq file - sequence header (%s)", line.c_str());
      valid = false;
      return dna;
    }
    dna.header_line = line.substr(1);
    //istringstream line_ss(dna.header_line);
    
    //line_ss >> dna.id;
    dna.id=dna.header_line;
    //printf("b4 id copy\n");
    memcpy(*id,dna.id.c_str(),7);
    //printf("b4 header copy\n");
    memcpy(*header_line,dna.header_line.c_str(),7);
    getline(file, dna.seq);
    //printf("b4 seq copy\n");
    memcpy(*seq,dna.seq.c_str(),101);
 
    getline(file, line);
    if (line.empty() || line[0] != '+') {
      if (line[0] != '\r')
        warnx("malformed fastq file - quality header (%s)", line.c_str());
      valid = false;
      return dna;
    }
    getline(file, dna.quals);
    //printf("b4 quals copy\n");
    memcpy(*quals,dna.quals.c_str(),101);
    
    return dna;
  }

  bool FastqReader::is_valid() {
    return valid;
  }
} // namespace
