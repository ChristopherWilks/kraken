CXX = g++
CXXFLAGS = -Wall -fopenmp -g #-O3
PROGS = db_sort set_lcas classify make_seqid_to_taxid_map db_shrink

#ifeq (1,$(WITH_THREAD_PROFILING))
CXXFLAGS += -DPER_THREAD_TIMING=1
#endif

.PHONY: all install clean

all: $(PROGS)

install: $(PROGS)
	cp $(PROGS) $(KRAKEN_DIR)/

clean:
	rm -f $(PROGS) *.o

db_shrink: krakendb.o quickfile.o

db_sort: krakendb.o quickfile.o

set_lcas: krakendb.o quickfile.o krakenutil.o seqreader.o

classify: krakendb.o quickfile.o krakenutil.o seqreader.o mockseqreader.o

make_seqid_to_taxid_map: quickfile.o

krakenutil.o: krakenutil.cpp krakenutil.hpp
	$(CXX) $(CXXFLAGS) -c krakenutil.cpp

krakendb.o: krakendb.cpp krakendb.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c krakendb.cpp

seqreader.o: seqreader.cpp seqreader.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c seqreader.cpp

mockseqreader.o: mockseqreader.cpp mockseqreader.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c mockseqreader.cpp

quickfile.o: quickfile.cpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c quickfile.cpp
