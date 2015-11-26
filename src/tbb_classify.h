#ifndef TBB_CLASSIFY_H_
#define TBB_CLASSIFY_H_

#include <tbb/tbb.h>
#include <tbb/task_group.h>

class pclassify {
    void* args_;     
public:
	pclassify(const pclassify& p): args_(p.args_) {};
	pclassify(void* args):args_(args) {};
	void operator()();
};
#endif
        
