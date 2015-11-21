#ifndef ThreadProfile_H_
#define ThreadProfile_H_

#include "timer.h"
using namespace std;

class ThreadProfile 
{
private:
	uint64_t ncpu_changeovers;
	uint64_t nnuma_changeovers;
	int current_cpu;
	int current_node;
	std::stringstream ss;
	std::string msg;
	Timer *timer;
	const char *prefix;

public:
		
/// Based on http://stackoverflow.com/questions/16862620/numa-get-current-node-core
	static void get_cpu_and_node(int& cpu, int& node) {
		unsigned long a,d,c;
		__asm__ volatile("rdtscp" : "=a" (a), "=d" (d), "=c" (c));
		node = (c & 0xFFF000)>>12;
		cpu = c & 0xFFF;
	}

	ThreadProfile(const char *prefix_)
	{
#ifdef PER_THREAD_TIMING
		prefix=prefix_;
		ncpu_changeovers = 0;
		nnuma_changeovers = 0;
		current_cpu = 0;
		current_node = 0;
		get_cpu_and_node(current_cpu, current_node);
		
		//ss << "thread: " << tid << " time: ";
		ss << prefix << " thread: " << " time: ";
		msg = ss.str();
    	  	cout << prefix << " initialized timer/changeover tracker\n";
                timer=new Timer(std::cout, msg.c_str());
#endif
	}
	
	~ThreadProfile()
	{
		delete timer;
	}

	void update()
	{
#ifdef PER_THREAD_TIMING
			int cpu = 0, node = 0;
			get_cpu_and_node(cpu, node);
			if(cpu != current_cpu) {
				ncpu_changeovers++;
				current_cpu = cpu;
			}
			if(node != current_node) {
				nnuma_changeovers++;
				current_node = node;
			}
#endif
	}

	void finish()
	{
#ifdef PER_THREAD_TIMING
		ss.str("");
		ss.clear();
		ss << prefix << " thread: " << " cpu_changeovers: " << ncpu_changeovers << std::endl
		   << prefix << " thread: " << " node_changeovers: " << nnuma_changeovers << std::endl;
		std::cout << ss.str();
#endif
	}
};
#endif /* ThreadProfile_H_ */
