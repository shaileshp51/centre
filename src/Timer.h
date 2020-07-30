#ifndef INC_CENTRE_TIMER_H
#define INC_CENTRE_TIMER_H

#include "configcentre.h"

/// A class for timing the execution time of execution of a section of code
/// precision is milliseconds
class Timer {
public:
	Timer() {
		_start = Clock::now();
		_stop = Clock::now();
	}
	void start() {
		_start = Clock::now();
	}
	void stop() {
		_stop = Clock::now();
		_total = std::chrono::duration_cast<std::chrono::milliseconds>(
               _stop - _start).count();
	}
	double total_milliseconds() const {
		return _total;
	}
	double total() const {
		return _total/1000;
	}

private:
	Clock::time_point _start;
	Clock::time_point _stop;
	double _total;
};
#endif
