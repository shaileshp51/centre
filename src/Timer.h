#ifndef INC_CENTRE_TIMER_H
#define INC_CENTRE_TIMER_H
/// Class used to get timing information.
/** Under the hood the function gettimeofday is used by default, which has
 * microsecond resolution. If TIMER is defined, clock_gettime is used instead
 * which has nanosecond resolution and is superior to gettimeofday in many 
 * ways, but requires linking to another library for some versions of GLIB
 * an so is less portable.
 */
class Timer {
public:
	Timer();
	void start() {
		getWallTime(start_sec_, start_ns_);
	}
	void stop();
	double total() const {
		return total_;
	}
	void writeTiming(int, const char*, double) const;
	void writeTiming(int i, const char *h) const {
		return writeTiming(i, h, 0.0);
	}
private:
	void getWallTime(int&, int&);
	int start_sec_;
	int start_ns_;
	double total_;
};
#endif
