#ifndef _PROGRESS_BAR_H__
#define _PROGRESS_BAR_H__

#include <time.h>


class ProgressBar
{
	public:

		ProgressBar();

		void start();
		void stop();

		void status(unsigned int current, unsigned int max);

	private:

		void set_buffer(float percent);
		const float time_elapsed();
		int my_clock_gettime(int b, struct timespec *ts);

		struct timespec *timer;          // timer (global)
		char buffer[80];
		unsigned int size;

};
#endif
