#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include"progressbar.h"

ProgressBar::ProgressBar()
{
	timer = NULL;
	memset(buffer, ' ', 80);
	size = 50;
};

void ProgressBar::start()
{
	if ( timer == NULL )
	{
		timer = (struct timespec*) malloc (sizeof(struct timespec));
	}
	my_clock_gettime(CLOCK_REALTIME, timer);
}

void ProgressBar::status(unsigned int current, unsigned int max)
{
	float percent = (float) current / (float) max * 100.0;

	set_buffer(percent);
	float time = time_elapsed();
	float torun = (100.0 / percent * time) - time;
	printf("%u \t %6.2f%% %s %7.2f seconds left\r",current, percent, buffer, torun);
	fflush(NULL);
}

void ProgressBar::stop()
{
	set_buffer(100.0);
	printf("%6.2f%% %s took %.2f seconds              \n", 100.0, buffer, time_elapsed());
	fflush(NULL);
}


void ProgressBar::set_buffer(float percent)
{
	buffer[0] = '[';
	unsigned int steps = floor( percent / 100 * size);
	unsigned int pos = 1;
	for (; pos <= steps; pos++)
	{
		buffer[pos] = '-';
	}
	pos++;
	for (; pos < size; pos++)
	{
		buffer[pos] = ' ';
	}

	buffer[size] = ']';
	buffer[size +1] = '\0';


}

const float ProgressBar::time_elapsed()
{
	struct timespec *timer2 = (struct timespec*) malloc (sizeof(struct timespec));

	// save current time
	my_clock_gettime(CLOCK_REALTIME, timer2);

	// compute right time difference
	int sec = timer2->tv_sec - timer->tv_sec;
	int nsec = (timer2->tv_nsec - timer->tv_nsec);

	if (nsec < 0)
	{
		sec--;
		nsec = 1000000000 +nsec;
	}

	float ret = sec + (((float) nsec) / (1000*1000*1000));
	free(timer2);
	return ret;
}

int ProgressBar::my_clock_gettime(int b, struct timespec *ts)
{
	return clock_gettime(b,ts);
}



