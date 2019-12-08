#ifndef TIMER_H
#define TIMER_H

#ifdef _WIN32

#include <Windows.h>

double get_wall_time(void)
{
    LARGE_INTEGER time,freq;
	//  Handle error
    if (!QueryPerformanceFrequency(&freq)) { 
		return 0;
	}
    if (!QueryPerformanceCounter(&time)) {
		return 0;
	}
    return (double)time.QuadPart / freq.QuadPart;
};

double get_cpu_time(void) {
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0) {
        return (double)(d.dwLowDateTime | ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    }
	else {
        return 0;
    }
};

//  Posix/Linux
#else

#include <time.h>
#include <sys/time.h>

double get_wall_time(void) 
{
    struct timeval time;
    if (gettimeofday(&time,NULL)) {
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
};

double get_cpu_time(void)
{
    return (double)clock() / CLOCKS_PER_SEC;
}
#endif

#endif