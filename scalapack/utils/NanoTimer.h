#pragma once
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>

class NanoTimer {
public:
   struct timespec start;

   NanoTimer() {
      clock_gettime(CLOCK_MONOTONIC,  &start);

   }
   double elapsedSeconds() {
      struct timespec now;
      clock_gettime(CLOCK_MONOTONIC,  &now);
      double time = (now.tv_sec - start.tv_sec) + (double) (now.tv_nsec - start.tv_nsec) * 1e-9;
      start = now;
      return time;
   }
    void toc(string label) {
        double elapsed = elapsedSeconds();
        cout << label << ": " << elapsed << "s" << endl;        
    }
};
