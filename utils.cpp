#include <iostream>
#include <sstream>

#include "utils.h"

//Split a string with the given delimiter
std::vector<std::string> split(std::string &s, char delim) {
    
    std::vector<std::string> elements;
    
    std::stringstream stream(s);
    std::string element;
    while (std::getline(stream, element, delim)) {
        elements.push_back(element);
    }

    return elements;
}

//Print time elapsed in seconds
void print_time_elapsed(std::string desc, struct timeval* start, struct timeval* end) {
    
    /*
    struct timeval {
        time_t      tv_sec;
        suseconds_t tv_usec;    
    }*/
    struct timeval elapsed;
    
    if(start->tv_usec > end->tv_usec) {
        end->tv_usec += 1000000;
        end->tv_sec--;
    }
    elapsed.tv_usec = end->tv_usec - start->tv_usec;
    elapsed.tv_sec  = end->tv_sec  - start->tv_sec;
    float time_elapsed = (elapsed.tv_sec*1000000 + elapsed.tv_usec)/1000000.f;
    std::cout << desc << " Total Time Elapsed = " << time_elapsed << std::endl;

    return;
}