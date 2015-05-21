#include <iostream>
#include <sstream>

#include "utils.h"

//Split a string with the given delimiter
std::vector<std::string> split(std::string &s, char delim) {

    /*
    std::vector<std::string> elements;

    std::stringstream stream(s);
    std::string element;
    while (std::getline(stream, element, delim)) {
        elements.push_back(element);
    }

    return elements;
    */

    // implementation 2.0
    std::vector<std::string> elements2;
    auto it = s.begin();
    size_t start = 0;
    size_t next = 0;
    while ((next = s.find_first_of(delim, start)) != std::string::npos) {
        elements2.emplace_back(s, start, (next-start));
        //std::cerr << "elements.back() = " << elements.back() << ", delim = [" << delim << "]\n";
        //std::cerr << "start = " << start << ", next = " << next << "\n";
        start = next + 1;
    }
    elements2.emplace_back(s, start, (s.length()-start));
    /*
    if ( elements != elements2 ) {
        std::cerr << "elements vectors differ\n";
        std::cerr << "sizes = " << elements.size() << ", " << elements2.size() << "\n";
    }
    */
    return elements2;
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
