#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>

#include <time.h>
#include <sys/resource.h>
#include <sys/time.h>

std::vector<std::string> split(std::string &, char);
void print_time_elapsed(std::string, struct timeval*, struct timeval*);

#endif