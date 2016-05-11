//
// Created by giovanni on 07.03.16.
//
#ifndef EX2_PROFILER_H
#define EX2_PROFILER_H
#include <unistd.h>
#include <ctype.h>
#include <stdlib.h>
#include <string>
#include<iostream>

class Profiler{
public:
    Profiler(std::string n0="",bool frozen=0);
    inline ~Profiler();
    void start();
    double stop(bool print=true);
private:
    std::string name;
    timespec t0,t1;
    bool mute=false;
    float elapsed();
};


Profiler::Profiler(std::string n0,bool frozen):name(n0){
    if(not frozen) clock_gettime(CLOCK_MONOTONIC, &t0);
}

inline Profiler::~Profiler() {
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if(not mute)
    std::cout<<"function "<<name<<"\t took "<<elapsed()<< "\t s to complete"<<std::endl;
}

float Profiler::elapsed() {
    return (t1.tv_nsec-t0.tv_nsec)*1e-9+(t1.tv_sec-t0.tv_sec);
}

void Profiler::start() {
    clock_gettime(CLOCK_MONOTONIC, &t0);
}

double Profiler::stop(bool print) {
    clock_gettime(CLOCK_MONOTONIC, &t1);
    mute=true;
    double time=elapsed();
    if(print)
        std::cout<<"function "<<name<<"\t took "<<time<< "\t s to complete"<<std::endl;
    return time;
}


#endif //EX2_PROFILER_H
