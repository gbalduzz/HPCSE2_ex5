const static int N=1e6;
const static int k=20;
#include <iostream>
#include <string>
#include "include/Performance_measurer.h"
#include "include/file_IO.h"
#include "include/particles.h"
#include "include/profiler.h"
#include "include/expansion/P2E.h"
#include "include/expansion/E2P.h"
using std::cout; using std::endl;
using std::vector;


int main(int argc, char** argv) {
    Particles p(N);
    std::string filename="sources.bin";
    read_from_file(p,filename);
    cout<<"\nExpansion order: "<<k<<endl;
    double c_re[k+1],c_im[k+1];

    std::pair<double,double> tts = Performance(10,P2E<k>,p,c_re,c_im);
    cout<<"P2E time to solution: "<<tts.first<<" +- "<<tts.second<<endl;

    double z_re=1,z_im=1;
    double s;
    {
        Profiler prf("E2P");
        s=E2P<k>(z_re,z_im,c_re,c_im);
    }

    bool correct = check_from_file(s,k,"output.dat");
    if(correct) cout<<"Correct result!"<<endl;
}

