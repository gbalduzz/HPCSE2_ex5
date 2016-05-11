//
// Created by giovanni on 07.03.16.
//
#ifndef EX2_MORTON_H
#define EX2_MORTON_H
using uint=unsigned int;

inline uint integer_mantissa(float x,float xmin,float xmax);
uint   interleave(uint x,uint y);

void morton(const int N,const float * const x,const float * const y,const float xmin,const float ymin, const float ext, uint * index)
//takes array or coords and return array of square's label where the particle resides
{
  #pragma omp parallel for schedule(static)
for(int i=0;i<N;i++){
    index[i]=interleave(integer_mantissa(x[i],xmin,ext+xmin),integer_mantissa(y[i],ymin,ext+ymin));
}
}


uint   interleave(uint x,uint y)
{
    // Interleave lower 16 bits of x and y, so the bits of x
    // are in the even positions and bits from y in the odd;
    // z gets the resulting 32-bit Morton Number.
    // x and y must initially be less than 65536.
    //From  https://graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN
    static const unsigned int B[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};
    static const unsigned int S[] = {1, 2, 4, 8};

    x = (x | (x << S[3])) & B[3];
    x = (x | (x << S[2])) & B[2];
    x = (x | (x << S[1])) & B[1];
    x = (x | (x << S[0])) & B[0];

    y = (y | (y << S[3])) & B[3];
    y = (y | (y << S[2])) & B[2];
    y = (y | (y << S[1])) & B[1];
    y = (y | (y << S[0])) & B[0];

    return x | (y << 1);
}

inline uint integer_mantissa(float x,float xmin,float xmax){
    //map x,y from [x/ymin,x/ymin+ext] to [0, 65536)
    return static_cast<uint>( (x-xmin)/(xmax-xmin)*65535 );
}

uint integer_mantissa(float x);
#endif //EX2_MORTON_H
