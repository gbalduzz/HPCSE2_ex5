include(src/unroll.m4)

#pragma once

template<int n, int k>
struct Binomial
{
  constexpr static int value =  (Binomial<n-1,k-1>::value + Binomial<n-1,k>::value);
};

template<>
struct Binomial<0,0>
{
  constexpr static int value = 1;
};

template<int n>
struct Binomial<n,0>
{
  constexpr static int value = 1;
};

template<int n>
struct Binomial<n,n>
{
  constexpr static int value = 1;
};

template<int n, int k>
inline int binomial()
{
  return Binomial<n,k>::value;
}

void e2e(const  double cr[], const  double ci[],
                const  double zr_1, const  double zi_1,
                 double newr[],  double newi[]){
  newr[0] += cr[0];

  LUNROLL(pow,1,eval(ORDER-1),`

    newr[pow] += cr[0]*TMP(zr,pow)/pow;
    newi[pow] += cr[0]*TMP(zi,pow)/pow;

    LUNROLL(l,eval(pow+1),ORDER,`
        newr[l] += (TMP(zr,pow) * cr[eval(l-pow)] - TMP(zi,pow * ci[eval(l-pow)])) * binomial<l - 1, eval(l-pow) - 1>();
        newi[l] += (TMP(zr,pow) * ci[eval(l-pow)] + TMP(zi,pow) * cr[eval(l-pow)]) * binomial<l - 1, eval(l-pow) - 1>();

    ')
    const  double TMP(zr,eval(pow+1)) = TMP(zr, pow)*zr_1 - TMP(zi, pow)*zi_1;
    const  double TMP(zi,eval(pow+1)) = TMP(zi, pow)*zr_1 + TMP(zr, pow)*zi_1;
  ')

  newr[ORDER] += cr[0]*TMP(zr,ORDER)/ORDER;
  newi[ORDER] += cr[0]*TMP(zi,ORDER)/ORDER;
}
