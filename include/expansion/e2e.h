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


//do we need lim?
template<int l, int z_exp>
struct sum_z_term_to_l_coeff{
  static void inline execute(const double z_re, const double z_im,
                      const double *c_re, const double *c_im,
                      double *new_re, double *new_im) {
    constexpr int k = l - z_exp;
    new_re[l] += (z_re * c_re[k] - z_im * c_im[k]) * binomial<l - 1, k - 1>();
    new_im[l] += (z_re * c_im[k] + z_im * c_re[k]) * binomial<l - 1, k - 1>();
    sum_z_term_to_l_coeff<l - 1, z_exp>::execute(z_re, z_im, c_re, c_im, new_re, new_im);
  }
};

//end recursion sum_to_l_coeff and add first term of eq(4.15)
template<int l>
struct sum_z_term_to_l_coeff<l,l>{
  static void inline execute(const double z_re,const double z_im,
                    const double* c_re,const double* c_im,
                      double* new_re,double* new_im){
    new_re[l] += c_re[0]*z_re /l;
    new_im[l] += c_re[0]*z_im /l;
  }
};

template<>
struct sum_z_term_to_l_coeff<0,0>{
  static void inline execute(const double z_re,const double z_im,
                             const double* c_re,const double* c_im,
                             double* new_re,double* new_im){
    new_re[0] += c_re[0];
  }
};

template<int z_exp, int lim>
struct sum_z_term {
  static void inline execute(const double z_re, const double z_im,
                      const double zpow_re, const double zpow_im, const double *c_re,
                      const double *c_im, double *new_re, double *new_im) {
    sum_z_term_to_l_coeff<lim, z_exp>::execute(zpow_re, zpow_im, c_re, c_im, new_re, new_im);
    const double new_pow_re = z_re * zpow_re - z_im * zpow_im;
    const double new_pow_im = z_re * zpow_im - z_im * zpow_re;
    sum_z_term<z_exp + 1, lim>::execute(z_re, z_im,
                                        new_pow_re, new_pow_im,
                                        c_re, c_im, new_re, new_im);
  }
};

//end recursion sum_z_term
template<int lim>
struct sum_z_term<lim,lim> {
  static void inline execute(const double z_re, const double z_im,
                      const double zpow_re, const double zpow_im, const double *c_re,
                      const double *c_im, double *new_re, double *new_im) {
    sum_z_term_to_l_coeff<lim, lim>::execute(zpow_re, zpow_im, c_re, c_im, new_re, new_im);
  }
};

template<int ord>
void e2e(const double* c_re,const double* c_im, double z_re,double z_im,double* new_re,double* new_im){
  sum_z_term<0,ord>::execute(z_re,z_im,1.,0.,c_re,c_im,new_re,new_im);
}
