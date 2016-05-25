export uniform double e2p(const uniform double z_rev[],const uniform double z_imv[],
       const  uniform double** const uniform c_rev,
       const  uniform double** const uniform c_imv,
       const uniform int nexp, const uniform int order){
  double result=0;
  foreach(j=0 ... nexp){
    const double z_re = z_rev[j];
    const double z_im = z_imv[j];
    const double* c_re = c_rev[j];
    const double* c_im = c_imv[j];
    
    const double r2 = z_re*z_re+z_im*z_im;	   
    result+= c_re[0]*0.5*log(r2);
    const double inv_re = z_re / r2;
    const double inv_im = -z_im/r2;
    double pow_re = inv_re;
    double pow_im = inv_im;
    for(uniform int i=1; i<= order;i++){
      result += c_re[i]*pow_re-c_im[i]*pow_im;
      const double temp = pow_re*inv_re - pow_im*inv_im;
      pow_im = pow_re*inv_im + pow_im*inv_re;
      pow_re = temp;
    }
  }
 return reduce_add(result);
  }      
  
