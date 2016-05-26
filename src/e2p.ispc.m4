include(src/unroll.m4)

export uniform double e2p(const uniform double z_rev[],const uniform double z_imv[],
       const uniform double* const uniform c_rev[],
       const uniform double* const uniform c_imv[],
       const uniform int nexp){
  double result=0;
  foreach(j=0 ... nexp){
    const double z_re = z_rev[j];
    const double z_im = z_imv[j];
    const double* c_re = c_rev[j];
    const double* c_im = c_imv[j];
    
    const double r2 = z_re*z_re+z_im*z_im;	   
    result+= c_re[0]*0.5*log(r2);
    const double pow_re_1 = z_re / r2;
    const double pow_im_1 = -z_im/r2;

    result += c_re[1]*pow_re_1-c_im[1]*pow_im_1;
    LUNROLL(i,2,eval(ORDER),`
    // loop iteration
	    const double TMP(pow_re,i) = TMP(pow_re,eval(i-1))*pow_re_1 - TMP(pow_im, eval(i-1))*pow_im_1;
	    const double TMP(pow_im,i) = TMP(pow_im, eval(i-1))*pow_re_1 + TMP(pow_re, eval(i-1))*pow_im_1;
	    result += c_re[i]*TMP(pow_re,i)-c_im[i]*TMP(pow_im,i);
	    ')
  }
 return reduce_add(result);
  }      
  
