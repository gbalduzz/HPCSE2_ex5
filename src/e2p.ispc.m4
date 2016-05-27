include(src/unroll.m4)

//non vectorized, either buffer or try to pack expansions and z^k in a vector
export uniform double e2p(const uniform double z_re,const uniform double z_im,
       const uniform double c_re[],
       const uniform double c_im[]){
  const uniform double r2 = z_re*z_re+z_im*z_im;
  uniform double result = c_re[0]*0.5*log(r2);
  const uniform double pow_re_1 = z_re / r2;
  const uniform double pow_im_1 = -z_im/r2;

  result += c_re[1]*pow_re_1-c_im[1]*pow_im_1;
  LUNROLL(i,2,eval(ORDER),`
  // loop iteration
    const uniform double TMP(pow_re,i) = TMP(pow_re,eval(i-1))*pow_re_1 - TMP(pow_im, eval(i-1))*pow_im_1;
    const uniform double TMP(pow_im,i) = TMP(pow_im, eval(i-1))*pow_re_1 + TMP(pow_re, eval(i-1))*pow_im_1;
    result += c_re[i]*TMP(pow_re,i)-c_im[i]*TMP(pow_im,i);
    ')

 return result;//reduce_add(result);
  }      
  
