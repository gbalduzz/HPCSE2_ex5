include(src/unroll.m4)

export void p2e(const uniform double xv[],const uniform double yv[],const uniform double wv[],
	       const uniform double xcom, const uniform double ycom,
	       uniform double cr_out[],uniform double ci_out[],const uniform int N){
  double cr[eval(ORDER+1)];
  double ci[eval(ORDER+1)];
 
  foreach(i=0 ... N){
    const double zr_1=xv[i]-xcom;
    const double zi_1=yv[i]-ycom;
    const double w = wv[i];

    cr[0]+=w;
    cr[1] -= w*zr_1;
    ci[1] -= w*zi_1;
    LUNROLL(k,2,eval(ORDER),`
      //loop iteration
      const double TMP(zr,k) = TMP(zr,eval(k-1))*zr_1 - TMP(zi,eval(k-1))*zi_1;
      const double TMP(zi,k) = TMP(zr,eval(k-1))*zi_1 + TMP(zi,eval(k-1))*zr_1;
      cr[k] -= w*TMP(zr,k);
      ci[k] -= w*TMP(zi,k);
     ')
    }	     
    cr_out[0] = reduce_add(cr[0]);
    LUNROLL(k,2,ORDER,`
      cr_out[k] = reduce_add(cr[k])/k;
      ci_out[k] = reduce_add(ci[k])/k;
    ')
  }
