export uniform double p2p(const uniform double x[], const uniform double y[], const uniform double w[],
			  const uniform double xt, const uniform double yt, const uniform int N){
  double result=0;
  foreach(i=0 ... N){
    const double delta_x = x[i]-xt;
    const double delta_y = y[i]-yt;
    result += w[i] * log(delta_x*delta_x+delta_y*delta_y);
  }
  return reduce_add(result)*0.5;
}
