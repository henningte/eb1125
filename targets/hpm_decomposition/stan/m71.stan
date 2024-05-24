vector m71(vector t, vector k_1, vector l_1, real s) {

  int N = size(t);
  vector[N] res;

  res = (1.0 - l_1) .* exp(-k_1 .* t);
  for(n in 1:N) {
    if(res[n] >= 1.0) {
      res[n] = 1.0 - s;
    }
  }

  return(res);

}
