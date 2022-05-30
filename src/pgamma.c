#include <math.h>

double p_gamma(double a, double x, double loggamma_a);
double q_gamma(double a, double x, double loggamma_a);

double p_gamma(double a, double x, double loggamma_a) {
  int k;
  double result, term, previous;
  if (x >= 1+a) return 1 - q_gamma(a, x, loggamma_a);
  if (x == 0)   return 0;
  result = term = exp(a * log(x) - x - loggamma_a) / a;
  for (k=1; k<1000; k++) {
    term *= x / (a+k);
    previous = result;
    result += term;
    if (result == previous) return result;
  }
  return result;
}

double q_gamma(double a, double x, double loggamma_a) {
  int k;
  double result, w, temp, previous;
  double la, lb;
  la = 1; lb = 1 + x - a;
  if (x < 1+a) return 1 - p_gamma(a, x, loggamma_a);
  w = exp(a * log(x) - x - loggamma_a);
  result = w/lb;
  for (k=2; k<1000; k++) {
    temp = ((k-1-a)*(lb-la)+(k+x)*lb)/k;
    la = lb;
    lb = temp;
    w *= (k-1-a)/k;
    temp = w/(la*lb);
    previous = result;
    result += temp;
    if (result == previous) return result;
  }
  return result;
}

