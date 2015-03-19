/*
  poisson pmf

 */


#ifndef _POISSON_H_
#define _POISSON_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
  Description & parameters :
  IN
    lambda: Poisson parameter (mean)
    left, right: left and right bounds
  OUT
    prob: Poisson probabilities from left to right bounds
    weight: weight value, i.e., exact pmf is given by prob[i]/weight
 */

void poisson_prob(double lambda, int left, int right,
		double *prob, double *weight);

/*
  Description: compute the right bound of Poisson range 
               for a given tolerance error
  Parameters:
   IN
     lambda: Poisson rate (mean)
     eps: error tolerance
   OUT
     right bound is a return value
 */

int poisson_rightbound(double lambda, double eps);

#ifdef __cplusplus
}
#endif

#endif
