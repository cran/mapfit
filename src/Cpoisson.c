/*
  poisson 
 */

#include <stdlib.h>
#include <math.h>

static double NORMALQ_LOWER_Q = 3.0;
static double NORMALQ_UPPER_Q = 37.0;

static double NORMALQ_LOWER_LOGP = -689.0;
static double NORMALQ_UPPER_LOGP = -6.6;

static double NORMALQ_EPS = 1.0e-8;

static double LOG2PIOVER2 = 0.9189385332046727417803297364; // log(2pi) / 2

// POISSON_LAMBDA_MIN is lower bound to use tail approximation. If
// lambda is less than POISSON_LAMBDA_MIN, the right bound of Poisson
// can be founded by the simple iteration for computing cumulative
// Poisson prob. POISSON_RIGHT_MAX is used to determine the upper
// bound to search the Poisson right bound by the simple iteration. If
// n of Poi(n) attans POISSON_RIGHT_MAX, the iteration is stopped. That
// is, the probability
//    sum_{n=POISSON_RIGHT_MAX+1}^Infinity Poi(n)
// is ignored.

#define POISSON_LAMBDA_MIN 3.0
#define POISSON_RIGHT_MAX  23

#define PROB(x) prob[x-left]

/*
    Description
       return a tail probability of standard normal distribution
    Parameters
     IN
       x: input value
     OUT
       return value
*/

double normalt(double x) {

	// work variables
	double x2, tmp, sum;
	x2 = x*x;
	tmp = x;
	sum = 1.0 / tmp;
	tmp = tmp * x2;
	sum = sum - 1.0 / tmp;
	tmp = tmp * x2;
	sum = sum + 3.0 / tmp;
	tmp = tmp * x2;
	sum = sum - 15.0 / tmp;
	tmp = tmp * x2;
	sum = sum + 105.0 / tmp;
	return (log(sum) - x2/2.0 - LOG2PIOVER2);
}

/*
    Description
        return a quantile of standard normal distribution
    Parameters
      IN
       p: probability for the quantile
      OUT
       return value
*/

double normalq(double p) {

	// work variables
	double leps;
	double l, u, m, fm;

	leps = log(p);
	if (leps > NORMALQ_UPPER_LOGP || leps < NORMALQ_LOWER_LOGP) {
		//	  exit(EXIT_FAILURE);
		return 0.0;
	}
	l = NORMALQ_LOWER_Q;
	u = NORMALQ_UPPER_Q;
	m = (l + u) / 2.0;
	fm = normalt(m) - leps;
	while (fabs(fm) > NORMALQ_EPS) {
		if (fm > 0.0) {
			l = m;
		} else {
			u = m;
		}
		m = (l + u)/2.0;
		fm = normalt(m) - leps;
	}
	return m;
}

/*
  ! Description & parameters :
  !  IN
  !    lambda: Poisson parameter (mean)
  !    left, right: left and right bounds
  !  OUT
  !    prob: Poisson probabilities from left to right bounds
  !    weight: weight value, i.e., exact pmf is prob[i]/weight
 */

void poisson_prob(double lambda, int left, int right,
		double *prob, double *weight) {

	// work variables
	int mode, j, t, s;

	mode = (int) lambda;
	if (mode >= 1) {
		PROB(mode) = exp(-lambda + mode * log(lambda)
		- LOG2PIOVER2
		- (mode + 1.0/2.0)
		* log(mode) + mode);
	} else {
		PROB(mode) = exp(-lambda);
	}
	// -- down --
	for (j=mode; j>=left+1; j--) {
		PROB(j-1) = PROB(j) * j / lambda;
	}
	// -- up --
	for (j=mode; j<=right-1; j++) {
		PROB(j+1) = PROB(j) * lambda / (j+1);
	}
	// -- compute W --
	*weight = 0.0;
	s = left;
	t = right;
	while (s < t) {
		if (PROB(s) <= PROB(t)) {
			*weight += PROB(s);
			s++;
		} else {
			*weight += PROB(t);
			t--;
		}
	}
	*weight += PROB(s);
}

/*
  ! Description: compute the right bound of Poisson range
  !              for a given error tolerance
  !
  ! Parameters:
  !   IN
  !    lambda: Poisson rate (mean)
  !    epsi: error tolerance
  !   OUT
  !    right bound is a return value
 */

int poisson_rightbound(double lambda, double epsi) {

	// work variables
	int k;
	int right;
	double z, tmp;
	double total;

	if (lambda < POISSON_LAMBDA_MIN) {
		tmp = exp(-lambda);
		total = tmp;
		right = 0;
		for (k=1; k<=POISSON_RIGHT_MAX; k++) {
			right++;
			tmp = tmp * lambda / right;
			total += tmp;
			if (total + epsi >= 1.0)
				break;
		}
	} else {
		z = normalq(epsi);
		tmp = z + sqrt(4.0 * lambda - 1.0);
		right = (int) (tmp * tmp / 4.0 + 1.0);
	}
	return right;
}

