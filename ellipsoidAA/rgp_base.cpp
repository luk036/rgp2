#include "rgp_base.h"
#include <map>

typedef std::map<unsigned int, int> pmap;

double max(const AAF& af, pmap& pol)
{
  double sum = af.get_center();
  
  for (unsigned i=0; i< af.get_length(); i++) {
    //assert(coefficients[i] != 0);
    if (af.get_coeff(i) > 0.0) {
      sum += af.get_coeff(i);
      pol[af.get_index(i)] = 1;
    }
    else if (af.get_coeff(i) < 0.0) {
      sum -= af.get_coeff(i);
      pol[af.get_index(i)] = -1;
    }
    else {
      pol[af.get_index(i)] = 0;
    }
  }
  return sum;
}

double eval(const AAF& af, const pmap& pol)
{
    double sum = af.get_center();
    
    for (unsigned i=0; i< af.get_length(); i++) {
      unsigned int idx1 = af.get_index(i);
      // Note: if corresponding index cannot be found, throw an exception
      pmap::const_iterator pit = pol.find(idx1);
      if (pit == pol.end()) throw;
      int polarity = pit->second;
      sum += polarity * af.get_coeff(i);
    }
    return sum;
}

 
/** Constructor (for AAF -> double) */
template <> template <>
monomial<double>::monomial(const monomial<AAF>& mon, 
			   const pmap& pol)
  : _a(mon._a.size()), _b(eval(mon._b, pol))
{
  for (size_t i=0; i<_a.size(); ++i) {
    _a[i] = eval(mon._a[i], pol);
  }
}

template <class Arr>
void rgp_base::assess(const Arr& x)
{
  for (size_t i= 1; i<_M.size(); ++i) {
    pmap pol;
    _f_value = max(_M[i](x), pol);
    if (_f_value > 0) {
      posynomial<double> P(_M[i], pol);
      _subgradient = P.log_exp_gradient(x);
      _is_violated = true;
      return;
    }
  }
  
  _is_violated = false;
  pmap pol;
  _f_value = max(_M[0](x), pol);
  posynomial<double> P(_M[0], pol);
  //_f_value = _M[0].fvalue_with_gradient(x, _subgradient);
  _subgradient = P.log_exp_gradient(x);
}

template
void rgp_base::assess(const std::valarray<double>& x);
