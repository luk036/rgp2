#ifndef RGP_BASE_H
#define RGP_BASE_H

#include "posynomial.h"
#include <valarray>
#include <vector>
#include <aaf/aa.h>

class rgp_base
{
  typedef std::valarray<double> Vec;
  
public:
  rgp_base(size_t n) : _subgradient(n) {}
  ~rgp_base() {}
  
  bool is_violated() const { return _is_violated; }
  double f_value() const { return _f_value; }
  const Vec& subgradient() const { return _subgradient; }
  
  template <class Arr>
  void assess(const Arr& x);
  
private:
  double _f_value;
  Vec	 _subgradient;
  bool   _is_violated;

protected:
  std::vector<posynomial<AAF> > _M; 
};

#endif 
