#include <algorithm>
#include <cmath>

namespace jama{
namespace util{
class Maths {
/** sqrt(a^2 + b^2) without under/overflow. **/

public:
   static double hypot(double a, double b) {
	  double r;
	  if (std::abs(a) > std::abs(b)) {
		 r = b/a;
		 r = std::abs(a)*std::sqrt(1+r*r);
	  } else if (b != 0) {
		 r = a/b;
		 r = std::abs(b)*std::sqrt(1+r*r);
	  } else {
		 r = 0.0;
	  }
	  return r;
   }   

};
}
}