#ifndef random_
#define random_
#include <cstdlib>
#include <cmath>

/**
 * !\brief Random numbers for a number of distributions
 *
 * Random number generation using the C++ random number function
 */


class CRandom
{
 protected:
  static bool one; //!< used for Gaus
  static float angle; //!< used for Gaus
  static float x; //!< parameter
  static float const pi; //!< 3.14159
 public:
  static double Rndm();
  static float Gaus(float mean,float sigma);
  static float expDecayTime(float width);
  static float BreitWigner(float mean ,float width);
};


#endif
