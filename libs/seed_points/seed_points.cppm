module;

#include <concepts>
#include <random>

#include <xoshiro.h>

#include <UBK/macros.hpp>

export module UBKLib.seed_points;

import UBKLib.utils;

/* 
 * Generators for starting points for field lines
*/

export namespace ubk {

template<std::floating_point T, class Generator = xso::rng>
class UniformEquatorGenerator {
public:

  [[nodiscard]] constexpr T
  minDistSquared(void) const {
    return m_minDist * m_minDist;
  }

  [[nodiscard]] constexpr T
  maxDistSquared(void) const {
    return m_maxDist * m_maxDist;
  }

  UniformEquatorGenerator(T minDist, T maxDist) 
  : m_minDist(minDist),
    m_maxDist(maxDist),
    m_generator(std::uniform_real_distribution<T>(minDist, maxDist)) {
    check(minDist < maxDist);
  }

  [[nodiscard]] Vector3<T>
  gen(void) {
    while(true) {
      Vector3<T> retVal = m_generator.gen();
      retVal.z = 0;
      T distSquared = retVal.ampSquared();
      if (distSquared < maxDistSquared() &&
          distSquared > minDistSquared()) {
        return retVal;
      }
    }
  }

private:
  T m_minDist = static_cast<T>(1.01);
  T m_maxDist = 15;
  RandomVector3Generator<T, 
    std::uniform_real_distribution<T>,
    std::uniform_real_distribution<T>,
    std::uniform_real_distribution<T>,
    Generator
  > m_generator;
  
};

}
