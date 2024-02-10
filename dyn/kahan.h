#ifndef GRASS_KAHAN_H
#define GRASS_KAHAN_H

namespace dyn {
template <typename T = float> class Kahan {
  T a{}, e{};

public:
  Kahan() = default;
  explicit Kahan(T a) : a{a} {};
  constexpr T operator()() const { return a; }
  constexpr Kahan &add(T v) {
    auto y = v - e, t = a + y;
    return e = t - a - y, a = t, *this;
  }
  constexpr Kahan &operator+=(T v) { return add(v); }
};
} // namespace dyn

#endif // GRASS_KAHAN_H
