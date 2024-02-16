#ifndef GRASS_TREE_H
#define GRASS_TREE_H

#include <complex>

namespace dyn::tree {

template <std::floating_point F> class Node {
  std::complex<F> ll, gg, cm;
  F m{};
  bool in(std::complex<F> xy);
  Node *quadrant(unsigned char q);
};

} // namespace dyn::tree

#endif // GRASS_TREE_H
