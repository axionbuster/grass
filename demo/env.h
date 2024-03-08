#ifndef GRASS_ENV_H
#define GRASS_ENV_H

#include <cstdlib>
#include <optional>
#include <string>
#include <string_view>

namespace env {

std::optional<std::string> get(std::string_view sv) {
  if (auto s = std::getenv(sv.begin()))
    return s;
  else
    return {};
}

} // namespace env

#endif // GRASS_ENV_H
