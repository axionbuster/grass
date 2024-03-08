#ifndef GRASS_ENV_H
#define GRASS_ENV_H

#include <cstdlib>
#include <optional>
#include <string>

namespace env {

std::optional<std::string> get(char const *sv) {
#ifdef _WIN32
  char *s;
  size_t n;
  // On Windows, _dupenv_s is the Microsoft-recommended way of reading an
  // environment variable. When the variable isn't found, the error (e) is 0 and
  // the count (n) is also 0 (also, the buffer [s] is set to NULL).
  if (auto e = _dupenv_s(&s, &n, sv); e || !n)
    return {};
  std::string t{s, n};
  // According to Microsoft, the buffer s must be freed using a call to `free`.
  std::free(s);
  return t;
#else
  // On non-Windows platforms, s may exist as an empty character if the
  // environment variable exists but the value is empty. In this case, the first
  // branch below is taken.
  if (char const *s = std::getenv(sv))
    return s;
  else
    return {};
#endif
}

} // namespace env

#endif // GRASS_ENV_H
