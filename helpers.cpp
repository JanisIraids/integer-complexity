#include <cmath>
// for memory info
#include <sys/sysinfo.h>
#include <unistd.h>

#include "helpers.h"


double lc(const long long n, const int complexity)
{
  static const double l3 = log(3);
  return static_cast<double>(complexity)*l3/log(n);
}

bool is_prime(const long long n)
{
  if (n < 2)
    return false;
  for (long long i = 0; i*i <= n; i++)
    if (n % i == 0)
      return false;
  return true;
}

ull isqrt(const ull n)
{
  ull l = 0, u = 0x100000000;
  while (u - l > 1)
  {
    ull m = (l + u)/2;
    if (m*m < n)
      l = m;
    else
      u = m;
  }
  if (u*u == n)
    return u;
  else
    return l;
}


unsigned long long get_total_system_memory()
{
  long pages = get_avphys_pages();
  long page_size = getpagesize();
  return pages * page_size;
}

ull pmultiple(const ull n, const ull p)
{
  return (n-1)/p+1;
}
