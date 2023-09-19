#pragma once

typedef unsigned int ui;
typedef unsigned long long ull;

template <typename T>
T min(const T left, const T right) __attribute__((const));

template <typename T>
T min(const T left, const T right)
{
  if (left < right)
    return left;
  else
    return right;
}

template <typename T>
T max(const T left, const T right) __attribute__((const));

template <typename T>
T max(const T left, const T right)
{
  if (left > right)
    return left;
  else
    return right;
}

template <typename T>
T gcd(const T a, const T b) __attribute__((const));

template <typename T>
T gcd(const T a, const T b)
{
  if (b != 0)
    return gcd(b, a%b);
  else
    return a;
}


double lc(const long long n, const int complexity) __attribute__((const));

bool is_prime(const long long n) __attribute__((const));

ull isqrt(const ull n) __attribute__((const));

// smallest k such that k*p >= n
ull pmultiple(const ull n, const ull p) __attribute__((const));


unsigned long long get_total_system_memory();
