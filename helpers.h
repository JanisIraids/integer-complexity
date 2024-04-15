#pragma once

typedef unsigned int ui;
typedef unsigned long long ull;

template <typename T>
T min(T left, T right) __attribute__((const));

template <typename T>
T min(const T left, const T right)
{
  if (left < right)
    return left;
  else
    return right;
}

template <typename T>
T max(T left, T right) __attribute__((const));

template <typename T>
T max(const T left, const T right)
{
  if (left > right)
    return left;
  else
    return right;
}

template <typename T>
T gcd(T a, T b) __attribute__((const));

template <typename T>
T gcd(const T a, const T b)
{
  if (b != 0)
    return gcd(b, a%b);
  else
    return a;
}


double lc(long long n, int complexity) __attribute__((const));

bool is_prime(ull n) __attribute__((const));

ull isqrt(ull n) __attribute__((const));

// smallest k such that k*p >= n
ull pmultiple(ull n, ull p) __attribute__((const));


unsigned long long get_total_system_memory();
