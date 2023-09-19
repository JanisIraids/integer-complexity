#pragma once

#define _FILE_OFFSET_BITS 64
#include <fstream>
#include <string>

#include <gmp.h>
#include <gmpxx.h>

#include "helpers.h"

class FRetriever
{
  const static std::string fname;
  std::ifstream in;
  unsigned char** caches;
  ull * cachestart;
  ull cachesize;
  ui cachecount;
  void init(const ui cachecount, const ui cachesize);

public:

  ull size;
  FRetriever();
  FRetriever(const ui cachecount);
  FRetriever(const ui cachecount, const ui cachesize);
  ~FRetriever();
  unsigned char f(const ull x);
  unsigned char mbf(const ull x, const ui cacheid); // for multiple caches
  unsigned char bf(const ull x); // buffered f, assume increasing values of x
  unsigned char mbrf(const ull x, const ui cacheid); // for multiple caches
  unsigned char brf(const ull x); // buffered f, assume decreasing values of x
  unsigned char fbig(const mpz_t bignum);
  unsigned char fbig(const mpz_class& bignum);
};
