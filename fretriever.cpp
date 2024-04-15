#include "fretriever.h"
#include <iostream>
const std::string FRetriever::fname = "/mnt/e/f.bin";

void FRetriever::init(const ui cachecount, const ui cachesize)
{
  this->cachecount = cachecount;
  this->cachesize = cachesize;
  caches = new unsigned char*[cachecount];
  cachestart = new ull[cachecount];
  for (ui i = 0; i < cachecount; i++)
  {
    cachestart[i] = -cachesize;
    caches[i] = new unsigned char[cachesize];
  }
  in.open(fname.c_str(), std::ios::in|std::ios::binary);
  in.seekg(0, std::ios_base::end);
  size = in.tellg();
}

FRetriever::FRetriever()
{
  init(1, 1024*1024*1024);
}

FRetriever::FRetriever(const ui cachecount)
{
  init(cachecount, 1024*1024*1024);
}

FRetriever::FRetriever(const ui cachecount, const ui cachesize)
{
  init(cachecount, cachesize);
}

FRetriever::~FRetriever()
{
  for (ui i = 0; i < cachecount; i++)
    delete[] caches[i];
  delete[] caches;
  delete[] cachestart;
  in.close();
}

unsigned char FRetriever::f(const ull x)
{
  unsigned char result;
  in.seekg(x, std::ios_base::beg);
  in.read(reinterpret_cast<char*>(&result), sizeof(unsigned char));
  return result;
}

unsigned char FRetriever::mbf(const ull x, const ui cacheid) // for multiple caches
{
  if (x - cachestart[cacheid] < 0 || x - cachestart[cacheid] >= cachesize)
  {
    in.seekg(x, std::ios_base::beg);
    ull readable = min(cachesize, size-x);
    in.read(reinterpret_cast<char*>(caches[cacheid]), sizeof(unsigned char)*readable);
    cachestart[cacheid] = x;
  }
  return caches[cacheid][x-cachestart[cacheid]];
}

unsigned char FRetriever::bf(const ull x) // buffered f, assume increasing values of x
{
  return mbf(x, 0);
}

unsigned char FRetriever::mbrf(const ull x, const ui cacheid) // for multiple caches
{
  if (x - cachestart[cacheid] < 0 || x - cachestart[cacheid] >= cachesize)
  {
    in.seekg(max(x-cachesize+1,0ULL), std::ios_base::beg);
    ull readable = min(cachesize, size-max(x-cachesize+1,0ULL));
    in.read(reinterpret_cast<char*>(caches[cacheid]), sizeof(unsigned char)*readable);
    cachestart[cacheid] = max(x-cachesize+1,0ULL);
  }
  return caches[cacheid][x-cachestart[cacheid]];
}

unsigned char FRetriever::brf(const ull x) // buffered f, assume decreasing values of x
{
  return mbrf(x, 0);
}

unsigned char FRetriever::fbig(const mpz_t bignum)
{
  char bignumstr[64];
  mpz_get_str(bignumstr, 10, bignum);
  ull n;
  sscanf(bignumstr, "%llu", &n);
  if (n < cachesize)
    return bf(n);
  return f(n);
}

unsigned char FRetriever::fbig(const mpz_class& bignum)
{
  return fbig(bignum.get_mpz_t());
}
