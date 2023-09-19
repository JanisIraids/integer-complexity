#include "experiment.h"

Experiment::Experiment(std::ostream& outs) : out(outs) {}
void Experiment::perform(const std::vector<std::string> args) {}

ExperimentWRet::ExperimentWRet(std::ostream& outs) : Experiment(outs) {}
void ExperimentWRet::perform(const std::vector<std::string> args)
{
  ret = new FRetriever();
  _perform(args);
  delete ret;
}

ShortestExpr::ShortestExpr(std::ostream& outs) : Experiment(outs) {}
void ShortestExpr::printShortest(const long long n)
{
  if (n < 6)
  {
    out << n;
    return;
  }
  unsigned char fn = ret->f(n);
  for (long long i = 2; i * i <= n; i++)
  {
    if (n % i == 0 && ret->f(n/i)+ret->mbf(i,0) <= fn)
    {
      out << "(";
      printShortest(n/i);
      out << ")*(";
      printShortest(i);
      out << ")";
      return;
    }
  }
  for (long long i = 1; i < n; i++)
  {
    if (ret->mbf(i,0) + ret->mbrf(n-i,1) <= fn)
    {
      printShortest(i);
      out << "+";
      printShortest(n-i);
      return;
    }
  }
}

void ShortestExpr::printUsage()
{
  std::cout << "Print one shortest expression for n" << std::endl;
}

std::string ShortestExpr::id()
{
  return "shortest";
}

void ShortestExpr::perform(const std::vector<std::string> args)
{
  ret = new FRetriever(2);
  long long n;
  std::istringstream(args[0]) >> n;
  printShortest(n);
  out << std::endl;
  delete ret;
}

bool grdbl(const double& lhs, const double& rhs)
{
  return rhs < lhs;
}

bool lll(const unsigned int& lhs, const unsigned int& rhs)
{
  return lhs < rhs;
}

bool cmpll(const unsigned long long& lhs, const unsigned long long& rhs)
{
  return lhs < rhs;
}
