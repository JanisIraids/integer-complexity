#include "exp_container.h"

ExpCont::ExpCont(std::ostream& outs)
{
  experiments.push_back(std::make_pair(new GenHypo<true>(outs), new GenHypo<false>(outs)));
  experiments.push_back(std::make_pair(new SecondSmallest<true>(outs), new SecondSmallest<false>(outs)));
  experiments.push_back(std::make_pair(new ShortestExpr(outs), new ShortestExpr(outs)));
  experiments.push_back(std::make_pair(new GoodNums<true>(outs), new GoodNums<false>(outs)));
  experiments.push_back(std::make_pair(new Collapse<true>(outs), new Collapse<false>(outs)));
  experiments.push_back(std::make_pair(new Collapse2<true>(outs), new Collapse2<false>(outs)));
  experiments.push_back(std::make_pair(new Defect<true>(outs), new Defect<false>(outs)));
  experiments.push_back(std::make_pair(new SophiePrimes<true>(outs), new SophiePrimes<false>(outs)));
  experiments.push_back(std::make_pair(new AverageDigit<true>(outs), new AverageDigit<false>(outs)));
  experiments.push_back(std::make_pair(new ComplexityBruteForce<true>(outs), new ComplexityBruteForce<false>(outs)));
  experiments.push_back(std::make_pair(new FCalculator<true>(outs), new FCalculator<false>(outs)));
}

ExpCont::~ExpCont()
{
  for (std::vector<std::pair<Experiment*, Experiment*> >::iterator it = experiments.begin(); it < experiments.end(); it++)
  {
    delete it->first;
    delete it->second;
  }
}

void ExpCont::performExperiment(const std::string eid, const bool prog, const std::vector<std::string> args)
{
  for (std::vector<std::pair<Experiment*, Experiment*> >::iterator it = experiments.begin(); it < experiments.end(); it++)
    if (it->first->id() == eid)
    {
      if (prog)
	it->first->perform(args);
      else
	it->second->perform(args);
    }
}

void ExpCont::printAll()
{
  for (std::vector<std::pair<Experiment*, Experiment*> >::iterator it = experiments.begin(); it < experiments.end(); it++)
    printOne(it->first->id());
}

void ExpCont::printOne(const std::string eid)
{
  std::cout << eid << ":" << std::endl;
  for (std::vector<std::pair<Experiment*, Experiment*> >::iterator it = experiments.begin(); it < experiments.end(); it++)
    if (it->first->id() == eid)
      it->first->printUsage();
}
