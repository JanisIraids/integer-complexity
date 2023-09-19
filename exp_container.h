#pragma once

#include <string>
#include <vector>
#include <utility>
#include <ostream>
#include "experiment.h"

class ExpCont
{
  std::vector<std::pair<Experiment*, Experiment*> > experiments; // first is with PROG, second is with no PROG

public:
  ExpCont(std::ostream& outs);
  ~ExpCont();

  void performExperiment(const std::string eid, const bool prog, const std::vector<std::string> args);
  void printAll();
  void printOne(const std::string eid);
};
