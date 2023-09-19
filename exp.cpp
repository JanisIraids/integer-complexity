#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "exp_container.h"

int main(int argc, char** argv)
{
  bool outToFile = false;
  std::string oFileName;
  bool progress = true;
  const std::string NOID = "";
  std::string expID = NOID;
  bool printUsage = true;
  bool printHelp = false;
  std::vector<std::string> params;
  int argi = 1;
  bool badargs = false;
  // read args
  while (argi < argc && !badargs)
  {
    printUsage = false;
    if (argv[argi][0] == '-')
    {
      switch (argv[argi][1])
      {
      case 'h':
	printHelp = true;
	break;
      case 'f':
	argi++;
	outToFile = true;
	oFileName = std::string(argv[argi]);
	break;
      case 's':
	progress = false;
	break;
      default:
	badargs = true;
	argi = argc;
	break;
      }
    }
    else
    {
      expID = std::string(argv[argi]);
      argi++;
      while (argi < argc)
      {
	params.push_back(std::string(argv[argi]));
	argi++;
      }
    }
    argi++;
  }
  if (printUsage)
  {
    std::cout << "Usage: exp [OPTIONS] [EXP_ID [PARAMS]]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "    -h Print description of experiment/s" << std::endl;
    std::cout << "    -f <filename> Print results to the file (otherwise printed to console)" << std::endl;
    std::cout << "    -s Do not show progress messages" << std::endl;
    return 0;
  }

  // do work
  std::ostream* stream;
  if (outToFile)
    stream = new std::ofstream(oFileName.c_str(), std::ios::out);
  else
    stream = &(std::cout);
  ExpCont cont(*stream);
  if (printHelp)
  {
    if (expID == NOID)
      cont.printAll();
    else
      cont.printOne(expID);
  }
  else
    cont.performExperiment(expID, progress, params);
  if (outToFile)
  {
    (dynamic_cast<std::ofstream*>(stream))->close();
    delete stream;
  }
  return 0;
}
