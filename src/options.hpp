/**
 * \file options.hpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 */
#pragma once

#include <string>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

struct standardOptions {
  std::string instance;
  unsigned seed;
};

void add_std_options(boost::program_options::options_description &, standardOptions &);
unsigned get_terminal_width();

struct std_description : public boost::program_options::options_description {
  typedef boost::program_options::options_description Base;
  std_description(const std::string &caption, standardOptions &);
};
