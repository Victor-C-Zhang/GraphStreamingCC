//
// Created by victor on 4/28/22.
//
#include <string>
#include "../../include/test/graph_gen.h"

void run_gen_test() {
  std::string fname = __FILE__;
  size_t pos = fname.find_last_of("\\/");
  std::string curr_dir = (std::string::npos == pos) ? "" : fname.substr(0, pos);
  generate_stream();
}

int main() {
  run_gen_test();
}