#pragma once
#include <utility>

typedef struct genSet {
  long n;            // number of nodes
  int conn;          // maximum edge connectivity
  double p;          // prob of edge between nodes
  double r;          // geometric insertion/removal
  int max_appearances;  // the maximum number of times an edge can show up
                            // in the stream. 0 for no limit.
  std::string out_file; // file to write stream
  std::string cumul_out_file; // file to write cumul graph
  genSet(long n, int conn, double p, double r, int max_appearances,
         std::string out_file, std::string cumul_out_file)
         : n(n), conn(conn), p(p), r(r), max_appearances
         (max_appearances), out_file(std::move(out_file)), cumul_out_file
         (std::move(cumul_out_file)) {}
} GraphGenSettings;

/**
 * Generates a 1024-node graph with approximately 60,000 edge insert/deletes.
 * Writes stream output to sample.txt
 * Writes cumulative output to cumul_sample.txt
 */
void generate_stream(const GraphGenSettings& settings =
      {512,3,0.03,0.5,0,"./sample.txt", "./cumul_sample.txt"});
