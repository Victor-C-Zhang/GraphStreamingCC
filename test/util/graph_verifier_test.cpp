#include <gtest/gtest.h>
#include "../../include/test/file_graph_verifier.h"

const std::string fname = __FILE__;
size_t pos = fname.find_last_of("\\/");
const std::string curr_dir = (std::string::npos == pos) ? "" : fname.substr(0,pos);

TEST(DeterministicToolsTestSuite, TestKruskal) {
  ASSERT_EQ(78,FileGraphVerifier::kruskal(curr_dir+"/../res/multiples_graph_1024.txt").size());
}

TEST(DeterministicToolsTestSuite, TestEdgeVerifier) {
  FileGraphVerifier verifier(curr_dir+"/../res/multiples_graph_1024.txt");
  // add edges of the form {i,2i}
  Edge edge[3];
  edge[0] = 2;
  for (int i = 2; i < 512; ++i) {
    edge[1] = i;
    edge[2] = 2*i;
    verifier.verify_edge(edge);
  }

  // throw on nonexistent edge
  edge[1] = 69;
  edge[2] = 420;
  ASSERT_THROW(verifier.verify_edge(edge), BadEdgeException);

  // throw on already-included edge
  edge[1] = 120;
  edge[2] = 240;
  ASSERT_THROW(verifier.verify_edge(edge), BadEdgeException);

  // throw on edge within the same set
  edge[1] = 250;
  edge[2] = 1000;
  ASSERT_THROW(verifier.verify_edge(edge), BadEdgeException);
}

TEST(DeterministicToolsTestSuite, TestCCVerifier) {
  FileGraphVerifier verifier (curr_dir+"/../res/multiples_graph_1024.txt");
  // {0}, {1}, and primes \in [521,1021] are CCs
  verifier.verify_cc(0);
  verifier.verify_cc(1);
  verifier.verify_cc(911);
  // add edges of the form {i,2i}
  Edge edge[3];
  edge[0] = 2;
  for (int i = 2; i < 512; ++i) {
    edge[1] = i;
    edge[2] = 2*i;
    verifier.verify_edge(edge);
  }
  // nothing else is currently a CC
  for (int i = 2; i < 512; ++i) {
    ASSERT_THROW(verifier.verify_cc(i), NotCCException);
  }
}

// TODO: test on graphs with edge connectivity >2