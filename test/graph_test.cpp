#include <gtest/gtest.h>
#include <fstream>
#include <algorithm>
#include "../include/graph.h"
#include "../graph_worker.h"
#include "../include/test/file_graph_verifier.h"
#include "../include/test/graph_gen.h"
#include "../include/test/write_configuration.h"

/**
 * For many of these tests (especially for those upon very sparse and small graphs)
 * we allow for a certain number of failures per test.
 * This is because the responsibility of these tests is to quickly alert us 
 * to “this code is very wrong” whereas the statistical testing is responsible 
 * for a more fine grained analysis.
 * In this context a false positive is much worse than a false negative.
 * With 2 failures allowed per test our entire testing suite should fail 1/5000 runs.
 */

// We create this class and instantiate a paramaterized test suite so that we
// can run these tests both with the GutterTree and with StandAloneGutters
class GraphTest : public testing::TestWithParam<bool> {

};
INSTANTIATE_TEST_SUITE_P(GraphTestSuite, GraphTest, testing::Values(false));

TEST_P(GraphTest, SmallGraphConnectivity) {
  write_configuration(GetParam());
  const std::string fname = __FILE__;
  size_t pos = fname.find_last_of("\\/");
  const std::string curr_dir = (std::string::npos == pos) ? "" : fname.substr(0, pos);
  std::ifstream in{curr_dir + "/res/multiples_graph_1024.txt"};
  node_id_t num_nodes;
  edge_id_t m;
  int r;
  in >> num_nodes >> m >> r;
  Edge edge[r+1];
  Graph g{num_nodes, r};
  while (m--) {
    in >> edge[0];
    for (int i = 1; i <= edge[0]; ++i) {
      in >> edge[i];
    }
    g.update({edge, INSERT});
  }
  g.set_verifier(std::make_unique<FileGraphVerifier>(curr_dir + "/res/multiples_graph_1024.txt"));
  ASSERT_EQ(78, g.connected_components().size());
}

TEST(GraphTest, IFconnectedComponentsAlgRunTHENupdateLocked) {
  write_configuration(false);
  const std::string fname = __FILE__;
  size_t pos = fname.find_last_of("\\/");
  const std::string curr_dir = (std::string::npos == pos) ? "" : fname.substr(0, pos);
  std::ifstream in{curr_dir + "/res/multiples_graph_1024.txt"};
  node_id_t num_nodes;
  edge_id_t m;
  int r;
  in >> num_nodes >> m >> r;
  Edge edge[r+1];
  Graph g{num_nodes, r};
  while (m--) {
    in >> edge[0];
    for (unsigned i = 1; i <= edge[0]; ++i) {
      in >> edge[i];
    }
    g.update({edge, INSERT});
  }
  g.set_verifier(std::make_unique<FileGraphVerifier>(curr_dir + "/res/multiples_graph_1024.txt"));
  g.connected_components();

  Edge test_edge[] = {2, 1, 2};
  ASSERT_THROW(g.update({test_edge, INSERT}), UpdateLockedException);
  ASSERT_THROW(g.update({test_edge, DELETE}), UpdateLockedException);
}

TEST_P(GraphTest, TestSupernodeRestoreAfterCCFailure) {
  write_configuration(false, GetParam());
  const std::string fname = __FILE__;
  size_t pos = fname.find_last_of("\\/");
  const std::string curr_dir = (std::string::npos == pos) ? "" : fname.substr(0, pos);
  std::ifstream in{curr_dir + "/res/multiples_graph_1024.txt"};
  node_id_t num_nodes;
  edge_id_t m;
  int r;
  in >> num_nodes >> m >> r;
  Edge edge[r+1];
  Graph g{num_nodes, r};
  while (m--) {
    in >> edge[0];
    for (unsigned i = 1; i <= edge[0]; ++i) {
      in >> edge[i];
    }
    g.update({edge, INSERT});
  }
  g.set_verifier(std::make_unique<FileGraphVerifier>(curr_dir + "/res/multiples_graph_1024.txt"));
  g.should_fail_CC();

  // flush to make sure copy supernodes is consistent with graph supernodes
  g.gts->force_flush();
  GraphWorker::pause_workers();
  Supernode* copy_supernodes[num_nodes];
  for (node_id_t i = 0; i < num_nodes; ++i) {
    copy_supernodes[i] = Supernode::makeSupernode(*g.supernodes[i]);
  }

  ASSERT_THROW(g.connected_components(true), OutOfQueriesException);
  for (node_id_t i = 0; i < num_nodes; ++i) {
    for (int j = 0; j < copy_supernodes[i]->get_num_sktch(); ++j) {
      ASSERT_TRUE(*copy_supernodes[i]->get_sketch(j) ==
                *g.supernodes[i]->get_sketch(j));
    }
  }
}

TEST_P(GraphTest, TestCorrectnessOnSmallRandomGraphs) {
  write_configuration(GetParam());
  int num_trials = 5;
  int edge_conn = 2; // TODO: try on graphs with higher connectivity
  while (num_trials--) {
    generate_stream({1024, edge_conn, 0.002,0.5,0,"./sample.txt",""
                                                                 "./cumul_sample.txt"});
    std::ifstream in{"./sample.txt"};
    node_id_t n;
    edge_id_t m;
    int r;
    in >> n >> m >> r;
    Graph g{n, edge_conn};
    UpdateType type;
    Edge buf[edge_conn + 1];
    while (m--) {
      in >> buf[0];
      type = static_cast<UpdateType>(buf[0] & 1);
      buf[0] >>= 1;
      for (int i = 1; i <= buf[0]; ++i) {
        in >> buf[i];
      }
      g.update({buf, type});
    }

    g.set_verifier(std::make_unique<FileGraphVerifier>("./cumul_sample.txt"));
    g.connected_components();
  }
}

TEST_P(GraphTest, TestCorrectnessOnSmallSparseGraphs) {
  write_configuration(GetParam());
  int num_trials = 5;
  int edge_conn = 2; // TODO: try on edge conn > 2
  while(num_trials--) {
    generate_stream({1024, edge_conn, 0.002,0.5,0,"./sample.txt",""
                                                                 "./cumul_sample.txt"});
    std::ifstream in{"./sample.txt"};
    node_id_t n;
    edge_id_t m;
    int r;
    in >> n >> m >> r;
    Graph g{n, edge_conn};
    UpdateType type;
    Edge buf[edge_conn + 1];
    while (m--) {
      in >> buf[0];
      type = static_cast<UpdateType>(buf[0] & 1);
      buf[0] >>= 1;
      for (int i = 1; i <= buf[0]; ++i) {
        in >> buf[i];
      }
      g.update({buf, type});
    }

    g.set_verifier(std::make_unique<FileGraphVerifier>("./cumul_sample.txt"));
    g.connected_components();
  } 
}

TEST_P(GraphTest, TestCorrectnessOfReheating) {
  write_configuration(GetParam());
  int num_trials = 5;
  int edge_conn = 2;
  while (num_trials--) {
    generate_stream({1024,edge_conn,0.002,0.5,0,"./sample.txt","./cumul_sample"
                                                          ".txt"});
    std::ifstream in{"./sample.txt"};
    node_id_t n;
    edge_id_t m;
    int r;
    in >> n >> m >> r;
    auto *g = new Graph (n, edge_conn);
    printf("number of updates = %lu\n", m);
    UpdateType type;
    Edge buf[edge_conn + 1];
    while (m--) {
      in >> buf[0];
      type = static_cast<UpdateType>(buf[0] & 1);
      buf[0] >>= 1;
      for (int i = 1; i <= buf[0]; ++i) {
        in >> buf[i];
      }
      g->update({buf, type});
    }
    g->write_binary("./out_temp.txt");
    g->set_verifier(std::make_unique<FileGraphVerifier>("./cumul_sample.txt"));
    std::vector<std::set<node_id_t>> g_res;
    g_res = g->connected_components();
    printf("number of CC = %lu\n", g_res.size());
    delete g; // delete g to avoid having multiple graphs open at once. Which is illegal.

    Graph reheated {"./out_temp.txt"};
    reheated.set_verifier(std::make_unique<FileGraphVerifier>("./cumul_sample.txt"));
    auto reheated_res = reheated.connected_components();
    printf("number of reheated CC = %lu\n", reheated_res.size());
    ASSERT_EQ(g_res.size(), reheated_res.size());
    for (unsigned i = 0; i < g_res.size(); ++i) {
      std::vector<node_id_t> symdif;
      std::set_symmetric_difference(g_res[i].begin(), g_res[i].end(),
          reheated_res[i].begin(), reheated_res[i].end(),
          std::back_inserter(symdif));
      ASSERT_EQ(0, symdif.size());
    }
  }
}

// Test the multithreaded system by specifiying multiple
// Graph Workers of size 2. Ingest a stream and run CC algorithm.
TEST_P(GraphTest, MultipleInserters) {
  write_configuration(GetParam(), false, 4, 2);
  int num_trials = 5;
  int edge_conn = 2; // TODO: try on edge conn > 2
  while(num_trials--) {
    generate_stream({1024, edge_conn,0.002,0.5,0,"./sample.txt","./cumul_sample"
                                                          ".txt"});
    std::ifstream in{"./sample.txt"};
    node_id_t n;
    edge_id_t m;
    int r;
    in >> n >> m >> r;
    Graph g{n, edge_conn};
    UpdateType type;
    Edge buf[edge_conn + 1];
    while (m--) {
      in >> buf[0];
      type = static_cast<UpdateType>(buf[0] & 1);
      buf[0] >>= 1;
      for (int i = 1; i <= buf[0]; ++i) {
        in >> buf[i];
      }
      g.update({buf, type});
    }

    g.set_verifier(std::make_unique<FileGraphVerifier>("./cumul_sample.txt"));
    g.connected_components();
  } 
}

/* DISABLED because MatGraphVerifier is not implemented (read: a headache to
 * re-implement for hypergraphs)
TEST(GraphTest, TestQueryDuringStream) {
  write_configuration(false, false);
  { // test copying to disk
    generate_stream({1024, 0.002, 0.5, 0, "./sample.txt", "./cumul_sample.txt"});
    std::ifstream in{"./sample.txt"};
    node_id_t n;
    edge_id_t m;
    in >> n >> m;
    Graph g(n);
    MatGraphVerifier verify(n);

    int type;
    node_id_t a, b;
    edge_id_t tenth = m / 10;
    for(int j = 0; j < 9; j++) {
      for (edge_id_t i = 0; i < tenth; i++) {
        in >> type >> a >> b;
        g.update({{a,b}, (UpdateType)type});
        verify.edge_update(a, b);
      }
      verify.reset_cc_state();
      g.set_verifier(std::make_unique<MatGraphVerifier>(verify));
      g.connected_components(true);
    }
    m -= 9 * tenth;
    while(m--) {
      in >> type >> a >> b;
      g.update({{a,b}, (UpdateType)type});
      verify.edge_update(a, b);
    }
    verify.reset_cc_state();
    g.set_verifier(std::make_unique<MatGraphVerifier>(verify));
    g.connected_components();
  }

  write_configuration(false, true);
  { // test copying in memory
    generate_stream({1024, 0.002, 0.5, 0, "./sample.txt", "./cumul_sample.txt"});
    std::ifstream in{"./sample.txt"};
    node_id_t n;
    edge_id_t m;
    in >> n >> m;
    Graph g(n);
    MatGraphVerifier verify(n);

    int type;
    node_id_t a, b;
    edge_id_t tenth = m / 10;
    for(int j = 0; j < 9; j++) {
      for (edge_id_t i = 0; i < tenth; i++) {
        in >> type >> a >> b;
        g.update({{a,b}, (UpdateType)type});
        verify.edge_update(a, b);
      }
      verify.reset_cc_state();
      g.set_verifier(std::make_unique<MatGraphVerifier>(verify));
      g.connected_components(true);
    }
    m -= 9 * tenth;
    while(m--) {
      in >> type >> a >> b;
      g.update({{a,b}, (UpdateType)type});
      verify.edge_update(a, b);
    }
    verify.reset_cc_state();
    g.set_verifier(std::make_unique<MatGraphVerifier>(verify));
    g.connected_components();
  }
}
 */
