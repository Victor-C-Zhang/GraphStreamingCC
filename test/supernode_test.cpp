#include <gtest/gtest.h>
#include <cmath>
#include <chrono>
#include <thread>
#include "../include/supernode.h"
#include "../include/graph_worker.h"
#include "update_list.h"

const long seed = 7000000001;

class GraphUpdateList {
private:
  std::vector<Edge> _edge;
  int edge_conn;
public:
  GraphUpdateList(int edge_conn) : edge_conn(edge_conn) {}
  void insert(Edge* edge) {
    assert(edge[0] <= edge_conn);
    for (int i = 0; i <= edge[0]; ++i) {
      _edge.push_back(edge[i]);
    }
    // back-pad nonmaximal edges
    for (int i = edge[0]; i < edge_conn; ++i) {
      _edge.push_back(0);
    }
  }

  struct Iterator {
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = Edge*;
    using pointer = Edge**;
    using reference = Edge*;
    // no references since an edge is by default a pointer type

  private:
    value_type _ptr;
    difference_type _stride;


  public:
    Iterator(value_type ptr, difference_type stride) : _ptr(ptr), _stride(stride) {}

    reference operator*() const { return _ptr; }
    pointer operator->() { return &_ptr; }
    Iterator& operator++() { _ptr += _stride; return *this; }
    Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
    friend bool operator==(const Iterator& a, const Iterator& b) {
      return (a._ptr == b._ptr) && (a._stride == b._stride);
    }
    friend bool operator != (const Iterator& a, const Iterator& b) {
      return (a._ptr != b._ptr) || (a._stride != b._stride);
    }
  };

  Iterator begin() { return {_edge.data(), edge_conn + 1}; }
  Iterator end() {return {_edge.data() + _edge.size(), edge_conn + 1};}
};

class SupernodeTestSuite : public testing::Test {
protected:
  const static int edge_conn = 2;
  const static uint64_t num_nodes = 2000;
  static GraphUpdateList* graph_edges;
  static GraphUpdateList* odd_graph_edges;
  static bool* prime;
  static void SetUpTestSuite() {
    srand(1000000007);
    graph_edges = new GraphUpdateList(edge_conn);
    odd_graph_edges = new GraphUpdateList(edge_conn);
    Edge edge[edge_conn + 1];
    edge[0] = 2;
    for (unsigned i=2;i<num_nodes;++i) {
      edge[1] = i;
      for (unsigned j = i*2; j < num_nodes; j+=i) {
        edge[2] = j;
        graph_edges->insert(edge);
        if ((j/i)%2) odd_graph_edges->insert(edge);
      }
    }

    // sieve
    prime = (bool*) malloc(num_nodes*sizeof(bool));
    std::fill(prime,prime+num_nodes,true);
    for (unsigned i = 2; i < num_nodes; i++) {
      if (prime[i]) {
        for (unsigned j = i * i; j < num_nodes; j += i) prime[j] = false;
      }
    }

    // configure
    Supernode::configure(num_nodes, edge_conn);
  }
  static void TearDownTestSuite() {
    delete graph_edges;
    delete odd_graph_edges;
    free(prime);
  }

  void SetUp() override {Supernode::configure(num_nodes, edge_conn);}
  void TearDown() override {}
};

GraphUpdateList* SupernodeTestSuite::graph_edges;
GraphUpdateList* SupernodeTestSuite::odd_graph_edges;
bool* SupernodeTestSuite::prime;

TEST_F(SupernodeTestSuite, TestConcatPairingFn) {
  uint32_t buf[5];
  const int n = 1000;
  const int edge_conn = 4;
  Supernode::configure(n, edge_conn);
  Supernode s = Supernode(12345);

  // right-pad with the last value if necessary
  auto pairing_fn = [](int x1, int x2, int x3, int x4) {
    uint128_t e = x1 + x2 * n + x3 * n * n + x4 * n * n * n;
    return e;
  };

  for (int i = 0; i < 100; ++i) {
    for (int j = i + 1; j < 100; ++j) {
      for (int k = j + 1; k < 100; ++k) {
        for (int l = k + 1; l < n; ++l) {
          s.inv_concat_tuple_fn(pairing_fn(i,j,k,l), buf);
          ASSERT_EQ(buf[0], 4);
          ASSERT_EQ(buf[1], i);
          ASSERT_EQ(buf[2], j);
          ASSERT_EQ(buf[3], k);
          ASSERT_EQ(buf[4], l);
        }
      }
    }
  }

  for (int i = 0; i < 100; ++i) {
    for (int j = i + 1; j < 100; ++j) {
      for (int k = j + 1; k < n; ++k) {
          s.inv_concat_tuple_fn(pairing_fn(i,j,k,k), buf);
          ASSERT_EQ(buf[0], 3);
          ASSERT_EQ(buf[1], i);
          ASSERT_EQ(buf[2], j);
          ASSERT_EQ(buf[3], k);
      }
    }
  }
}

TEST_F(SupernodeTestSuite, GIVENnoEdgeUpdatesIFsampledTHENnoEdgeIsReturned) {
  Supernode* s = Supernode::makeSupernode(seed);
  Edge edge[3];
  SampleSketchRet retcode;
  s->sample(edge, &retcode);
  ASSERT_EQ(retcode, ZERO) << "Did not get ZERO when sampling empty vector";
}

TEST_F(SupernodeTestSuite, IFsampledTooManyTimesTHENthrowOutOfQueries) {
  Supernode* s = Supernode::makeSupernode(seed);
  Edge edge[3];
  SampleSketchRet retcode;
  for (int i = 0; i < s->get_num_sktch(); ++i) {
    s->sample(edge, &retcode);
  }
  ASSERT_THROW(s->sample(edge, &retcode), OutOfQueriesException);
}

TEST_F(SupernodeTestSuite, TestSampleInsertGrinder) {
  std::vector<Supernode*> snodes;
  snodes.reserve(num_nodes);
  for (unsigned i = 0; i < num_nodes; ++i)
    snodes[i] = Supernode::makeSupernode(seed);

  // insert all edges
  for (auto edge : *graph_edges) {
    const auto encoded = Supernode::concat_tuple_fn(edge);
    Update upd = {encoded, 1};
    for (int i = 1; i <= edge[0]; ++i) {
      snodes[edge[i]]->update(upd);
    }
  }

  // must have at least logn successes per supernode
  int successes = 0;

  Edge sampled[edge_conn + 1];
  SampleSketchRet retcode;
  for (unsigned i = 2; i < num_nodes; ++i) {
    for (int j = 0; j < (int) snodes[i]->get_num_sktch(); ++j) {
      snodes[i]->sample(sampled, &retcode);
      if (retcode == FAIL) continue;
    
      successes++;
      if (i >= num_nodes/2 && prime[i]) {
        ASSERT_EQ(retcode, ZERO) << "False positive in sample " << i;
      } else {
        ASSERT_NE(retcode, ZERO) << "False negative in sample " << i;

        // hardcode for 2-edges
        ASSERT_EQ(sampled[0], 2) << "Failed on edge with connectivity " << sampled[0];
        ASSERT_TRUE(sampled[2] % sampled[1] == 0) << "Failed on {" << sampled[1] << "," <<
                                                  sampled[2] << "} with i = " << i;
        ASSERT_TRUE(i == sampled[1] || i == sampled[2]) << "Failed on {" <<
                                                  sampled[1] << "," <<
                                                  sampled[2] << "} with i = " << i;
      }
    }
    ASSERT_GE(successes, (int) log2(num_nodes)) << "Fewer than logn successful queries: supernode " << i;
  }
  for (unsigned i = 0; i < num_nodes; ++i)
    delete snodes[i];
}

TEST_F(SupernodeTestSuite, TestSampleDeleteGrinder) {
  std::vector<Supernode*> snodes;
  snodes.reserve(num_nodes);
  for (unsigned i = 0; i < num_nodes; ++i)
    snodes[i] = Supernode::makeSupernode(seed);

  // insert all edges
  for (auto edge : *graph_edges) {
    const auto encoded = Supernode::concat_tuple_fn(edge);
    Update upd = {encoded, 1};
    for (int i = 1; i <= edge[0]; ++i) {
      snodes[edge[i]]->update(upd);
    }
  }
  // then remove half of them (odds)
  for (auto edge : *odd_graph_edges) {
    const auto encoded = Supernode::concat_tuple_fn(edge);
    Update upd = {encoded, -1};
    for (int i = 1; i <= edge[0]; ++i) {
      snodes[edge[i]]->update(upd);
    }
  }

  // must have at least logn successes per supernode
  int successes = 0;

  Edge sampled[edge_conn + 1];
  SampleSketchRet retcode;
  for (unsigned i = 2; i < num_nodes; ++i) {
    for (int j = 0; j < (int) snodes[i]->get_num_sktch(); ++j) {
      snodes[i]->sample(sampled, &retcode);
      if (retcode == FAIL) continue;
    
      successes++;
      if (i >= num_nodes/2 && i % 2) {
        ASSERT_EQ(retcode, ZERO) << "False positive in sample " << i;
      } else {
        ASSERT_NE(retcode, ZERO) << "False negative in sample " << i;
        // hardcode for 2-edges
        ASSERT_EQ(sampled[0], 2) << "Failed on edge with connectivity " << sampled[0];
        ASSERT_TRUE(sampled[2] % sampled[1] == 0) << "Failed on {" << sampled[1] << "," <<
                                                  sampled[2] << "} with i = " << i;
        ASSERT_TRUE((sampled[2] / sampled[1]) % 2 == 0) << "Failed on {" << sampled[1] << "," <<
                                                  sampled[2] << "} with i = " << i;
        ASSERT_TRUE(i == sampled[1] || i == sampled[2]) << "Failed on {" <<
                                                  sampled[1] << "," <<
                                                  sampled[2] << "} with i = " << i;
      }
    }
    ASSERT_GE(successes, (int) log2(num_nodes)) << "Fewer than logn successful queries: supernode " << i;
  }
  for (unsigned i = 0; i < num_nodes; ++i)
    delete snodes[i];
}

void inline apply_delta_to_node(Supernode* node, std::vector<Update>& updates) {
  auto* loc = (Supernode*) malloc(Supernode::get_size());
  Supernode::delta_supernode(node->seed, updates, loc);
  node->apply_delta_update(loc);
  free(loc);
}

TEST_F(SupernodeTestSuite, TestBatchUpdate) {
  unsigned long vec_size = 1000000000, num_updates = 100000;
  srand(time(nullptr));
  UpdateList ulist(vec_size, num_updates);
  auto seed = rand();
  Supernode::configure(vec_size, edge_conn);
  Supernode* supernode = Supernode::makeSupernode(seed);
  Supernode* supernode_batch = Supernode::makeSupernode(seed);
  auto start_time = std::chrono::steady_clock::now();
  for (auto& upd : ulist) {
    supernode->update(upd);
  }
  std::cout << "One by one updates took " << static_cast<std::chrono::duration<long double>>(std::chrono::steady_clock::now() - start_time).count() << std::endl;
  start_time = std::chrono::steady_clock::now();
  apply_delta_to_node(supernode_batch, ulist.updates());
  std::cout << "Batched updates took " << static_cast<std::chrono::duration<long double>>(std::chrono::steady_clock::now() - start_time).count() << std::endl;

  ASSERT_EQ(supernode->get_num_sktch(), supernode_batch->get_num_sktch());
  ASSERT_EQ(supernode->idx, supernode_batch->idx);
  for (int i=0;i<supernode->get_num_sktch();++i) {
    ASSERT_EQ(*supernode->get_sketch(i), *supernode_batch->get_sketch(i));
  }
}

TEST_F(SupernodeTestSuite, TestConcurrency) {
  int num_threads_per_group = 2;
  unsigned num_threads =
       std::thread::hardware_concurrency() / num_threads_per_group - 1; // hyperthreading?
  unsigned vec_len = 1000000;
  unsigned num_updates = 100000;
  Supernode::configure(vec_len, edge_conn);

  std::vector<UpdateList> test_list(num_threads, UpdateList(vec_len, num_updates));

  int seed = rand();
  Supernode *supernode = Supernode::makeSupernode(seed);
  Supernode *piecemeal = Supernode::makeSupernode(seed);

  GraphWorker::set_config(0,num_threads_per_group); // set number of threads per omp parallel

  // concurrently run batch_updates
  std::thread thd[num_threads];
  for (unsigned i = 0; i < num_threads; ++i) {
    thd[i] = std::thread(apply_delta_to_node, piecemeal, std::ref(test_list[i].updates()));
  }

  // do single-update sketch in the meantime
  for (auto& list : test_list) {
    for (auto upd : list) {
      supernode->update(upd);
    }
  }

  for (unsigned i = 0; i < num_threads; ++i) {
    thd[i].join();
  }

  for (int i = 0; i < supernode->get_num_sktch(); ++i) {
    ASSERT_EQ(*supernode->get_sketch(i), *piecemeal->get_sketch(i));
  }
}

TEST_F(SupernodeTestSuite, TestSerialization) {
  std::vector<Supernode*> snodes;
  snodes.reserve(num_nodes);
  for (unsigned i = 0; i < num_nodes; ++i)
    snodes[i] = Supernode::makeSupernode(seed);

  // insert all edges
  for (auto edge : *graph_edges) {
    const auto encoded = Supernode::concat_tuple_fn(edge);
    Update upd = {encoded, 1};
    for (int i = 1; i <= edge[0]; ++i) {
      snodes[edge[i]]->update(upd);
    }
  }

  std::cout << "Finished inserting" << std::endl;
  auto file = std::fstream("./out_supernode.txt", std::ios::out | std::ios::binary);
  snodes[num_nodes/2]->write_binary(file);
  file.close();
  std::cout << "Finished writing" << std::endl;

  auto in_file = std::fstream("./out_supernode.txt", std::ios::in | std::ios::binary);

  Supernode* reheated = Supernode::makeSupernode(seed, in_file);

  for (int i = 0; i < snodes[num_nodes / 2]->get_num_sktch(); ++i) {
    ASSERT_EQ(*snodes[num_nodes / 2]->get_sketch(i), *reheated->get_sketch(i));
  }
}
