#pragma once
#include <cstdlib>
#include <exception>
#include <set>
#include <fstream>
#include <atomic>  // REMOVE LATER

#include <guttering_system.h>
#include "supernode.h"

#ifdef VERIFY_SAMPLES_F
#include "test/graph_verifier.h"
#endif

// forward declarations
class GraphWorker;

typedef std::pair<Edge*, UpdateType> GraphUpdate;

// Exceptions the Graph class may throw
class UpdateLockedException : public std::exception {
  virtual const char* what() const throw() {
    return "The graph cannot be updated: Connected components algorithm has "
           "already started";
  }
};

class MultipleGraphsException : public std::exception {
  virtual const char * what() const throw() {
    return "Only one Graph may be open at one time. The other Graph must be deleted.";
  }
};

/**
 * Undirected graph object with n nodes labelled 0 to n-1, no self-edges,
 * multiple edges, or weights.
 */
class Graph {
protected:
  node_id_t num_nodes;
  uint64_t seed;
  bool update_locked = false;
  bool modified = false;
  // a set containing one "representative" from each supernode
  std::set<node_id_t>* representatives;
  Supernode** supernodes;
  // DSU representation of supernode relationship
  node_id_t* parent;
  node_id_t* size;
  node_id_t get_parent(node_id_t node);

  // Guttering system for batching updates
  GutteringSystem *gts;

  void backup_to_disk(const std::vector<node_id_t>& ids_to_backup);
  void restore_from_disk(const std::vector<node_id_t>& ids_to_restore);

  uint128_t concat_tuple_fn(const uint32_t* edge_buf) const;

  /**
   * Update the query array with new samples
   * @param edge_query      an array of supernode query results
   * @param retcode_query   an array of corresponding query return codes
   * @param reps            an array containing node indices for the
   *                        representative of each supernode
   */
  virtual void sample_supernodes(Edge *edge_query, SampleSketchRet *retcode_query,
                    std::vector<node_id_t> &reps);

  /**
   * @param copy_supernodes  an array to be filled with supernodes
   * @param to_merge         an list of lists of supernodes to be merged
   *
   */
  void merge_supernodes(Supernode** copy_supernodes, std::vector<node_id_t> &new_reps,
                        std::vector<std::vector<node_id_t>> &to_merge, bool make_copy);

  /**
   * Run the disjoint set union to determine what supernodes
   * Should be merged together.
   * Map from nodes to a vector of nodes to merge with them
   * @param edge_query      an array of supernode query results
   * @param retcode_query   an array of corresponding query retcodes
   * @param reps            an array containing node indices for the
   *                        representative of each supernode
   */
  std::vector<std::vector<node_id_t>> supernodes_to_merge(Edge *edge_query,
                      SampleSketchRet *retcode_query, std::vector<node_id_t> &reps);

  /**
   * Main parallel algorithm utilizing Boruvka and L_0 sampling.
   * @return a vector of the connected components in the graph.
   */
  std::vector<std::set<node_id_t>> boruvka_emulation(bool make_copy);


  std::string backup_file; // where to backup the supernodes
  bool copy_in_mem = false; // should backups be made in memory or on disk

  FRIEND_TEST(GraphTestSuite, TestCorrectnessOfReheating);
  FRIEND_TEST(GraphTest, TestSupernodeRestoreAfterCCFailure);

  static bool open_graph;
public:
  explicit Graph(node_id_t num_nodes, int edge_connectivity);
  explicit Graph(const std::string &input_file);

  ~Graph();

  /**
   * Update a graph with a hyperedge
   * @param upd a "hyperedge buffer":
   *                a buffer that contains the number of hyperedge vertices
   *                in index 0, and fills the remaining buffer values with
   *                vertices, starting at index 1.
   */
  void update(const GraphUpdate& upd);

  /**
   * Update all the sketches in supernode, given a batch of updates.
   * @param src        The supernode where the edges originate.
   * @param edges      A vector of destinations.
   * @param delta_loc  Memory location where we should initialize the delta
   *                   supernode.
   */
  void batch_update(node_id_t src, const std::vector<delta_update_t> &edges, Supernode *delta_loc);

  /**
   * Main parallel algorithm utilizing Boruvka and L_0 sampling.
   * If cont is true, allow for additional updates when done.
   * @param cont
   * @return a vector of the connected components in the graph.
   */
  std::vector<std::set<node_id_t>> connected_components(bool cont=false);

#ifdef VERIFY_SAMPLES_F
  std::unique_ptr<GraphVerifier> verifier;
  void set_verifier(std::unique_ptr<GraphVerifier> verifier) {
    this->verifier = std::move(verifier);
  }

  // to induce a failure mid-CC
  bool fail_round_2 = false;
  void should_fail_CC() { fail_round_2 = true; }
#endif

  // temp to verify number of updates -- REMOVE later
  std::atomic<uint64_t> num_updates;

  /**
   * Generate a delta node for the purposes of updating a node sketch
   * (supernode).
   * @param node_n     the total number of nodes in the graph.
   * @param node_seed  the seed of the supernode in question.
   * @param src        the src id.
   * @param edges      a list of edge indices to which src is connected and
   *                   their respective deltas.
   * @param delta_loc  the pre-allocated memory where the delta_node should be
   *                   placed. this allows memory to be reused by the same
   *                   calling thread.
   * @returns nothing (supernode delta is in delta_loc).
   */
  static void generate_delta_node(node_id_t node_n, uint64_t node_seed, node_id_t src,
                                  const std::vector<delta_update_t> &edges, Supernode
                                  *delta_loc);

  /**
   * Serialize the graph data to a binary file.
   * @param filename the name of the file to (over)write data to.
   */
  void write_binary(const std::string &filename);

  // time hooks for experiments
  std::chrono::steady_clock::time_point flush_start;
  std::chrono::steady_clock::time_point flush_end;
  std::chrono::steady_clock::time_point cc_alg_start;
  std::chrono::steady_clock::time_point cc_alg_end;
};
