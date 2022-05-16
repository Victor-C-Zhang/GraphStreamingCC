#include <stdexcept>
#include <cmath>
#include "../include/supernode.h"
#include "../include/graph_worker.h"
#include "types.h"

size_t Supernode::bytes_size;
int Supernode::edge_connectivity;
uint64_t Supernode::num_nodes;

Supernode::Supernode(uint64_t seed): idx(0), num_sketches(log2(num_nodes)/
              (log2(3)-1)), seed(seed), sketch_size(Sketch::sketchSizeof()) {

  size_t sketch_width = guess_gen(Sketch::get_failure_factor());
  // generate num_sketches sketches for each supernode (read: node)
  for (int i = 0; i < num_sketches; ++i) {
    Sketch::makeSketch(get_sketch(i), seed);
    seed += sketch_width;
  }
}

Supernode::Supernode(uint64_t seed, std::istream &binary_in) :
  idx(0), num_sketches(log2(num_nodes)/(log2(3)-1)), seed(seed), sketch_size
  (Sketch::sketchSizeof()) {

  size_t sketch_width = guess_gen(Sketch::get_failure_factor());
  // read num_sketches sketches from file for each supernode (read: node)
  for (int i = 0; i < num_sketches; ++i) {
    Sketch::makeSketch(get_sketch(i), seed, binary_in);
    seed += sketch_width;
  }
}

Supernode::Supernode(const Supernode& s) : idx(s.idx), num_sketches(s.num_sketches),
    seed(s.seed), sketch_size(s.sketch_size) {
  for (int i = 0; i < num_sketches; ++i) {
    Sketch::makeSketch(get_sketch(i), *s.get_sketch(i));
  }
}

Supernode* Supernode::makeSupernode(long seed, void *loc) {
  return new (loc) Supernode(seed);
}

Supernode* Supernode::makeSupernode(long seed, std::istream &binary_in, void *loc) {
  return new (loc) Supernode(seed, binary_in);
}

Supernode* Supernode::makeSupernode(const Supernode& s, void *loc) {
  return new (loc) Supernode(s);
}

Supernode::~Supernode() {
}


// pad {x_1, \dots, x_k} to {x_1, \dots, x_k, x_k, \dots, x_k}
// return (x_k, x_k, \dots, x_k, x_{k-1}, \dots, x_1)_{n}
uint128_t Supernode::concat_tuple_fn(const uint32_t* edge_buf) {
  const auto num_vals = edge_buf[0];
  uint128_t retval = 0;
  for (int i = 0; i < edge_connectivity - num_vals; ++i) {
    retval *= num_nodes;
    retval += edge_buf[num_vals];
  }
  for (int i = num_vals; i > 0; --i) {
    retval *= num_nodes;
    retval += edge_buf[i];
  }
  return retval;
}

inline void Supernode::inv_concat_tuple_fn(uint128_t catted, Edge* edge_buf) {
  edge_buf[1] = catted % num_nodes;
  catted /= num_nodes;
  int i = 2;
  while (catted > 0) {
    edge_buf[i] = catted % num_nodes;
    if (edge_buf[i] == edge_buf[i-1]) break;
    catted /= num_nodes;
    ++i;
  }
  edge_buf[0] = i - 1;
}


void Supernode::sample(Edge *edge_buf, SampleSketchRet *ret_buf) {
  if (idx == num_sketches) throw OutOfQueriesException();

  auto query_ret = get_sketch(idx++)->query();
  inv_concat_tuple_fn(query_ret.first, edge_buf);
  *ret_buf = query_ret.second;
}

void Supernode::merge(Supernode &other) {
  idx = std::max(idx, other.idx);
  for (int i=idx;i<num_sketches;++i) {
    (*get_sketch(i))+=(*other.get_sketch(i));
  }
}

void Supernode::update(Update& upd) {
  for (int i = 0; i < num_sketches; ++i)
    get_sketch(i)->update(upd);
}

void Supernode::apply_delta_update(const Supernode* delta_node) {
  std::unique_lock<std::mutex> lk(node_mt);
  for (int i = 0; i < num_sketches; ++i) {
    *get_sketch(i) += *delta_node->get_sketch(i);
  }
  lk.unlock();
}

/*
 * Consider fiddling with environment vars
 * OMP_DYNAMIC: whether the OS is allowed to dynamically change the number
 * of threads employed for each parallel section
 * OMP_NUM_THREADS (or set_omp_num_threads): how many threads to spin up for
 * each parallel section. the default is (probably) one per CPU core
 * available, but we may want to set it lower if num_sketches is a nice multiple of
 * a lower number.
 *
 * We may want to use omp option schedule(dynamic) or schedule(guided) if
 * there are very many more iterations of loop than threads. Dynamic
 * scheduling is good if loop iterations are expected to take very much
 * different amounts of time. Refer to
 * http://www.inf.ufsc.br/~bosco.sobral/ensino/ine5645/OpenMP_Dynamic_Scheduling.pdf
 * for a detailed explanation.
 */
/*
 * Current impl uses default threads and parallelism within batched_update.
 * Considered using spin-threads and parallelism within sketch::update, but
 * this was slow (at least on small graph inputs).
 */
void Supernode::delta_supernode(uint64_t seed, const std::vector<Update>
      &updates, void *loc) {
  auto delta_node = makeSupernode(seed, loc);
#pragma omp parallel for num_threads(GraphWorker::get_group_size()) default(shared)
  for (int i = 0; i < delta_node->num_sketches; ++i) {
    delta_node->get_sketch(i)->batch_update(updates);
  }
}

void Supernode::write_binary(std::ostream& binary_out) {
  for (int i = 0; i < num_sketches; ++i) {
    get_sketch(i)->write_binary(binary_out);
  }
}
