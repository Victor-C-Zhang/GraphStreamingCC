#include <fstream>
#include <algorithm>
#include <random>
#include <iostream>
#include <cassert>
#include "../../include/test/graph_gen.h"
#include "../../include/graph.h"

#define endl '\n'

typedef uint32_t ul;
typedef uint64_t ull;

const ull ULLMAX = std::numeric_limits<ul>::max();


std::ofstream& operator<< (std::ofstream &os, const std::pair<ull,ull> p) {
  os << p.first << " " << p.second;
  return os;
}

int conn;
node_id_t n;

// pad {x_1, \dots, x_k} to {x_1, \dots, x_k, x_k, \dots, x_k}
// return (x_k, x_k, \dots, x_k, x_{k-1}, \dots, x_1)_{n}
uint128_t concat_tuple_fn(const uint32_t* edge_buf) {
  const auto num_vals = edge_buf[0];
  uint128_t retval = 0;
  for (int i = 0; i < conn - num_vals; ++i) {
    retval *= n;
    retval += edge_buf[num_vals];
  }
  for (int i = num_vals; i > 0; --i) {
    retval *= n;
    retval += edge_buf[i];
  }
  return retval;
}

void inv_concat_tuple_fn(uint128_t catted, Edge* edge_buf) {
  edge_buf[1] = catted % n;
  catted /= n;
  int i = 2;
  while (catted > 0) {
    edge_buf[i] = catted % n;
    if (edge_buf[i] == edge_buf[i-1]) break;
    catted /= n;
    ++i;
  }
  edge_buf[0] = i - 1;
}

// simulate k-nested for loop with the inner loop doing something
void loop(unsigned idx, std::vector<uint32_t>& state, uint128_t* arr, ul& e) {
  if (idx == state[0]) {
    for (auto i = state[idx - 1] + 1; i < n; ++i) {
      state[idx] = i;
      arr[e++] = concat_tuple_fn(state.data());
    }
    return;
  } else if (idx == 1) {
    for (int i = 0; i < n; ++i) {
      state[idx] = i;
      loop(idx + 1, state, arr, e);
    }
  } else {
    for (auto i = state[idx - 1] + 1; i < n; ++i) {
      state[idx] = i;
      loop(idx + 1, state, arr, e);
    }
  }
}

void write_edges(double p, const std::string& out_f) {
  // really bad calculation
  ul num_edges = 0;
  for (int i = 2; i <= conn; ++i) {
    ull temp = n;
    for (int j = 1; j < i; ++j) {
      temp *= n - j;
    }
    for (int j = 2; j <= i; ++j) {
      temp /= j;
    }
    num_edges += temp;
  }

  ul num_idxs = pow(n,conn);
  auto* arr = static_cast<uint128_t *>(malloc(num_idxs * sizeof(uint128_t)));
  std::fill(arr, arr + num_idxs, 0);
  ul e = 0;

  std::vector<node_id_t> state(conn + 1);
  // TODO: do all edges with non-max connectivity too
  // do only edges with max connectivity
  state[0] = conn;
  loop(1, state, arr, e);

  std::shuffle(arr,arr+num_idxs, std::mt19937(std::random_device()()));
  std::ofstream out(out_f);
  ul desired_m = (ul) (num_edges*p);
  ul m = 0;
  for (int i = 0; i < desired_m; ++m) {
    if (arr[m] != 0) ++i;
  }
  out << n << " " << desired_m << endl;
  Edge buf[conn + 1];

  while (m--) {
    if (arr[m] == 0) continue; // may have some zeroes for unaccessed indices
    inv_concat_tuple_fn(arr[m], buf);
    out << buf[0];
    for (int i = 1; i <= buf[0]; ++i) {
      out << "\t" << buf[i];
    }
    out << "\n";
  }
  out.flush();
  out.close();
  free(arr);
}

void insert_delete(double p, int max_appearances, const std::string& in_file,
                   const std::string& out_file, const std::string&
                   cumul_out_file) {
  std::ifstream in(in_file);
  std::ofstream out(out_file);
  std::ofstream cumul_out(cumul_out_file);
  ul m; in >> m >> m;
  long long full_m = m;
  ull ins_del_arr[(ul)log2(m)+2];
  std::fill(ins_del_arr,ins_del_arr + (ul)log2(m)+2,0);
  ins_del_arr[0] = m;
  if (max_appearances == 0) {
    for (unsigned i = 0; ins_del_arr[i] > 1; ++i) {
      ins_del_arr[i + 1] = (ull) (ins_del_arr[i] * p);
      full_m += ins_del_arr[i + 1];
    }
  } else {
    for (int i = 0; i < max_appearances - 1; ++i) {
      ins_del_arr[i + 1] = (ull) (ins_del_arr[i] * p);
      full_m += ins_del_arr[i + 1];
    }
  }
  unsigned stopping = 1;
  if (max_appearances == 0) {
    for (; ins_del_arr[stopping] >= 1; ++stopping);
  } else {
    stopping = max_appearances;
  }
  ull cumul_m = m;
  for (int i = 1; i < stopping; ++i) {
    auto mult = (i % 2 == 1) ? -1 : 1;
    cumul_m += mult * ins_del_arr[i];
  }

  out << n << " " << full_m << endl;
  cumul_out << n << " " << cumul_m << endl;

  auto* memoized = static_cast<uint128_t *>(malloc( ins_del_arr[1] * sizeof(uint128_t)));
  Edge buf[conn + 1];

  // write all only-one-insertion edges
  for (unsigned i=ins_del_arr[1];i<m;++i) {
    in >> buf[0];
    for (int j = 1; j <= buf[0]; ++j) {
      in >> buf[j];
    }
    out << ((buf[0] << 1) | 0);
    for (int j = 1; j <= buf[0]; ++j) {
      out << "\t" << buf[j];
    }
    out << "\n";

    cumul_out << buf[0];
    for (int j = 1; j <= buf[0]; ++j) {
      cumul_out << "\t" << buf[j];
    }
    cumul_out << "\n";
  }

  // write other edges to output
  for (unsigned i=0;i<ins_del_arr[1];++i) {
    in >> buf[0];
    for (int j = 1; j <= buf[0]; ++j) {
      in >> buf[j];
    }
    out << ((buf[0] << 1) | 0);
    for (int j = 1; j <= buf[0]; ++j) {
      out << "\t" << buf[j];
    }
    out << "\n";
    memoized[i] = concat_tuple_fn(buf);
  }

  in.close();

  // write ins/out to output
  for (unsigned i = 1; i < stopping; ++i) {
    int temp = i % 2;
    for (unsigned j = 0; j < ins_del_arr[i]; ++j) {
      inv_concat_tuple_fn(memoized[j], buf);
      out << ((buf[0] << 1) | temp);
      for (int k = 1; k <= buf[0]; ++k) {
        out << "\t" << buf[k];
      }
      out << "\n";
    }
  }
  out.flush();
  out.close();

  // write rest of cumul output
  for (int i = 2; i < stopping; i+=2) {
    ull to_stop;
    if (i + 1 == stopping) {
      to_stop = 0;
    } else {
      to_stop = ins_del_arr[i + 1];
    }
    for (ull j = ins_del_arr[i] - 1ull; j >= to_stop; --j) {
      inv_concat_tuple_fn(memoized[j], buf);
      cumul_out << buf[0];
      for (int k = 1; k <= buf[0]; ++k) {
        cumul_out << "\t" << buf[k];
      }
      cumul_out << "\n";
      if (j == 0) break; // since we can't decrement to -1
    }
  }
  cumul_out.flush();
  cumul_out.close();

  free(memoized);
}

void write_cumul(const std::string& stream_f, const std::string& cumul_f) {
  std::ifstream in(stream_f);
  std::ofstream out(cumul_f);
  ull m; in >> m >> m;
  std::vector<std::vector<bool>> adj(n,std::vector<bool>(pow(n,conn),false));
  Edge buf[conn + 1];
  bool parity;
  for (ull i=1;i<=m;++i) {
    in >> buf[0];
    parity = buf[0] & 1;
    buf[0] >>= 1;
    for (int j = 1; j <= buf[0]; ++j) {
      in >> buf[j];
    }
    auto a = buf[1];
    auto b = concat_tuple_fn(buf);
    if ((parity == INSERT && adj[a][b] == 1) || (parity == DELETE &&
    adj[a][b] == 0)) {
      std::cerr << "Insertion/deletion error at line " << i
                << " in " << stream_f;
      return;
    }
    adj[a][b] = !adj[a][b];
  }

  in.close();

  // write cumul output
  ull m_cumul = 0;
  for (node_id_t i = 0; i < n; ++i) {
    for (uint128_t j = 0; j < pow(n,conn); ++j) {
      if (adj[i][j]) ++m_cumul;
    }
  }
  out << n << " " << m_cumul << endl;
  for (int i = 0; i < n; ++i) {
    for (uint128_t j = 0; j < pow(n,conn); ++j) {
      if (adj[i][j]) {
        inv_concat_tuple_fn(j, buf);
        out << buf[0];
        for (int k = 1; k <= buf[0]; ++k) {
          out << "\t" << buf[k];
        }
        out << "\n";
      }
    }
  }
  out.flush();
  out.close();
}

void generate_stream(const GraphGenSettings& settings) {
  n = settings.n;
  conn = settings.conn;
  write_edges(settings.p, "./TEMP_F");
  insert_delete(settings.r, settings.max_appearances, "./TEMP_F", settings
  .out_file, settings.cumul_out_file);
//  write_cumul(settings.out_file,settings.cumul_out_file);
}
