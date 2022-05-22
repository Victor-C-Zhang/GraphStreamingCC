#include "../../include/l0_sampling/sketch.h"
#include "prime_generator.h"
#include <cassert>
#include <cstring>
#include <iostream>

vec_t Sketch::failure_factor = 100;
vec_t Sketch::n;
size_t Sketch::num_elems;
size_t Sketch::num_buckets;
size_t Sketch::num_guesses;
ubucket_t Sketch::large_prime;

/*
 * Static functions for creating sketches with a provided memory location.
 * We use these in the production system to keep supernodes virtually contiguous.
 */
Sketch* Sketch::makeSketch(void* loc, uint64_t seed) {
  return new (loc) Sketch(seed);
}

Sketch* Sketch::makeSketch(void* loc, uint64_t seed, std::istream &binary_in) {
  return new (loc) Sketch(seed, binary_in);
}

Sketch* Sketch::makeSketch(void* loc, const Sketch& s) {
  return new (loc) Sketch(s);
}

Sketch::Sketch(uint64_t seed): seed(seed) {
  // zero buckets
  std::memset(_bucket_data, 0, num_elems * sizeof(Bucket_Boruvka));

  buckets = reinterpret_cast<Bucket_Boruvka *>(_bucket_data);
}

Sketch::Sketch(uint64_t seed, std::istream &binary_in): seed(seed) {
  binary_in.read(_bucket_data, num_elems * sizeof(Bucket_Boruvka));
  buckets = reinterpret_cast<Bucket_Boruvka *>(_bucket_data);
}

Sketch::Sketch(const Sketch& s) : seed(s.seed) {
  std::memcpy(_bucket_data, s._bucket_data, num_elems * sizeof(Bucket_Boruvka));
  buckets = reinterpret_cast<Bucket_Boruvka *>(_bucket_data);
}

void Sketch::update(const Update& update) {
  auto cbucket_seed = Bucket_Boruvka::gen_bucket_seed(num_elems - 1, seed);
  auto cr = Bucket_Boruvka::gen_r(cbucket_seed, large_prime);
  buckets[num_elems - 1].update(update, large_prime, cr);

  for (unsigned i = 0; i < num_buckets; ++i) {
    for (unsigned j = 0; j < num_guesses; ++j) {
      unsigned bucket_id = i * num_guesses + j;
      auto bucket_seed = Bucket_Boruvka::gen_bucket_seed(bucket_id, seed);
      auto r = Bucket_Boruvka::gen_r(bucket_seed, large_prime);
      auto& bucket = buckets[bucket_id];
      if (bucket.contains(update.index, bucket_seed, 2 << j)) {
        bucket.update(update, large_prime, r);
      }
    }
  }
}

void Sketch::batch_update(const std::vector<Update>& updates) {
  for (const auto& upd : updates) {
    update(upd);
  }
}

std::pair<bucket_t , SampleSketchRet> Sketch::query() {
  if (already_queried) {
    throw MultipleQueryException();
  }
  already_queried = true;

  auto& determ_bucket = buckets[num_elems - 1];
  if (determ_bucket == BUCKET_ZERO) {
    return {0, ZERO}; // the "first" bucket is deterministic so if it is all zero then there are no edges to return
  }

  auto cbucket_seed = Bucket_Boruvka::gen_bucket_seed(num_elems - 1, seed);
  auto cr = Bucket_Boruvka::gen_r(cbucket_seed, large_prime);
  if (determ_bucket.is_good(n, large_prime, cbucket_seed, cr)) {
    return {determ_bucket.b / determ_bucket.a - 1, GOOD};
  }
  for (unsigned i = 0; i < num_buckets; ++i) {
    for (unsigned j = 0; j < num_guesses; ++j) {
      unsigned bucket_id = i * num_guesses + j;
      auto bucket_seed = Bucket_Boruvka::gen_bucket_seed(bucket_id, seed);
      auto r = Bucket_Boruvka::gen_r(bucket_seed, large_prime);
      auto& bucket = buckets[bucket_id];
      if (bucket.is_good(n, large_prime, bucket_seed, r, 2 << j)) {
        return {bucket.b / bucket.a - 1, GOOD};
      }
    }
  }
  return {0, FAIL};
}

Sketch &operator+= (Sketch &sketch1, const Sketch &sketch2) {
  assert (sketch1.seed == sketch2.seed);
  assert (sketch1.large_prime == sketch2.large_prime);
  for (unsigned i = 0; i < Sketch::num_elems; i++) {
    auto& bucket1 = sketch1.buckets[i];
    const auto& bucket2 = sketch2.buckets[i];
    bucket1.a += bucket2.a;
    bucket1.b += bucket2.b;
    bucket1.c += bucket2.c;
    bucket1.c %= Sketch::large_prime;
  }
  sketch1.already_queried = sketch1.already_queried || sketch2.already_queried;
  return sketch1;
}

bool operator== (const Sketch &sketch1, const Sketch &sketch2) {
  if (sketch1.seed != sketch2.seed || sketch1.already_queried != sketch2.already_queried) 
    return false;

  for (size_t i = 0; i < Sketch::num_elems; ++i) {
    if (sketch1.buckets[i] != sketch2.buckets[i]) return false;
  }

  return true;
}

std::ostream& operator<< (std::ostream &os, const Sketch &sketch) {
  unsigned cbucket_id = Sketch::num_buckets * Sketch::num_guesses;
  auto cbucket_seed = Bucket_Boruvka::gen_bucket_seed(cbucket_id, sketch.seed);
  auto cr = Bucket_Boruvka::gen_r(cbucket_seed, Sketch::large_prime);
  auto& cbucket = sketch.buckets[cbucket_id];
  for (unsigned k = 0; k < Sketch::n; k++) {
    os << '1';
  }
  os << std::endl
     << "a:" << std::hex << (uint64_t)(cbucket.a >> 64) << (uint64_t) cbucket.a << std::endl
     << "c:" << std::hex << (uint64_t)(cbucket.a >> 64) << (uint64_t) cbucket.a << std::endl
     << (cbucket.is_good(1, Sketch::large_prime, cbucket_seed, cr, 1) ? "good" : "bad") << std::endl;

  for (unsigned i = 0; i < Sketch::num_buckets; ++i) {
    for (unsigned j = 0; j < Sketch::num_guesses; ++j) {
      unsigned bucket_id = i * Sketch::num_guesses + j;
      auto bucket_seed = Bucket_Boruvka::gen_bucket_seed(bucket_id, sketch.seed);
      auto r = Bucket_Boruvka::gen_r(bucket_seed, Sketch::large_prime);
      auto& bucket = sketch.buckets[bucket_id];

      for (unsigned k = 0; k < Sketch::n; k++) {
        os << (bucket.contains(k, bucket_seed, 1<<j) ? '1' : '0');
      }
      os << std::endl
         << "a:" << std::hex << (uint64_t)(bucket.a >> 64) << (uint64_t) bucket.a << std::endl
         << "c:" << std::hex << (uint64_t)(bucket.a >> 64) << (uint64_t) bucket.a << std::endl
         << (bucket.is_good(i, Sketch::large_prime, bucket_seed, r, 1 << j) ? "good" : "bad") << std::endl;
    }
  }
  return os;
}

void Sketch::write_binary(std::ostream& binary_out) {
  const_cast<const Sketch*>(this)->write_binary(binary_out);
}

void Sketch::write_binary(std::ostream &binary_out) const {
  binary_out.write(_bucket_data, num_elems * sizeof(Bucket_Boruvka));
}
