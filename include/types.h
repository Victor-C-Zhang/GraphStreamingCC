#pragma once
#include <xxhash.h>
#include <graph_zeppelin_common.h>

typedef uint64_t col_hash_t;
static const auto& vec_hash = XXH32;
static const auto& col_hash = XXH64;

enum UpdateType {
  INSERT = 0,
  DELETE = 1,
  END_OF_FILE = 2
}; // special type to indicate that there is no more data in the stream

typedef int128_t bucket_t;
typedef uint128_t ubucket_t;
typedef int128_t bucket_prod_t;
typedef uint128_t ubucket_prod_t;
