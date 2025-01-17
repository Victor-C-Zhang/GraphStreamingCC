#include <gtest/gtest.h>
#include "../include/util.h"

TEST(UtilTestSuite, TestNonDirectionalNonSEPairingFn) {
  std::pair<uint32_t,uint32_t> exp;
  for (int i = 0; i < 1000; ++i) {
    for (int j = 0; j < 1000; ++j) {
      if (i==j) continue;
      exp = {std::min(i,j),std::max(i,j)};
      ASSERT_EQ(exp, inv_nondir_non_self_edge_pairing_fn
      (nondirectional_non_self_edge_pairing_fn(i,j)));
    }
  }
}
