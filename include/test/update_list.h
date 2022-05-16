#pragma once
#include <vector>
#include "../l0_sampling/update.h"

/**
 * This class takes as input N and P and generates a random stream of P
 * updates to a vector of length N.
 */
class UpdateList {
  const unsigned long num_updates;
  std::vector<Update> stream;

public:
  //n is size of vector and m is number of updates
  UpdateList(unsigned long vector_length, unsigned long num_updates) :
        num_updates(num_updates), stream(num_updates) {
    //Initialize the stream, and finalize the input vector.
    for (unsigned int i = 0; i < num_updates; i++){
      vec_t index = rand() % vector_length;
      long delta = rand() % 9 - 4;
      if (delta == 0) delta = 1;
      stream[i] = {index, delta};
    }
  }

  std::vector<Update>& updates(){
    return stream;
  }

  std::vector<Update>::iterator begin() {
    return stream.begin();
  }

  std::vector<Update>::iterator end() {
    return stream.end();
  }

  friend std::ostream& operator<< (std::ostream &os, const UpdateList &vec);
};

inline std::ostream& operator<<(std::ostream& os, const UpdateList& vec) {
  for (unsigned long i = 0; i < vec.num_updates; i++) {
    os << vec.stream[i] << std::endl;
  }
  os << std::endl;
  return os;
}
