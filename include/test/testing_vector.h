#pragma once
#include <vector>

/**
 * This class takes as input N and P and generates a random stream of P
 * updates to a vector of length N. It then processes those updates to the
 * vector and makes both the updates and the processed vector visible.
 */
class Testing_Vector {
  const unsigned long vector_length, num_updates;
  std::vector<long> vect;
  std::vector<Update> stream;

public:
  //n is size of vector and m is number of updates
  Testing_Vector(unsigned long vector_length, unsigned long num_updates) :
      vector_length(vector_length), num_updates(num_updates),
      vect(vector_length, 0), stream(num_updates) {
    //Initialize the stream, and finalize the input vector.
    for (unsigned int i = 0; i < num_updates; i++){
      vec_t index = rand() % vector_length;
      long delta = rand() % 9 - 4;
      if (delta == 0) delta = 1;
      stream[i] = {index, delta};
      vect[index] += delta;
      //cout << "Index: " << stream[i].index << " Delta: " << stream[i].delta << endl;
    }
  }

  //get the ith update (starting from index 0)
  Update get_update(unsigned long i){
    return stream[i];
  }

  //get the ith entry of the processed Testing_Vector
  long get_entry(unsigned long i){
    return vect[i];
  }

  friend std::ostream& operator<< (std::ostream &os, const Testing_Vector &vec);
};

inline std::ostream& operator<<(std::ostream& os, const Testing_Vector& vec) {
  for (unsigned long i = 0; i < vec.num_updates; i++) {
    os << vec.stream[i] << std::endl;
  }
  for (const auto& val : vec.vect) {
    os << val << ',';
  }
  os << std::endl;
  return os;
}
