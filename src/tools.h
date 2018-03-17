#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

};

// Because auto is implemented with templating, these functions
// need to be defined in the header file.
void print(auto thing, bool doEndl=true) {
  std::cout << thing;
  if(doEndl) {
    std::cout << endl;
  }
}

void printt(auto thing1, auto thing2, bool doEndl=true) {
  print(thing1, false);
  print(thing2, doEndl);
}

#endif /* TOOLS_H_ */