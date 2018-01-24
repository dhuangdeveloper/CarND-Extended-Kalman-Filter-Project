#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

//Tools::Tools() {}

//Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // Calculate the RMSE here
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  // check the validity of the following inputs:
  // * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if ( (estimations.size()==0) || (ground_truth.size()==0) || (ground_truth.size()!= estimations.size()))
  {
    return rmse;
  }
  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
	  // ... your code here
	  VectorXd residual(4);
	  residual = estimations[i] - ground_truth[i];
	  residual = residual.array()*residual.array(); 
	  rmse += residual;
  }
  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {  
  //Calculate a Jacobian here.
  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
    
  float pxy2 = px * px + py * py;
  float pxy = sqrt(pxy2);
  float pxy3_2 = pxy * pxy2;
  //check division by zero 
  
  if ( (fabs(px) < 0.0001  && fabs(py) < 0.0001))
  {
	  cout << "Error: Both px and py are zero " << endl;
	  return Hj;
  }
  
  //compute the Jacobian matrix
  Hj << px / pxy, py / pxy, 0, 0,
  - py / pxy2, px / pxy2, 0, 0,
  py * (vx * py - vy * px) / pxy3_2, px * (vy *px - vx * py) / pxy3_2, px / pxy, py / pxy;
  
  return Hj;  
}
