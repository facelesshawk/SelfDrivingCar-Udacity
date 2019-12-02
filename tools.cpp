#include "tools.h"
#include <iostream>
#include<math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd res,resi;
  for(int i=0;i<estimations.size();i++)
  {
    resi = (estimations[i]-ground_truth[i]);
    resi = resi.array()*resi.array();
    res += resi;
  }
  	res = (res/estimations.size());
  res = res.array().sqrt();
  return res;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  double px = x_state[0];
  double py = x_state[1];
  double vx = x_state[2];
  double vy = x_state[3];
  double res = sqrt((px*px)+(py*py));
  double var1 = (py*((vx*py)-(vy*px)));
  double  var2 = (px*((vy*px)-(vx*py)));
  MatrixXd hjacob(3,4);
  hjacob << (px/res),(py/res),0,0,-(py/pow(res,2)),(px/pow(res,2)),0,0,-(var1/pow(res,1.5)),(var2/pow(res,1.5)),(px/res),(py/res);
  return hjacob;
}
