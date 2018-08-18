/*
 * Copyright (c) 2014-2016, Humanoid Lab, Georgia Tech Research Corporation
 * Copyright (c) 2014-2017, Graphics Lab, Georgia Tech Research Corporation
 * Copyright (c) 2016-2017, Personal Robotics Lab, Carnegie Mellon University
 * All rights reserved.
 *
 * This file is provided under the following "BSD-style" License:
 *   Redistribution and use in source and binary forms, with or
 *   without modification, are permitted provided that the following
 *   conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 *   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *   POSSIBILITY OF SUCH DAMAGE.
 */

#include "Controller.hpp"
#include <nlopt.hpp>
//==============================================================================
Controller::Controller(dart::dynamics::SkeletonPtr _robot,
                       dart::dynamics::BodyNode* _endEffector)
  : mRobot(_robot),
    mEndEffector(_endEffector)
{
  assert(_robot != nullptr);
  assert(_endEffector != nullptr);

  int dof = mRobot->getNumDofs();

  mForces.setZero(dof);

  mKp.setZero();
  mKv.setZero();

  for (int i = 0; i < 3; ++i)
  {
    mKp(i, i) = 750.0;
    mKv(i, i) = 250.0;

    mKpOr(i, i) = 150.0;
    mKvOr(i, i) = 50; 
  }

  // Remove position limits
  for (int i = 0; i < dof; ++i)
    _robot->getJoint(i)->setPositionLimitEnforced(false);

  // Set joint damping
  for (int i = 0; i < dof; ++i)
    _robot->getJoint(i)->setDampingCoefficient(0, 0.5);
}

//==============================================================================
Controller::~Controller()
{
}
//==============================================================================
struct OptParams{
  Eigen::Matrix<double, -1, 7> P;
  Eigen::VectorXd b;
};



//==============================================================================
void matprint(Eigen::MatrixXd A){
  for(int i=0; i<A.rows(); i++){
    for(int j=0; j<A.cols(); j++){
      std::cout << A(i,j) << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

//==============================================================================
double optFunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
  OptParams* optParams = reinterpret_cast<OptParams *>(my_func_data);
  Eigen::Matrix<double, 7, 1> X(x.data());

  if (!grad.empty()) {
    Eigen::Matrix<double, 7, 1> mGrad = optParams->P.transpose()*(optParams->P*X - optParams->b);
    Eigen::VectorXd::Map(&grad[0], mGrad.size()) = mGrad;
  }
  return (0.5 * pow((optParams->P*X - optParams->b).norm(), 2));
}

//==============================================================================
void Controller::update(const Eigen::Vector3d& _targetPosition, const Eigen::Vector3d& _targetRPY)
{
  double wPos = 1, wOr = 1;

  using namespace dart;

  Eigen::VectorXd dq = mRobot->getVelocities();                 // n x 1

  // End-effector Position
  Eigen::Vector3d x = mEndEffector->getTransform().translation();
  Eigen::Vector3d dx = mEndEffector->getLinearVelocity();
  Eigen::Vector3d ddxref = -mKp*(x - _targetPosition) - mKv*dx;
  math::LinearJacobian Jv = mEndEffector->getLinearJacobian();       // 3 x n
  math::LinearJacobian dJv = mEndEffector->getLinearJacobianDeriv();  // 3 x n
  Eigen::Matrix<double, 3, 7> PPos = Jv;
  Eigen::Vector3d bPos = -(dJv*dq - ddxref);

  // End-effector Orientation
  Eigen::Quaterniond quat(mEndEffector->getTransform().rotation());
  double quat_w = quat.w(); 
  Eigen::Vector3d quat_xyz(quat.x(), quat.y(), quat.z());
  if(quat_w < 0) {quat_w *= -1.0; quat_xyz *= -1.0; }
  Eigen::Quaterniond quatRef(Eigen::AngleAxisd(_targetRPY(0), Eigen::Vector3d::UnitX()) *
           Eigen::AngleAxisd(_targetRPY(1), Eigen::Vector3d::UnitY()) *
           Eigen::AngleAxisd(_targetRPY(2), Eigen::Vector3d::UnitZ()));
  double quatRef_w = quatRef.w(); 
  Eigen::Vector3d quatRef_xyz(quatRef.x(), quatRef.y(), quatRef.z());
  if(quatRef_w < 0) { quatRef_w *= -1.0; quatRef_xyz *= -1.0; }
  Eigen::Vector3d quatError_xyz = quatRef_w*quat_xyz - quat_w*quatRef_xyz + quatRef_xyz.cross(quat_xyz);
  // double quatError_w = quat_w*quatRef_w - quat_xyz.dot(quatRef_xyz);
  Eigen::Vector3d w = mEndEffector->getAngularVelocity();
  // Eigen::Vector3d dwref = -mKpOr*quatError_xyz/(2*quatError_w) - mKvOr*w;
  Eigen::Vector3d dwref = -mKpOr*quatError_xyz - mKvOr*w;
  math::AngularJacobian Jw = mEndEffector->getAngularJacobian();       // 3 x n
  math::AngularJacobian dJw = mEndEffector->getAngularJacobianDeriv();  // 3 x n
  Eigen::Matrix<double, 3, 7> POr = Jw;
  Eigen::Vector3d bOr = -(dJw*dq - dwref);


  // Optimizer stuff
  nlopt::opt opt(nlopt::LD_MMA, 7);
  OptParams optParams;
  std::vector<double> ddq_vec(7);
  double minf;

  // Perform optimization to find joint accelerations 
  Eigen::MatrixXd P(PPos.rows() + POr.rows(), PPos.cols() );
  P << wPos*PPos,
       wOr*POr;
  
  Eigen::VectorXd b(bPos.rows() + bOr.rows(), bPos.cols() );
  b << wPos*bPos,
       wOr*bOr;
       
  optParams.P = P;
  optParams.b = b;
  opt.set_min_objective(optFunc, &optParams);
  opt.set_xtol_rel(1e-4);
  opt.set_maxtime(0.005);
  opt.optimize(ddq_vec, minf);
  Eigen::Matrix<double, 7, 1> ddq(ddq_vec.data()); 
  

  //torques
  Eigen::MatrixXd M = mRobot->getMassMatrix();                   // n x n
  Eigen::VectorXd Cg   = mRobot->getCoriolisAndGravityForces();        // n x 1
  mForces = M*ddq + Cg;

  // Apply the joint space forces to the robot
  mRobot->setForces(mForces);
}

//==============================================================================
dart::dynamics::SkeletonPtr Controller::getRobot() const
{
  return mRobot;
}

//==============================================================================
dart::dynamics::BodyNode* Controller::getEndEffector() const
{
  return mEndEffector;
}

//==============================================================================
void Controller::keyboard(unsigned char /*_key*/, int /*_x*/, int /*_y*/)
{
}

