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

#ifndef EXAMPLES_OPERATIONALSPACECONTROL_CONTROLLER_HPP_
#define EXAMPLES_OPERATIONALSPACECONTROL_CONTROLLER_HPP_

#include <Eigen/Eigen>

#include <dart/dart.hpp>

#include <unistd.h> //for usleep thing
#include <somatic.h>
#include <somatic/daemon.h>
#include <somatic.pb-c.h>
#include <somatic/motor.h>
#include <ach.h>

#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
#include <string>
#include <time.h>
#include <sys/time.h>

#include "file_ops.hpp"

#include <config4cpp/Configuration.h>

#include <nlopt.hpp>



#define VELOCITY SOMATIC__MOTOR_PARAM__MOTOR_VELOCITY
#define POSITION SOMATIC__MOTOR_PARAM__MOTOR_POSITION
#define CURRENT SOMATIC__MOTOR_PARAM__MOTOR_CURRENT


using namespace std;
using namespace config4cpp;

/// \brief Operational space controller for 6-dof manipulator
class Controller
{
public:
  /// \brief Constructor
  Controller(dart::dynamics::SkeletonPtr _robot,
             dart::dynamics::BodyNode* _endEffector);

  /// \brief Destructor
  virtual ~Controller();

  /// \brief
  void update(const Eigen::Vector3d& _targetPosition, const Eigen::Vector3d& _targetRPY);

  /// \brief initializes the structures related to the arm hardware
  void initArm (somatic_d_t& daemon_cx, somatic_motor_t& arm, const char* armName);

  /// \brief Get robot
  dart::dynamics::SkeletonPtr getRobot() const;

  /// \brief Get end effector of the robot
  dart::dynamics::BodyNode* getEndEffector() const;

  /// \brief Keyboard control
  virtual void keyboard(unsigned char _key, int _x, int _y);

  unsigned long get_time();

  bool mEnable;

private:
  /// \brief Robot
  dart::dynamics::SkeletonPtr mRobot;

  /// \brief End-effector of the robot
  dart::dynamics::BodyNode* mEndEffector;

  /// \brief Control forces
  Eigen::VectorXd mForces;

  /// \brief Proportional gain for the virtual spring forces at the end effector
  double mWPos, mWOr;
  Eigen::Matrix3d mKp, mKpOr;

  /// \brief Derivative gain for the virtual spring forces at the end effector
  Eigen::Matrix3d mKv, mKvOr;

  Eigen::Matrix<double, 7, 7> mWReg;
  double mKvReg;

  somatic_d_t mDaemon_cx;
  somatic_motor_t mSomaticSinglearm;  
  // double mqInit[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double mqInit[7] = {0.0, -1.04, 0.0, 1.57, 0.0, 0.0, 0.0};
  Eigen::Matrix<double, 7, 7> mRotorInertia, mViscousFriction, mCoulombFriction, mKmInv, mGRInv;
  double mKm_array[7] = {31.4e-3, 31.4e-3, 38e-3, 38e-3, 16e-3, 16e-3, 16e-3};
  double mGR_array[7] = {596, 596, 625, 625, 552, 552, 552};

  double mPriorTime;
  long mStartTime;
  
  Eigen::Matrix<double, 7, 1> mq, mdq;

  Eigen::Matrix<double, 7, 1> mCurLim;

  int mSteps;

  bool mdtFixed;
  double mdt;

  bool mCompensateFriction;
  bool mPredictFriction;

  double mGainFactor;
};

#endif  // EXAMPLES_OPERATIONALSPACECONTROL_CONTROLLER_HPP_
