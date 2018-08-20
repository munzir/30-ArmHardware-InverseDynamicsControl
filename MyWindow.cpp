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

#include "MyWindow.hpp"

#include <iostream>

//==============================================================================
MyWindow::MyWindow(Controller* _controller)
  : SimWindow(),
    mController(_controller),
    mCircleTask(false)
{
  assert(_controller != nullptr);

  // Set the initial target positon to the initial position of the end effector
  mTargetPosition = mController->getEndEffector()->getTransform().translation();
  Eigen::Vector3d xEE = mController->getEndEffector()->getTransform().rotation().block<3, 1>(0, 0);
  Eigen::Vector3d yEE = mController->getEndEffector()->getTransform().rotation().block<3, 1>(0, 1);
  Eigen::Vector3d zEE = mController->getEndEffector()->getTransform().rotation().block<3, 1>(0, 2);
  std::cout << "xEE: " << xEE(0) << ", " << xEE(1) << ", " << xEE(2) << std::endl;
  std::cout << "yEE: " << yEE(0) << ", " << yEE(1) << ", " << yEE(2) << std::endl;
  std::cout << "zEE: " << zEE(0) << ", " << zEE(1) << ", " << zEE(2) << std::endl;
  std::cout << "PEE: " << mTargetPosition(0) << ", " << mTargetPosition(1) << ", " << mTargetPosition(2) << std::endl;
  mTargetRPY.setZero();
  mTargetRPY = dart::math::matrixToEulerXYZ(mController->getEndEffector()->getTransform().rotation());
  std::cout << "RPY: " << mTargetRPY(0) << ", " << mTargetRPY(1) << ", " << mTargetRPY(2) << std::endl;
}

//==============================================================================
MyWindow::~MyWindow()
{
   delete mController;
   std::cout << "closing window ..." << std::endl;
}

//==============================================================================
void MyWindow::timeStepping()
{
  if (mCircleTask)
  {
    static double time = 0.0;
    const double dt = 0.0005;
    const double radius = 0.6;
    Eigen::Vector3d center = Eigen::Vector3d(0.0, 0.1, 0.0);

    mTargetPosition = center;
    mTargetPosition[0] = radius * std::sin(time);
    mTargetPosition[1] = 0.25 * radius * std::sin(time);
    mTargetPosition[2] = radius * std::cos(time);

    time += dt;
  }

  // Update the controller and apply control force to the robot
  mController->update(mTargetPosition, mTargetRPY);

  // Step forward the simulation
  //mWorld->step();
}

//==============================================================================
void MyWindow::drawWorld() const
{
  // Draw the target position
  if (mRI)
  {
    // Draw Target Frame
    Eigen::Matrix3d mat;
    mat =  Eigen::AngleAxisd(mTargetRPY(0), Eigen::Vector3d::UnitX()) *
           Eigen::AngleAxisd(mTargetRPY(1), Eigen::Vector3d::UnitY()) *
           Eigen::AngleAxisd(mTargetRPY(2), Eigen::Vector3d::UnitZ());
    Eigen::Vector3d localTranslation; 
    double axisLength = 0.1;
    double axisWidth = 0.003;
    
    localTranslation << axisLength/2, 0, 0;
    mRI->setPenColor(Eigen::Vector3d(0.8, 0.2, 0.2));
    mRI->pushMatrix();
    mRI->translate(mTargetPosition + mat*localTranslation);
    mRI->rotate(Eigen::Vector3d::UnitX(), mTargetRPY(0)*180/M_PI);
    mRI->rotate(Eigen::Vector3d::UnitY(), mTargetRPY(1)*180/M_PI);
    mRI->rotate(Eigen::Vector3d::UnitZ(), mTargetRPY(2)*180/M_PI);
    mRI->drawCube(Eigen::Vector3d(axisLength, axisWidth, axisWidth));
    mRI->popMatrix();

    localTranslation << 0, axisLength/2, 0;
    mRI->setPenColor(Eigen::Vector3d(0.2, 0.8, 0.2));
    mRI->pushMatrix();
    mRI->translate(mTargetPosition + mat*localTranslation);
    mRI->rotate(Eigen::Vector3d::UnitX(), mTargetRPY(0)*180/M_PI);
    mRI->rotate(Eigen::Vector3d::UnitY(), mTargetRPY(1)*180/M_PI);
    mRI->rotate(Eigen::Vector3d::UnitZ(), mTargetRPY(2)*180/M_PI);
    mRI->drawCube(Eigen::Vector3d(axisWidth, axisLength, axisWidth));
    mRI->popMatrix();

    localTranslation << 0, 0, axisLength/2;
    mRI->setPenColor(Eigen::Vector3d(0.2, 0.2, 0.8));
    mRI->pushMatrix();
    mRI->translate(mTargetPosition + mat*localTranslation);
    mRI->rotate(Eigen::Vector3d::UnitX(), mTargetRPY(0)*180/M_PI);
    mRI->rotate(Eigen::Vector3d::UnitY(), mTargetRPY(1)*180/M_PI);
    mRI->rotate(Eigen::Vector3d::UnitZ(), mTargetRPY(2)*180/M_PI);
    mRI->drawCube(Eigen::Vector3d(axisWidth, axisWidth, axisLength));
    mRI->popMatrix();

    // Draw End-Effector Frame
    Eigen::Vector3d eeRPY = dart::math::matrixToEulerXYZ(mController->getEndEffector()->getTransform().rotation());
    Eigen::Vector3d eePosition = mController->getEndEffector()->getTransform().translation();
    
    mat =  Eigen::AngleAxisd(eeRPY(0), Eigen::Vector3d::UnitX()) *
           Eigen::AngleAxisd(eeRPY(1), Eigen::Vector3d::UnitY()) *
           Eigen::AngleAxisd(eeRPY(2), Eigen::Vector3d::UnitZ());
    axisLength = 0.08;
    axisWidth = 0.006;
    
    localTranslation << axisLength/2, 0, 0;
    mRI->setPenColor(Eigen::Vector3d(0.8, 0.2, 0.2));
    mRI->pushMatrix();
    mRI->translate(eePosition + mat*localTranslation);
    mRI->rotate(Eigen::Vector3d::UnitX(), eeRPY(0)*180/M_PI);
    mRI->rotate(Eigen::Vector3d::UnitY(), eeRPY(1)*180/M_PI);
    mRI->rotate(Eigen::Vector3d::UnitZ(), eeRPY(2)*180/M_PI);
    mRI->drawCube(Eigen::Vector3d(axisLength, axisWidth, axisWidth));
    mRI->popMatrix();

    localTranslation << 0, axisLength/2, 0;
    mRI->setPenColor(Eigen::Vector3d(0.2, 0.8, 0.2));
    mRI->pushMatrix();
    mRI->translate(eePosition + mat*localTranslation);
    mRI->rotate(Eigen::Vector3d::UnitX(), eeRPY(0)*180/M_PI);
    mRI->rotate(Eigen::Vector3d::UnitY(), eeRPY(1)*180/M_PI);
    mRI->rotate(Eigen::Vector3d::UnitZ(), eeRPY(2)*180/M_PI);
    mRI->drawCube(Eigen::Vector3d(axisWidth, axisLength, axisWidth));
    mRI->popMatrix();

    localTranslation << 0, 0, axisLength/2;
    mRI->setPenColor(Eigen::Vector3d(0.2, 0.2, 0.8));
    mRI->pushMatrix();
    mRI->translate(eePosition + mat*localTranslation);
    mRI->rotate(Eigen::Vector3d::UnitX(), eeRPY(0)*180/M_PI);
    mRI->rotate(Eigen::Vector3d::UnitY(), eeRPY(1)*180/M_PI);
    mRI->rotate(Eigen::Vector3d::UnitZ(), eeRPY(2)*180/M_PI);
    mRI->drawCube(Eigen::Vector3d(axisWidth, axisWidth, axisLength));
    mRI->popMatrix();

  }

  // Draw world
  SimWindow::drawWorld();
}

//==============================================================================
void MyWindow::keyboard(unsigned char _key, int _x, int _y)
{
  double incremental = 0.01;
  double angleInc = 3*M_PI/180;
  
  switch (_key)
  {
    case 'c':  // print debug information
      if (mCircleTask)
      {
        std::cout << "Circle task [off]." << std::endl;
        mCircleTask = false;
      }
      else
      {
        std::cout << "Circle task [on]." << std::endl;
        mCircleTask = true;
      }
      break;
    case 'q':
      mTargetPosition[0] -= incremental;
      break;
    case 'w':
      mTargetPosition[0] += incremental;
      break;
    case 'a':
      mTargetPosition[1] -= incremental;
      break;
    case 's':
      mTargetPosition[1] += incremental;
      break;
    case 'z':
      mTargetPosition[2] -= incremental;
      break;
    case 'x':
      mTargetPosition[2] += incremental;
      break;
    case 'y':
      mTargetRPY[0] -= angleInc;
      std::cout << "RPY: " << mTargetRPY(0) << ", " << mTargetRPY(1) << ", " << mTargetRPY(2) << std::endl;
      break;
    case 'u':
      mTargetRPY[0] += angleInc;
      std::cout << "RPY: " << mTargetRPY(0) << ", " << mTargetRPY(1) << ", " << mTargetRPY(2) << std::endl;
      break;
    case 'h':
      mTargetRPY[1] -= angleInc;
      std::cout << "RPY: " << mTargetRPY(0) << ", " << mTargetRPY(1) << ", " << mTargetRPY(2) << std::endl;
      break;
    case 'j':
      mTargetRPY[1] += angleInc;
      std::cout << "RPY: " << mTargetRPY(0) << ", " << mTargetRPY(1) << ", " << mTargetRPY(2) << std::endl;
      break;
    case 'n':
      mTargetRPY[2] -= angleInc;
      std::cout << "RPY: " << mTargetRPY(0) << ", " << mTargetRPY(1) << ", " << mTargetRPY(2) << std::endl;
      break;
    case 'm':
      mTargetRPY[2] += angleInc;
      std::cout << "RPY: " << mTargetRPY(0) << ", " << mTargetRPY(1) << ", " << mTargetRPY(2) << std::endl;
      break;
    case '&':
      mController->mEnable = !mController->mEnable;
      cout << "Hardware Enabled: " << (mController->mEnable?"true":"false") << endl;
      break;
    default:
      // Default keyboard control
      SimWindow::keyboard(_key, _x, _y);
      break;
  }

  // Keyboard control for Controller
  mController->keyboard(_key, _x, _y);

  glutPostRedisplay();
}

