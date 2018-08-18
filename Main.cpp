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

#include <dart/dart.hpp>
#include <dart/utils/urdf/urdf.hpp>
#include "file_ops.hpp"

#include "MyWindow.hpp"

int main(int argc, char* argv[])
{
  
  // load skeletons
  dart::utils::DartLoader dl;
  dart::dynamics::SkeletonPtr ground
      = dl.parseSkeleton("dart://sample/urdf/KR5/ground.urdf");
  dart::dynamics::SkeletonPtr robot
      = dl.parseSkeleton("/home/munzir/Documents/Software/09-URDF/7DOFArm/singlearm.urdf");

  // load dynamic parameters and set them in the robot
  Eigen::MatrixXd beta 
      = readInputFileAsMatrix("../../20c-RidgeRegression_arm/betaFull/betaFull.txt");
  int numBodies = robot->getNumBodyNodes();
  int paramsPerBody = 13;
  for(int i=1; i<numBodies; i++) {
    int ind = paramsPerBody*(i-1);
    double m = beta(ind + 0);
    Eigen::Vector3d mCOM = beta.block(ind+1, 0, 3, 1);
    double xx = beta(ind + 4);
    double yy = beta(ind + 5);
    double zz = beta(ind + 6);
    double xy = beta(ind + 7);
    double xz = beta(ind + 8);
    double yz = beta(ind + 9);

    robot->getBodyNode(i)->setMass(m);
    robot->getBodyNode(i)->setLocalCOM(mCOM/m);
    robot->getBodyNode(i)->setMomentOfInertia(xx, yy, zz, xy, xz, yz);
  }

  for(int i=1; i<numBodies; i++) {
    cout << "body node: " << i << ". " << robot->getBodyNode(i)->getName() << endl;
    cout << "===================================" << endl;
    cout << "mass: " << robot->getBodyNode(i)->getMass() << endl;
    cout << "Local mCOM: " << robot->getBodyNode(i)->getLocalCOM().transpose()*robot->getBodyNode(i)->getMass() << endl;
    double xx, yy, zz, xy, yz, xz;
    robot->getBodyNode(i)->getMomentOfInertia(xx, yy, zz, xy, xz, yz);
    cout << "xx: " << xx << endl;
    cout << "yy: " << yy << endl;
    cout << "zz: " << zz << endl;
    cout << "xy: " << xy << endl;
    cout << "xz: " << xz << endl;
    cout << "yz: " << yz << endl << endl;
  }
  
  // create and initialize the world
  dart::simulation::WorldPtr world(new dart::simulation::World);
  assert(world != nullptr);
  world->addSkeleton(ground);
  world->addSkeleton(robot);

  // create and initialize the world
  Eigen::Vector3d gravity(0.0, -9.81, 0.0);
  world->setGravity(gravity);
  world->setTimeStep(1.0/1000);

  // create a window and link it to the world
  MyWindow window(new Controller(robot, robot->getBodyNode("L7_EFF")));
  window.setWorld(world);

  glutInit(&argc, argv);
  window.initWindow(640, 480, "Forward Simulation");
  glutMainLoop();

  return 0;
}
