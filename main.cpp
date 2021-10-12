// CS591C1 Spring 2021 Final Project
// A "zero-twist" rod visualizer:
// Given a spline or helix, the program will show what the curve containing "zero-twist" will look like

// CONTRIBUTORS: DANIEL KEHR, PAUL MENEXAS
// CODE CITATION:
//     libigl's supplied example project and visualization tutorials aided in programming
//	the user-definable parameters of the libigl viewing window: https://libigl.github.io/tutorial/

#include "Eigen/Eigen"
#include <unsupported/Eigen/Splines>
#include <igl/floor.h>
#include <igl/readOFF.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui_demo.cpp>
#include <imgui/imgui.h>
#include <iostream>
#include <stdio.h>
#include <iostream>
#include <math.h>

// define spline type
typedef Eigen::Spline<double, 3> Spline3d;

using namespace Eigen;

// Eigen does not support a rowwise crossproduct function for a matrix, so I have implemented
// my own version here -DK
MatrixXd crossProduct(MatrixXd m1, MatrixXd m2)
{
  MatrixXd cross = MatrixXd::Constant(m1.rows(),3,0);
  cross(all,0) = (m1(all,1)*m2(all,2).transpose()-m1(all,2)*m2(all,1).transpose()).diagonal();
  cross(all,1) = (m1(all,2)*m2(all,0).transpose()-m1(all,0)*m2(all,2).transpose()).diagonal();
  cross(all,2) = (m1(all,0)*m2(all,1).transpose()-m1(all,1)*m2(all,0).transpose()).diagonal();
  return cross;
}

// Same as above, I have made a quick rowwise dot product function -DK
MatrixXd rowDot(MatrixXd m1, MatrixXd m2)
{
  return (m1*m2.transpose()).diagonal();
}

// rowwise norm calculation
MatrixXd rowNorm(MatrixXd m)
{
  return rowDot(m,m).array().sqrt();
}

// straighten line function for old approach (before spline). This is not used in the final
// product but still helps for debugging purposes -DK
MatrixXd straighten(MatrixXd vertices, double edgeLength)
{
  for (int i=0; i < vertices.rows(); ++i) {

    //place vertices
    vertices(i,0) = ((-(vertices.rows()-1)*(edgeLength))/2.0)+(i*edgeLength);
  }

  return vertices;
}

// obtain edges given vertex positions
MatrixXd updateEdges(MatrixXd v)
{
  // using Eigen slicing, subtract adjacent rows to obtain edges
  return v(seqN(1,v.rows()-1),all)-v(seqN(0,v.rows()-1),all);
}

// Function to obtain the angle between adjacent edges for discrete parallel transport when
// rotating around a curvature binormal
VectorXd getTheta(MatrixXd e)
{
  // set the endpoint vertices' cruvature to zero (since there is but one connected edge on each)
  VectorXd k = VectorXd::Constant(e.rows()+1,0);

  // return angle between adjacent edges
  k(seqN(1,e.rows()-1)) = ((e(seqN(0,e.rows()-1),all).rowwise().normalized()*
  			     e(seqN(1,e.rows()-1),all).rowwise().normalized().transpose()).array().acos()).matrix().diagonal();

  return k;
}


// get curvatures per each vertex
VectorXd getCurvatures(MatrixXd e)
{
  // set the endpoint vertices' cruvature to zero (since there is but one connected edge on each)
  VectorXd k = VectorXd::Constant(e.rows()+1,0);

  // for the vertices between the endpoints,take the normalized dotproduct of adjacent edges
  // convert to array so that the coefficient-wise arcos can be applied
  // calculate the tangent before converting back to a matrix and take the diagonal
  k(seqN(1,e.rows()-1)) = (((e(seqN(0,e.rows()-1),all).rowwise().normalized()*
  		     e(seqN(1,e.rows()-1),all).rowwise().normalized().transpose()).array().acos()/2.0)
  		     .tan()*2.0).matrix().diagonal();

  return k;
}

// simple function to obtain the tangent of an edge
MatrixXd getTangent(MatrixXd e)
{
  return e.rowwise().normalized();
}


// Function to obtain a matrix of curvature binormals for each vertex of a rod
// Rows correlate to a curvature binormal vector for a given vertex
MatrixXd getKb(MatrixXd e){
  // define a new curvature binormal matrix
  MatrixXd kb = MatrixXd::Constant((e.rows()+1),3,0);
  // obbtain scaled cross product of adjacent edges
  MatrixXd cross = (2*crossProduct(e(seqN(0,e.rows()-1),all),e(seqN(1,e.rows()-1),all)));
  // obtain norms and dot products of respective edges
  MatrixXd scaledNorm = ((rowNorm(e(seqN(0,e.rows()-1),all))*rowNorm(e(seqN(1,e.rows()-1),all)).transpose()).diagonal() + rowDot(e(seqN(0,e.rows()-1),all),e(seqN(1,e.rows()-1),all))).array().inverse().matrix();
  // variable scaling of a matrix is tricky in Eigen, so a diagonal matrix is constructed to obtain the curvature binormal
  MatrixXd normDiagonal = DiagonalMatrix<double,Eigen::Dynamic>(scaledNorm.asDiagonal());
  // multiply the inverse scaled norm by the cross product of adjacent edges and return
  kb(seqN(1,e.rows()-1),all) = normDiagonal*cross;

  return kb;
}

// Obtain the bending energy of a rod, assuming that it is naturally straight
double getBendingEnergy(MatrixXd e, MatrixXd kb, double alpha)
{
  double bend = 0;
  for (int i = 1; i < e.rows(); i++) {
    bend += ((alpha*kb.row(i).dot(kb.row(i)))*(1.0/(e.row(i-1).norm()+e.row(i).norm())));
  }

  return bend;
}

// Given a user-set number of vertices, interpolate a spline to start to create our polyline
MatrixXd interpolateSpline(MatrixXd &cps, int N, int D)
{
  MatrixXd v(N,3);
  // interpolate over B-spline with user's desired number of vertices
  Spline3d spline = SplineFitting<Spline3d>::Interpolate(cps,D);
  float timeStep = 0;
  for (int i = 0; i < N; i++) {
    timeStep += 1.0/N;
    // populate vertex matrix with interpolated spline points
    v.row(i) = spline(timeStep).transpose();
  }

  return v;
}

// Given user-defined parameters of height, cycles, radius, and vertex count, interpolate a 
// helix to gather vertices for our polyline
MatrixXd interpolateHelix(int N, int T, double timeStep, double radius)
{
  // N: how many vertices to place for each rotation of the unit circle
  // T: how many unit circle to interpolate (time value or z-coord)
  // timeStep: coordinate multiplier for T
  // radius: radius of the helix
  
  MatrixXd v(N*T,3);
  Vector3d p(3);

  v.row(0) << 0,0,0;
  int counter = 0;

  // find vertices along a helix, N points around each circle, T unit circles,
  // and each circle of the helix goes down in z direction by timeStep
  for (int i = 0; i < T; i++) {
    for (int j = 0; j < N; j++) {
      // populate vertex matrix with interpolated helix points
      p << radius*cos((j/double(N))*(2.0*M_PI)), radius*sin((j/double(N))*(2.0*M_PI)),
      (i*timeStep) + ((j+1.0)/N)*timeStep;
      v.row(counter) = p;
      counter++;
    }
  }

  return v;
}

// FUNCTION FOR UPDATING U & V VECTORS OF OUR TWIST-FREE BISHOP FRAME
// Our notion of discrete parallel transport is calculated and propogated forward here
void computeBishopFrame(MatrixXd e, VectorXd thetas, VectorXd k, MatrixXd kb, MatrixXd &u, MatrixXd &v)
{
  // define single vector to later derive vector perpendicular to tangent
  Vector3d zaxis(3);
  zaxis << 0,0,-1;
  // tangent of edge 0
  Vector3d t(3);
  t = e.row(0).normalized();
  // cross product of tangent with random vector
  Vector3d u0(3);
  u0 = t.cross(zaxis).normalized();
  // update u vector
  u.row(0) = u0;
  // obtain v vector as the cross product of t and u
  v.row(0) = t.cross(u0).normalized();

  // iterate over rod
  for (int i = 1; i < e.rows(); i++) {

    t = e.row(i).normalized();
    // don't change bishop frame if segment has no curve
    if(k(i) <= 0){
      u0 = u.row(i-1);
      u.row(i) = u0;
    }
    else{
      // define rotation matrix around curvature binormal
      AngleAxis<double> kbRot(thetas(i),kb.row(i));
      // in order to use .cross, u needs to be a hardcoded vector hence the u0
      u0 = u.row(i-1).normalized();
      u0 = (kbRot*u0).normalized();
      u.row(i) = u0;
    }
    // define v as the cross product of u and t
    v.row(i) = t.cross(u0).normalized();
  }
}

// Function to calculate midpoints of an edge to place an anchor point for mesh coordinate calculation
MatrixXd getMidpoints(MatrixXd vertices){

  MatrixXd mid = MatrixXd::Constant(vertices.rows()-1,3,0);

  // create a matrix that has the midpoint vertices of each edge
  mid = (vertices(seqN(0,vertices.rows()-1),all) + vertices(seqN(1,vertices.rows()-1),all))*(0.5);
    //for (int i = 0; i < vertices.rows()-1; ++i) {
    //mid(i, all) = (vertices(i, all) + vertices(i+1, all))/2;
    //}

  return mid;
}

MatrixXd getMeshVertices(MatrixXd A, MatrixXd B){

  // concatenate A and B where A are the vertices where we add a scaled version of u
  // and B are the vertices where we subtract a scale version of u
  Eigen::MatrixXd V= (Eigen::MatrixXd(A.rows()*2,3)<<
    A,B).finished();

  return V;
}

MatrixXi getMeshFaces(MatrixXd A){

  MatrixXi F = MatrixXi::Constant(A.rows()-2,3,0);

  // Generate faces matrix where we assign vertices counter-clockwise to constructed
  // the triangles that make up the mesh
  for (int i = 0, m = A.rows()/2, counter = 0; i < A.rows()/2-1; ++i, ++m, counter+=2) {
    F(counter, all) << i, m, m+1;
    F(counter+1, all) << i, m+1, i+1;
  }

  return F;
}

void updateFrame(igl::opengl::glfw::Viewer &viewer, MatrixXd &controlPoints,
  int &vertNum, int &splineDegree, double &thickness, double &bendingModulus,
  MatrixXd &kb, int &drawSpline, int &vertPerCircle, int &ringNum, double &ts, double &radius){

  MatrixXd vertices;

  if (drawSpline == 0){
    // obtain vertices
    vertices = interpolateSpline(controlPoints,vertNum, splineDegree);
  }
  else{
    vertices = interpolateHelix(vertPerCircle, ringNum, ts, radius);
  }

  MatrixXd edges = updateEdges(vertices);

  // curvature at each vertex -> magnitude of curvature binormal at said vertex
  VectorXd curvatures = getCurvatures(edges);
  VectorXd thetas = getTheta(edges);
  kb = getKb(edges); // curvature binormals

  MatrixXd u = MatrixXd::Constant(edges.rows(),3,0);
  MatrixXd v = MatrixXd::Constant(edges.rows(),3,0);

  computeBishopFrame(edges, thetas, curvatures, kb, u, v);

  MatrixXd midpoints = getMidpoints(vertices);

  MatrixXd vertpu = midpoints + u * thickness;
  MatrixXd vertmu = midpoints - u * thickness;

  MatrixXd VA = getMeshVertices(vertpu, vertmu);
  MatrixXi FA = getMeshFaces(VA);

  VectorXi K(static_cast<int>(FA.rows())/2);
  for (int i = 0; i < FA.rows(); i++) {
    if (i % 2 == 0) {
      K(i/2) = i;
    }
  }

  // default green for all faces
  MatrixXd C = RowVector3d(0.4,0.8,0.3).replicate(FA.rows(),1);
  // Red for each in K
  MatrixXd R = RowVector3d(1.0,0.3,0.3).replicate(K.rows(),1);
  // C(K,:) = R
  igl::slice_into(R,K,1,C);

  viewer.data().clear();
  viewer.data().set_mesh(VA, FA);
  viewer.data().set_colors(C);
}


int main(int argc, char *argv[])
{

  // beginning spline
  // define four default control points for B-Spline
  Vector3d p1(-121,29,59);
  Vector3d p2(-35,-54,-112);
  Vector3d p3(70,83,136);
  Vector3d p4(123,-37,134);
  // add control points to matrix
  static MatrixXd controlPoints(3,4);
  controlPoints << p1, p2, p3, p4;

  // Initialize number of vertices, thickness of mesh, bending modulus, and bending energy
  static int vertNum = 4;
  static int splineDegree = 2;
  static double bendingModulus = 0.5;
  static double thickness = 5;
  static double bendingEnergy = 0;
  static int drawSpline = 0;
  static int vertPerCircle = 20;
  static int ringNum = 5;
  static double ts = 10.0;
  static double radius = 5.0;

  // obtain vertices
  MatrixXd vertices = interpolateSpline(controlPoints,vertNum, splineDegree);\

  MatrixXd edges = updateEdges(vertices);

  // curvature at each vertex -> magnitude of curvature binormal at said vertex
  VectorXd curvatures = getCurvatures(edges);
  VectorXd thetas = getTheta(edges);
  MatrixXd kb = getKb(edges); // curvature binormals

  // Initialize u and v vectors which define normals perpendicular to the tangent
  MatrixXd u = MatrixXd::Constant(edges.rows(),3,0);
  MatrixXd v = MatrixXd::Constant(edges.rows(),3,0);

  // Compute the bishop frame
  computeBishopFrame(edges, thetas, curvatures, kb, u, v);

  // Given the spline interpolated vertices, find midpoints of each edge
  MatrixXd midpoints = getMidpoints(vertices);

  // Find the vertices of the mesh given the bishop frames along the centerline
  // and thickness
  MatrixXd vertpu = midpoints + u * thickness;
  MatrixXd vertmu = midpoints - u * thickness;

  // Define mesh vertices and mesh faces using helper functions
  MatrixXd VA = getMeshVertices(vertpu, vertmu);
  MatrixXi FA = getMeshFaces(VA);

  // Take the even rows of the full FA (Faces) matrix
  VectorXi K(static_cast<int>(FA.rows())/2);
  for (int i = 0; i < FA.rows(); i++) {
    if (i % 2 == 0) {
      K(i/2) = i;
    }
  }

  // default green for all faces
  MatrixXd C = RowVector3d(0.4,0.8,0.3).replicate(FA.rows(),1);
  // Red for each in K
  MatrixXd R = RowVector3d(1.0,0.3,0.3).replicate(K.rows(),1);

  // Change color matrix C so that even rows for FA (Faces) are now colored red
  igl::slice_into(R,K,1,C);

  // Initialize the viewer
  igl::opengl::glfw::Viewer viewer;

  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  // Customize the menu
  double doubleVariable = 0.1f;

  // Add content to the default menu window
  menu.callback_draw_viewer_menu = [&]()
  {
    // Draw parent menu content
    menu.draw_viewer_menu();

    // Add new group
    if (ImGui::CollapsingHeader("Model Parameters", ImGuiTreeNodeFlags_DefaultOpen))
    {
      // Update frame with either spline or helix mesh
      ImGui::RadioButton("Spline", &drawSpline, 0); ImGui::SameLine();
      ImGui::RadioButton("Helix", &drawSpline, 1);
      updateFrame(viewer, controlPoints, vertNum, splineDegree, thickness, bendingModulus,
        kb, drawSpline, vertPerCircle, ringNum, ts, radius);

      if (ImGui::TreeNode("Spline Parameters"))
      {
        // Update frame with new number of vertices through panel
        if (ImGui::InputInt("Num vertices", &vertNum))
        {
          vertNum = std::max(4, std::min(300, vertNum));
          updateFrame(viewer, controlPoints, vertNum, splineDegree, thickness, bendingModulus,
            kb, drawSpline, vertPerCircle, ringNum, ts, radius);
        }

        // Update frame with new degrees for spline
        if (ImGui::InputInt("Degree of the Spline", &splineDegree))
        {
          splineDegree = std::max(1, std::min(3, splineDegree));
          updateFrame(viewer, controlPoints, vertNum, splineDegree, thickness, bendingModulus,
            kb, drawSpline, vertPerCircle, ringNum, ts, radius);
        }

        ImGui::TreePop();
      }

      if (ImGui::TreeNode("Helix Parameters"))
      {

        // Update frame with new number of vertices per ring through panel
        if (ImGui::InputInt("Vertices per ring", &vertPerCircle))
        {
          vertPerCircle = std::max(2, std::min(50, vertPerCircle));
          updateFrame(viewer, controlPoints, vertNum, splineDegree, thickness, bendingModulus,
            kb, drawSpline, vertPerCircle, ringNum, ts, radius);
        }

        // Update frame with new z-step per ring through panel
        if (ImGui::InputInt("Number of rings", &ringNum))
        {
          ringNum = std::max(1, std::min(30, ringNum));
          updateFrame(viewer, controlPoints, vertNum, splineDegree, thickness, bendingModulus,
            kb, drawSpline, vertPerCircle, ringNum, ts, radius);
        }

        // Update frame with new z-step per ring through panel
        if (ImGui::InputDouble("Z-step per ring", &ts))
        {
          ts = std::max(5.0, std::min(50.0, ts));
          updateFrame(viewer, controlPoints, vertNum, splineDegree, thickness, bendingModulus,
            kb, drawSpline, vertPerCircle, ringNum, ts, radius);
        }

        // Update frame with new radius through panel
        if (ImGui::InputDouble("Helix Radius", &radius))
        {
          radius = std::max(1.0, std::min(30.0, radius));
          updateFrame(viewer, controlPoints, vertNum, splineDegree, thickness, bendingModulus,
            kb, drawSpline, vertPerCircle, ringNum, ts, radius);
        }

        ImGui::TreePop();
      }

      // Update frame with new thickness through panel
      if (ImGui::InputDouble("Thickness", &thickness))
      {
        thickness = std::max(0.1, std::min(100.0, thickness));
        updateFrame(viewer, controlPoints, vertNum, splineDegree, thickness, bendingModulus,
          kb, drawSpline, vertPerCircle, ringNum, ts, radius);
      }

      // Update frame with new alpha through panel
      if (ImGui::InputDouble("Bending modulus", &bendingModulus))
      {
        bendingModulus = std::max(0.0, bendingModulus);
        updateFrame(viewer, controlPoints, vertNum, splineDegree, thickness, bendingModulus,
          kb, drawSpline, vertPerCircle, ringNum, ts, radius);
      }

      // Calculate and display bending energy through panel
      ImGui::Text("Bending energy: %f",getBendingEnergy(edges, kb, bendingModulus));
    }

    // Add new group
    if (ImGui::CollapsingHeader("Control Points", ImGuiTreeNodeFlags_DefaultOpen))
    {
      // create sliders for each coordinate of each control point
      static ImGuiSliderFlags flags = ImGuiSliderFlags_None;
      ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.8f);

      static float p1x = static_cast<float>(controlPoints(0,0));
      ImGui::SliderFloat("P1 X", &p1x, -150.0, 150.0, "%.3f", flags);
      controlPoints(0,0) = static_cast<double>(p1x);

      static float p1y = static_cast<float>(controlPoints(1,0));
      ImGui::SliderFloat("P1 Y", &p1y, -150.0, 150.0, "%.3f", flags);
      controlPoints(1,0) = static_cast<double>(p1y);

      static float p1z = static_cast<float>(controlPoints(2,0));
      ImGui::SliderFloat("P1 Z", &p1z, -150.0, 150.0, "%.3f", flags);
      controlPoints(2,0) = static_cast<double>(p1z);

      static float p2x = static_cast<float>(controlPoints(0,1));
      ImGui::SliderFloat("P2 X", &p2x, -150.0, 150.0, "%.3f", flags);
      controlPoints(0,1) = static_cast<double>(p2x);

      static float p2y = static_cast<float>(controlPoints(1,1));
      ImGui::SliderFloat("P2 Y", &p2y, -150.0, 150.0, "%.3f", flags);
      controlPoints(1,1) = static_cast<double>(p2y);

      static float p2z = static_cast<float>(controlPoints(2,1));
      ImGui::SliderFloat("P2 Z", &p2z, -150.0, 150.0, "%.3f", flags);
      controlPoints(2,1) = static_cast<double>(p2z);

      static float p3x = static_cast<float>(controlPoints(0,2));
      ImGui::SliderFloat("P3 X", &p3x, -150.0, 150.0, "%.3f", flags);
      controlPoints(0,2) = static_cast<double>(p3x);

      static float p3y = static_cast<float>(controlPoints(1,2));
      ImGui::SliderFloat("P3 Y", &p3y, -150.0, 150.0, "%.3f", flags);
      controlPoints(1,2) = static_cast<double>(p3y);

      static float p3z = static_cast<float>(controlPoints(2,2));
      ImGui::SliderFloat("P3 Z", &p3z, -150.0, 150.0, "%.3f", flags);
      controlPoints(2,2) = static_cast<double>(p3z);

      static float p4x = static_cast<float>(controlPoints(0,3));
      ImGui::SliderFloat("P4 X", &p4x, -150.0, 150.0, "%.3f", flags);
      controlPoints(0,3) = static_cast<double>(p4x);

      static float p4y = static_cast<float>(controlPoints(1,3));
      ImGui::SliderFloat("P4 Y", &p4y, -150.0, 150.0, "%.3f", flags);
      controlPoints(1,3) = static_cast<double>(p4y);

      static float p4z = static_cast<float>(controlPoints(2,3));
      ImGui::SliderFloat("P4 Z", &p4z, -150.0, 150.0, "%.3f", flags);
      controlPoints(2,3) = static_cast<double>(p4z);

      updateFrame(viewer, controlPoints, vertNum, splineDegree, thickness, bendingModulus,
        kb, drawSpline, vertPerCircle, ringNum, ts, radius);
    }
  };

  // Plot the mesh
  viewer.data().set_mesh(VA, FA);
  viewer.data().set_colors(C);
  viewer.launch();

  return 0;
}
