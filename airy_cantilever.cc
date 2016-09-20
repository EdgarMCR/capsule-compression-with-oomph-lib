//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
//Driver for Airy cantilever beam problem

//Oomph-lib includes
#include "generic.h"
#include "solid.h"
#include "constitutive.h"

//The mesh
#include "meshes/rectangular_quadmesh.h"
//#include "axisym_spherical_solid.h"
#include "axisym_solid_elements.h"
//#include "axisym_solid_traction_elements.h"

using namespace std;

using namespace oomph;

//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Physical_Variables
{

 /// thickness
 double t=0.5;

 /// degrees
 double theta=0.78539816339;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt;

 /// Elastic modulus
 double E=1.0;

 /// Poisson's ratio
 double Nu=0.3;

 /// Uniform pressure
 double P = 0.0;
 
 //Compression Height H
 double H = 1.0;
 
 bool printDebugInfo = false;

 /// \short Constant pressure load. The arguments to this function are imposed
 /// on us by the SolidTractionElements which allow the traction to 
 /// depend on the Lagrangian and Eulerian coordinates x and xi, and on the 
 /// outer unit normal to the surface. Here we only need the outer unit
 /// normal.
 void constant_pressure(const Vector<double> &xi, const Vector<double> &x,
                        const Vector<double> &n, Vector<double> &traction)
 {
  unsigned dim = traction.size();
  for(unsigned i=0;i<dim;i++)
   {
    traction[i] = -P*n[i];
   }
 } // end traction


 /// Non-dim gravity
 double Gravity=0.0;

 /// Non-dimensional gravity as body force
 void gravity(const double& time, 
              const Vector<double> &xi, 
              Vector<double> &b)
 {
  b[0]=0.0;
  b[1]=-Gravity;
 }
 
} //end namespace



//=============begin_problem============================================ 
/// Problem class for the cantilever "beam" structure.
//====================================================================== 
template<class ELEMENT>
class CantileverProblem : public Problem
{

public:

 /// Constructor:
 CantileverProblem();
 
 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// Update function (empty)
 void actions_before_newton_solve() {}
 
  /// Update function (empty)
 void actions_before_adapt(){}

 /// Update function (empty)
 void actions_after_adapt(){}

 /// Access function for the solid mesh
 ElasticRectangularQuadMesh<ELEMENT>*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 
  /*
 /// Access function for the solid mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 
  */
  
 /// Access function to the mesh of surface traction elements
 SolidMesh*& traction_mesh_pt(){return Traction_mesh_pt;} 

 /// Doc the solution
 void doc_solution();

private:

 /// \short Pass pointer to traction function to the
 /// elements in the traction mesh
 void set_traction_pt();

 /// \short Create traction elements
 void create_traction_elements();

 /// Delete traction elements
 void delete_traction_elements();

 /// Trace file
 ofstream Trace_file;
 
 /// Pointers to node whose position we're tracing
 Node* Trace_node_pt;

 /*
 /// Pointer to solid mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>* Solid_mesh_pt;
 */
 
 /// Pointer to solid mesh
 ElasticRectangularQuadMesh<ELEMENT>* Solid_mesh_pt;

 /// Pointer to mesh of traction elements
 SolidMesh* Traction_mesh_pt;

 /// DocInfo object for output
 DocInfo Doc_info;
 
};


//===========start_of_constructor======================================= 
/// Constructor: 
//====================================================================== 
template<class ELEMENT>
CantileverProblem<ELEMENT>::CantileverProblem() 
{

 // Create the mesh
 //Max_residuals = 400.0;

 // # of elements in x-direction
 unsigned n_r=4;

 // # of elements in y-direction
 unsigned n_theta=20;

 // Domain length in x-direction
 double l_r= Global_Physical_Variables::t;

 // Domain length in y-direction
 double l_theta=Global_Physical_Variables::theta;

 // Shift mesh downwards so that centreline is at y=0:
 Vector<double> origin(2);
 origin[0]=0.0;
 origin[1]=0.0 ; //-0.5*l_y;
 
 /*
 //Now create the mesh 
 solid_mesh_pt() = new ElasticRefineableRectangularQuadMesh<ELEMENT>(
  n_r,n_theta,l_r,l_theta,origin);

 // Set error estimator
 solid_mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 */
 
  if(Global_Physical_Variables::printDebugInfo){
 	cout << "Now create the mesh " << endl;
  }
  //Now create the mesh 
 solid_mesh_pt() = new ElasticRectangularQuadMesh<ELEMENT>(
  n_r,n_theta,l_r,l_theta,origin);
  
     // Now set the Eulerian equal to the Lagrangian coordinates
   solid_mesh_pt()->set_lagrangian_nodal_coordinates();
  
   if(Global_Physical_Variables::printDebugInfo){
 	cout << "Created the mesh " << endl;
  }

 //Assign the physical properties to the elements before any refinement
 //Loop over the elements in the main mesh
 unsigned n_element =solid_mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(i));
   
   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Constitutive_law_pt;

   //Set the body force
   //el_pt->body_force_fct_pt() = Global_Physical_Variables::gravity;
  }

  if(Global_Physical_Variables::printDebugInfo){
 	cout << "Set constitutive law pointer " << endl;
  }

 // Choose a control node: The last node in the solid mesh
 unsigned nnod=solid_mesh_pt()->nnode();
 Trace_node_pt=solid_mesh_pt()->node_pt(nnod-1);
 
 /*
 // Refine the mesh uniformly
 solid_mesh_pt()->refine_uniformly();
 */
 
 // Construct the traction element mesh
 Traction_mesh_pt=new SolidMesh;
 create_traction_elements();
 
 // Pass pointer to traction function to the elements
 // in the traction mesh
 set_traction_pt();
 
 // Solid mesh is first sub-mesh
 add_sub_mesh(solid_mesh_pt());

 // Add traction sub-mesh
 add_sub_mesh(traction_mesh_pt());

 // Build combined "global" mesh
 build_global_mesh();
 
 if(Global_Physical_Variables::printDebugInfo){
 	cout << "build global mesh " << endl;
  }
 
 // Pin the  boundary on symmetry in theta directions
 unsigned n_side = mesh_pt()->nboundary_node(3);
 
 //Loop over the nodes
 for(unsigned i=0;i<n_side;i++)
  {
   //solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(0);
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(1);
  }
  
 n_side = mesh_pt()->nboundary_node(1);
 
 //Loop over the nodes
 for(unsigned i=0;i<n_side;i++)
  {
   //solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(0);
   solid_mesh_pt()->boundary_node_pt(1,i)->pin_position(1);
  }

   if(Global_Physical_Variables::printDebugInfo){
 	cout << "BC set " << endl;
  }
  
 /* 
 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  solid_mesh_pt()->element_pt());
 */
 
 //Attach the boundary conditions to the mesh
 cout << assign_eqn_numbers() << std::endl; 

 // Set output directory
 Doc_info.set_directory("RESLT");

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",Doc_info.directory().c_str());
 Trace_file.open(filename);
 
  if(Global_Physical_Variables::printDebugInfo){
 	cout << "constructor finished " << endl;
  }
 

} //end of constructor




//==================start_of_set_traction_pt==============================
/// Set pointer to traction function for the relevant
/// elements in the traction mesh
//========================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::set_traction_pt()
{
 // Loop over the elements in the traction element mesh
 // for elements on the top boundary (boundary 2)
 unsigned n_element=traction_mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid traction element
   SolidTractionElement<ELEMENT> *el_pt = 
    dynamic_cast<SolidTractionElement<ELEMENT>*>
    (traction_mesh_pt()->element_pt(i));

   //Set the traction function
   el_pt->traction_fct_pt() = Global_Physical_Variables::constant_pressure;
  }
 
}// end of set traction pt


 
//============start_of_create_traction_elements==============================
/// Create traction elements 
//=======================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::create_traction_elements()
{
 // Traction elements are located on boundary 2:
 unsigned b=2;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = solid_mesh_pt()->nboundary_element(b);
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    solid_mesh_pt()->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b
   int face_index = solid_mesh_pt()->face_index_at_boundary(b,e);
      
   // Create new element and add to mesh
   Traction_mesh_pt->add_element_pt(new SolidTractionElement<ELEMENT>
                                    (bulk_elem_pt,face_index));   
  }  

 // Pass the pointer to the traction function to the traction elements
 set_traction_pt();
 
} // end of create_traction_elements




//============start_of_delete_traction_elements==============================
/// Delete traction elements and wipe the  traction meshes
//=======================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::delete_traction_elements()
{
 // How many surface elements are in the surface mesh
 unsigned n_element = Traction_mesh_pt->nelement();
 
 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete Traction_mesh_pt->element_pt(e);
  }
 
 // Wipe the mesh
 Traction_mesh_pt->flush_element_and_node_storage();

} // end of delete_traction_elements



//==============start_doc===========================================
/// Doc the solution
//==================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::doc_solution()
{

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned n_plot = 5; 

 // Output shape of and stress in deformed body
 //--------------------------------------------
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 solid_mesh_pt()->output(some_file,n_plot);
 some_file.close();

 // Write trace file: Load/displacement characteristics
 Trace_file << Global_Physical_Variables::P  << " " 
            << Trace_node_pt->x(0) << " " 
            << Trace_node_pt->x(1) << " " 
            << std::endl;

 // Increment label for output files
 Doc_info.number()++;

} //end doc



//=======start_of_main==================================================
/// Driver for cantilever beam loaded by surface traction and/or
/// gravity
//======================================================================
int main()
{
 if(Global_Physical_Variables::printDebugInfo){
 	cout << "Starting Main " << endl;
 }
 // Create generalised Hookean constitutive equations
 Global_Physical_Variables::Constitutive_law_pt = 
  new GeneralisedHookean(&Global_Physical_Variables::Nu,
                         &Global_Physical_Variables::E);
 
 //Set up the problem
 //CantileverProblem<MySolidElement<RefineableQPVDElement<2,3> > > problem;
 
 CantileverProblem<AxisymQPVDElement> problem;
 
 if(Global_Physical_Variables::printDebugInfo){
 	cout << "Build Problem " << endl;
 }

 // Uncomment these as an exercise

 // CantileverProblem<MySolidElement<
 //  RefineableQPVDElementWithContinuousPressure<2> > > problem;

 // CantileverProblem<MySolidElement<
 //  RefineableQPVDElementWithPressure<2> > > problem;

 // Initial values for parameter values
 Global_Physical_Variables::P=0.0; 
 Global_Physical_Variables::Gravity=0.0;

 // Max. number of adaptations per solve
 //unsigned max_adapt=0;
 
  if(Global_Physical_Variables::printDebugInfo){
 	cout << "Initial Newton solve: " << endl;
 }
 //solve once in undeformed configuration, hopefully this should converge immediatly
 problem.newton_solve();
 
 //Parameter incrementation
 unsigned nstep=5; 
 double p_increment=1.0e-5;   
 for(unsigned i=0;i<nstep;i++)
  {
   // Increment pressure load
   Global_Physical_Variables::P+=p_increment;

   // Solve the problem with Newton's method, allowing
   // up to max_adapt mesh adaptations after every solve.
   //problem.newton_solve(max_adapt);
   
   problem.newton_solve();

   // Doc solution
   problem.doc_solution();

  }
 
} //end of main








