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
#include "axisym_spherical_solid.h"

//The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Physical_Variables
{

 /// thickness
 double t=0.1;

 /// degrees
 double theta=1.57079632679;

 //number of elements
 unsigned n_ele_r = 4;
 unsigned n_ele_theta = 40;
 /// Pointer to strain energy function
 StrainEnergyFunction*Strain_energy_function_pt;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt;

 /// Elastic modulus
 double E=1.0;

 /// Poisson's ratio
 double Nu=0.3;

 //Constants for Mooney-Rivlin Material constants
 double C1 = 1.3;
 double C2 = 1.1;

 /// Uniform pressure
 double P = 0.0;
 
 //Compression Height H
 double H = 1.0;
 
 bool printDebugInfo = false;
 bool printTractionInfo = true;

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
    if(printTractionInfo)
    {
        cout << "i = " << i <<"  P = " << P << "       ";
        cout << "n[" << i << "] = " << n[i] << endl;
    }
   }
   printTractionInfo = false;
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
 CantileverProblem(const bool& incompress);
 
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

 ///Run it with increasing pressure
 void run_function(const unsigned nstep, const double p_increment);

private:

 /// \short Pass pointer to traction function to the
 /// elements in the traction mesh
 void set_traction_pt();


 /// \short Create contact elements
 /// This will impose the boundary condition on the top
 void create_contact_elements();

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
CantileverProblem<ELEMENT>::CantileverProblem(const bool& incompress) 
{

 // Create the mesh
 //Max_residuals = 100.0;

 // # of elements in x-direction
 unsigned n_r=Global_Physical_Variables::n_ele_r;

 // # of elements in y-direction
 unsigned n_theta=Global_Physical_Variables::n_ele_theta;

 // Domain length in x-direction
 double l_r= Global_Physical_Variables::t; //0.2

 // Domain length in y-direction
 double l_theta=Global_Physical_Variables::theta; //pi/2

 // Shift origin in r direction
 Vector<double> origin(2);
 origin[0]=1.0-l_r; 
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
           
           //I need to loop over deformed coordinates and set them correctly
           // i.e. R=r but Theta = 0
           unsigned nnod_el=solid_mesh_pt()->finite_element_pt(i)->nnode();
           // Loop over all nodes in element
           for (unsigned j=0;j<nnod_el;j++)
            {
            	// set ri = r thetai=0
            	solid_mesh_pt()->finite_element_pt(i)->node_pt(j)->x(1) = 0;
            }
           
           //Set the body force
           //el_pt->body_force_fct_pt() = Global_Physical_Variables::gravity;
          

            // Is it incompressible 
           if (incompress)
            {     
                     AxisymmetricPVDEquationsWithPressure* test_pt = 
                      dynamic_cast<AxisymmetricPVDEquationsWithPressure*>(
                       solid_mesh_pt()->element_pt(i));
                     if (test_pt!=0)
                      {
                       test_pt->set_incompressible();
                      }
            }
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
 
 //Creating traction mesh, commented to see if newton solver converges without
   //Traction face elements
 // Construct the traction element mesh
 Traction_mesh_pt=new SolidMesh;
 if(Global_Physical_Variables::printDebugInfo){
 	cout << "Going into create_traction_elements()" << endl;
  }
 create_traction_elements();
 
 if(Global_Physical_Variables::printDebugInfo){
 	cout << "After create_traction_elements()" << endl;
  }
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


 unsigned boundary_num = 2;
 // Pin the  boundary on symmetry in theta directions
 unsigned n_side = mesh_pt()->nboundary_node(boundary_num);

 //Loop over the nodes
 for(unsigned i=0;i<n_side;i++)
  {
   //solid_mesh_pt()->boundary_node_pt(boundary_num,i)->pin_position(0);
   solid_mesh_pt()->boundary_node_pt(boundary_num,i)->pin_position(1);
  }

 boundary_num = 0;
 n_side = mesh_pt()->nboundary_node(boundary_num);
 
 //Loop over the nodes
 for(unsigned i=0;i<n_side;i++)
  {
   //solid_mesh_pt()->boundary_node_pt(boundary_num,i)->pin_position(0);
   solid_mesh_pt()->boundary_node_pt(boundary_num,i)->pin_position(1);
  } //

   if(Global_Physical_Variables::printDebugInfo){
 	cout << "BC set " << endl;
   }
  


 /* solid pressure not implemented in a general way
 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  solid_mesh_pt()->element_pt());
 */
 
 //Attach the boundary conditions to the mesh
 cout << "assign_eqn_numbers()  = " << assign_eqn_numbers() << std::endl; 

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
  if(Global_Physical_Variables::printDebugInfo){
 	cout << "Starting  set_traction_pt()" << endl;
  }
 // Loop over the elements in the traction element mesh
 // for elements on the top boundary 
 unsigned n_element=traction_mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid traction element
   AxisymmetricSolidTractionElement<ELEMENT> *el_pt = 
    dynamic_cast<AxisymmetricSolidTractionElement<ELEMENT>*>
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
 unsigned b=3;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = solid_mesh_pt()->nboundary_element(b);
 if(Global_Physical_Variables::printDebugInfo){
 	cout << "Going to create " << n_element << " elements." << endl;
  }
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
 //if(Global_Physical_Variables::printDebugInfo){
 //	cout << "Starting " << e << " ..." ;
  //}
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    solid_mesh_pt()->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b
   int face_index = solid_mesh_pt()->face_index_at_boundary(b,e);
      
   // Create new element and add to mesh
   //Traction_mesh_pt->add_element_pt(new SolidTractionElement<ELEMENT>
   //                                 (bulk_elem_pt,face_index));

   //Use axisymmetric traction elements
   Traction_mesh_pt->add_element_pt(new AxisymmetricSolidTractionElement<ELEMENT>
                                    (bulk_elem_pt,face_index));
  //if(Global_Physical_Variables::printDebugInfo){
 //	cout << " ... finished! " << endl ;
  //}
 
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
 sprintf(filename,"%s/t_%.2f_MR_C1%.2f_C2%.2f_P%.5f_soln%i.dat",Doc_info.directory().c_str(),
        Global_Physical_Variables::t,        
        Global_Physical_Variables::C1, Global_Physical_Variables::C2 , 
        Global_Physical_Variables::P, Doc_info.number());
 //sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
 //       Doc_info.number());
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

template<class ELEMENT>
void CantileverProblem<ELEMENT>::run_function(const unsigned nstep, const double p_increment){

 Global_Physical_Variables::P =0.0;
 //solve once in undeformed configuration, hopefully this should converge immediatly
 newton_solve();

  //Lets check that displacement is zero
 doc_solution();
 
 for(unsigned i=0;i<nstep;i++)
  {
   // Increment pressure load
   Global_Physical_Variables::P+=p_increment;
	
   cout << endl;
   cout << "Step " << i << " out of " << nstep << 
        ". Incrementing pressure to " << Global_Physical_Variables::P << endl;
   // Solve the problem    
   newton_solve();

   // Doc solution
   doc_solution();

  }

}

//=======start_of_main==================================================
/// Driver for cantilever beam loaded by surface traction and/or
/// gravity
//======================================================================
int main()
{

 bool incompressible = true;
 if(Global_Physical_Variables::printDebugInfo){
 	cout << "Starting Main " << endl;
 }

 //Global_Physical_Variables::n_ele_theta = 60;

 // Create generalised Hookean constitutive equations
 //Global_Physical_Variables::Constitutive_law_pt = 
 // new GeneralisedHookean(&Global_Physical_Variables::Nu,
 //                        &Global_Physical_Variables::E);

 /*
 Global_Physical_Variables::Strain_energy_function_pt = 
 new GeneralisedMooneyRivlin(&Global_Physical_Variables::C1,
              &Global_Physical_Variables::C2);
 */

 // Create MooneyRivlin constitutive equations
 Global_Physical_Variables::Strain_energy_function_pt = 
 new MooneyRivlin(&Global_Physical_Variables::C1,
                        &Global_Physical_Variables::C2);

  // Define a constitutive law (based on strain energy function)
 Global_Physical_Variables::Constitutive_law_pt = 
 new IsotropicStrainEnergyFunctionConstitutiveLaw(
  Global_Physical_Variables::Strain_energy_function_pt);


 

 //Set up the problem
 //CantileverProblem<MySolidElement<RefineableQPVDElement<2,3> > > problem;
 
// CantileverProblem<AxisymQPVDElement> problem;
// CantileverProblem<AxisymQPVDElementWithPressure> problem(incompressible);
 //CantileverProblem<AxisymDiagHermitePVDElement> problem;

 /*
 if(Global_Physical_Variables::printDebugInfo){
 	cout << "Build Problem  with incompressible = " << incompressible << endl;
 }
 // Output boundaries 
 problem.solid_mesh_pt()->output_boundaries("RESLT/boundaries.dat");
*/
 // Initial values for parameter values
 Global_Physical_Variables::P=0.0; 
 Global_Physical_Variables::Gravity=0.0;

 // Max. number of adaptations per solve
 //unsigned max_adapt=0;
 
  if(Global_Physical_Variables::printDebugInfo){
        cout << " " << endl;
 	cout << "Initial Newton solve: " << endl;
 }


 Global_Physical_Variables::t = 0.02;
 Global_Physical_Variables::n_ele_r = 2;
 Global_Physical_Variables::n_ele_theta =(int) (1.571/(Global_Physical_Variables::t/    Global_Physical_Variables::n_ele_r)  + 0.5);
 cout << "Thickness = " << Global_Physical_Variables::t << " with " <<
        Global_Physical_Variables::n_ele_r << ", " << 
        Global_Physical_Variables::n_ele_theta << 
        " elements in the r and theta direction" << endl;
 CantileverProblem<AxisymQPVDElementWithPressure> problem2(incompressible);


 Global_Physical_Variables::P =0.0;
 //solve once in undeformed configuration, hopefully this should converge immediatly
 problem2.newton_solve();

 //Lets check that displacement is zero
 problem2.doc_solution();
 
 unsigned nstep=2;
 for(unsigned jj=0;jj<nstep;jj++)
  {
   // Increment pressure load
   Global_Physical_Variables::P+=0.01;
	
   cout << endl;
   cout << "Step " << jj << " out of " << nstep << 
        ". Incrementing pressure to " << Global_Physical_Variables::P << endl;
   // Solve the problem    
   problem2.newton_solve();

   // Doc solution
   problem2.doc_solution();

  }


// problem2.run_function(400, 2.5e-3);


return 0;
 
} //end of main








