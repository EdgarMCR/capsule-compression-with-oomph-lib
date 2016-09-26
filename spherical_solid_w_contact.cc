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


#include "axisym_solid_traction_elements.h"
#include "axisym_solid_elements.h"


//Oomph-lib includes
#include "generic.h"

//Solid and constitutive stuff
#include "solid.h"
#include "constitutive.h"

//axis symmetric elements
//#include "axisym_spherical_solid.h"

//The mesh
#include "meshes/rectangular_quadmesh.h"

// Contact stuff
#include "contact_elements.h"

using namespace std;

using namespace oomph;

/// To redirect a copy of cout to a log file
// If this gives an error, need to install boost
// 'sudo apt-get install libboost-iostreams-dev' 
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

using std::ostream;
using std::ofstream;
//using std::cout;

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;


//=============================================================================
/// One-dimensional integration scheme that locates integration points
/// uniformly along boundary (if elements are equal sized, that is...)
//=============================================================================
class MyIntegral : public Integral
{

public:

 /// \short Constructor: Specify number of uniformly spaced sample points in 
 /// unit interval
 MyIntegral(const unsigned& n_knot){N_knot=n_knot;}

 /// Broken copy constructor
 MyIntegral(const MyIntegral& dummy) 
  { 
   BrokenCopy::broken_copy("MyIntegral");
  } 
 
 /// Broken assignment operator
 void operator=(const MyIntegral&) 
  {
   BrokenCopy::broken_assign("MyIntegral");
  }

 /// Number of integration points of the scheme
 unsigned nweight() const {return N_knot;}

 /// \short Return coordinate s[j] (j=0) of integration point i -- 
 double knot(const unsigned &i, const unsigned &j) const 
  {
    
   // Doesnt s go from -1 to 1 ????
   //double dx=1.0/double(N_knot-1);
   //return ((double(i)))*dx;
   
   //assume -1 < s < 1
   double dx=(2.0-1e-8)/double(N_knot-1);
   return -1.0 + ((double(i)))*dx+1e-8;
   
  }
 
 /// Return weight of integration point i
 double weight(const unsigned &i) const 
  {
   //return 2.0/double(N_knot);
   return 1.0/double(N_knot);
  }
 
private:
 
 /// Number of integration points
 unsigned N_knot;
};


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
 
 ///Desired Volume
 double Volume = 0.243;
 
 //Compression Height H
 double H = 1.0;
 
  //Preinflation of the capsule, 1 being the undeformed state
  double lambda = 1.0; 

  // 0 =false, != 0  true 
 int printDebugInfo = 0;  
 bool printTractionInfo = false;

  /// Set Newton Tolerance via command line
  double newton_tol = 1e-8;

  /// Consitiutive Law string
  /// MR - Moony-Rivlin
  /// GH - Generalized Hookian
  std::string  constitutive_law = "MR"; 

  /// String to hold path to solution if the program is
  /// desired to restart from a saved state
  /// N.B. All parameters need to be set to same state manually
  std::string loadSol = ""; 


  /// string that allows different programs to be excuted
  /// Current options are:
  /// standard - lowers height adaptivly by 60%
  /// lowerC1 - lowers the value of C1 in the MoonyRivlin constitutive law
  std::string program = "standard";

 // Switches code between either using traction elements that impose a
 // set pressure,P, defined above or enforces a volume constraint set by
 // the variable Volume defined above.
 int enforce_volume_constraint = false;

  //Boolean values to set the options for the contact elements
  // Default is false = 0
  ///Set whether or not to use isoparametric basis function for pressure
 int Contact_use_isoparametric_flag;

  ///Set options for basis/test functions for penetration and pressure
  int Contact_use_collocated_penetration_flag;

  ///Set options for basis/test functions for penetration and pressure
  int Contact_use_collocated_contact_pressure_flag;

  // bool that is true, set the bool pointers to this when setting them to true
  // TODO: this should be constant but when I set it to constant that causes
  // problems with the setter functions in the contact_elements.h
  // Long story short, if you are reading this because something modified this
  // variable without you noticing, I'm sorry. 
  bool true_flag = true;

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
    traction[i] = -P*n[i]*-1;   //need the negative sign because normal calculation
                                // can be wrong by a minus sign 
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

 /// Penetrator
 Penetrator* Penetrator_pt=0;
 
 
 unsigned Nintpt = 100;
 bool addMoreGaussPoints = true;
 
 unsigned MaxNSteps= 10;
 double MaxResidual = 10000;
 
 
 // Max number of adaptations per newton solve
 unsigned max_adapt = 1;

  // double to hold unde-relaxation factor to be set via
  // command line
  double under_relaxation=0.4;

  //Pressure pre-factor for the contact equations
  double contact_pressure_prefactor = 1;
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
 void actions_before_newton_solve() {
   newton_step=0;

   //Update Volume in case lambda has been changed
   Global_Physical_Variables::Volume = calc_inflated_vol(Global_Physical_Variables::t, Global_Physical_Variables::lambda);
    /*
   std::cout << "solid mesh" << std::endl;
    /// Loop over all nodes in the mesh
    unsigned n_node = Solid_mesh_pt->nnode();
    for(unsigned inod=0;inod<n_node;inod++)
    {
     SolidNode* nod_pt=dynamic_cast<SolidNode*>(Solid_mesh_pt->node_pt(inod));
     
     std::cout << "(" <<nod_pt->x(0) << ", " << nod_pt->xi(1) + nod_pt->x(1) << "): ";
     
      std::cout <<" | ";
      unsigned n_value=nod_pt->nvalue();
      for(unsigned ival=0;ival<n_value;ival++){
        std::cout<<nod_pt->eqn_number(ival)<<" ";
      }
      std::cout <<" | ";
      
      unsigned n_position_type=nod_pt->nposition_type();
      for(unsigned k=0;k<n_position_type;k++)
      {
        for(unsigned i=0;i<2;i++)
        {
          std::cout << nod_pt->position_eqn_number(k,i)<<" ";
        }
      }

      std::cout<<std::endl;
    }
    
    
    if(Global_Physical_Variables::enforce_volume_constraint){
            /// Loop over all nodes in the mesh
        unsigned n_node = Vol_const_mesh_pt->nnode();
        
        cout << " Vol_const_mesh_pt->nnode()  = " << n_node << 
          " and Vol_const_mesh_pt->nelement() = " << 
          Vol_const_mesh_pt->nelement() << endl;
           

        AxisymmetricSolidTractionVolumeConstraintElement<ELEMENT>* el_pt = dynamic_cast<AxisymmetricSolidTractionVolumeConstraintElement<ELEMENT>* >(
            Vol_const_mesh_pt->element_pt(0));
        cout << "Vol_const_mesh_pt->element_pt(0)->node_pt(0)->x_position_pt(0) = " << el_pt->node_pt(0)->x(0) << endl;
        
        //Loop over all elements and then all nodes
        // This will mean nodes on the egde will eb duplicated but
        // for some reason the mesh function nnode() returns 0 
        // even when the mesh is filled with elements? 
        // TODO: track down issue with nnode() called from mesh
        unsigned nel = Vol_const_mesh_pt->nelement();
        for(unsigned e=0; e<nel;e++){
          AxisymmetricSolidTractionVolumeConstraintElement<ELEMENT>* el_pt = dynamic_cast<AxisymmetricSolidTractionVolumeConstraintElement<ELEMENT>* >(
            Vol_const_mesh_pt->element_pt(e));
          std::cout << "ele " << e << " # Ndof = " << el_pt->ndof() << " int: " 
            << el_pt->ninternal_data() << " ext: " <<  el_pt->nexternal_data() << std::endl; 
          unsigned nr_node = el_pt->nnode();
          for(unsigned inod=0;inod<nr_node;inod++)
              {
                
              SolidNode* nod_pt=dynamic_cast<SolidNode*>(el_pt->node_pt(inod));
              std::cout << nod_pt->x(0) << " " << nod_pt->xi(1) + nod_pt->x(1) << " ";
//              std::cout << nod_pt->eqn_number(0) << " ";
              
              std::cout <<" | ";
              unsigned n_value=nod_pt->nvalue();
              for(unsigned ival=0;ival<n_value;ival++){
                std::cout<<nod_pt->eqn_number(ival)<<" ";
              }
              std::cout <<" | ";
          
               unsigned n_position_type=nod_pt->nposition_type();
               for(unsigned k=0;k<n_position_type;k++)
                {
                 for(unsigned i=0;i<2;i++)
                  {
          //        if(nod_pt->position_eqn_number(k,i)>=0)
          //         {
                   std::cout << nod_pt->position_eqn_number(k,i)<<" ";
          //         }
                  }
                }

              std::cout<<std::endl;
              }
        }
    }
    
    AxisSymSolidVolumeConstraintElement<ELEMENT>* el_pt = 
      dynamic_cast<AxisSymSolidVolumeConstraintElement<ELEMENT>* >(
      Vol_const_master_mesh_pt->element_pt(0));
    
    std::cout << "Master Vol Control Ele p = " << el_pt->internal_data_pt(0)->value(0) <<
      " with eqn nr " << el_pt->internal_data_pt(0)->eqn_number(0) << ".";
        std::cout<<std::endl;
    
        std::cout << " el_pt->ndof() = " << el_pt->ndof() << 
          " el_pt->ninternal_data() = " <<
          el_pt->ninternal_data() << " el_pt->nexternal_data() " <<  
          el_pt->nexternal_data() << std::endl; 
          
//     DenseMatrix<double> ele_jacobian;
//     Vector<double> ele_residuals;          
//     el_pt->get_jacobian(ele_residuals,ele_jacobian);
// 
//     char filename_0[100];
//     sprintf(filename_0,"%s/ele_jacobian%i.dat",
//             Doc_info.directory().c_str(),
//             Doc_info.number()); 

          
    std::cout << "Contact mesh" << std::endl;
   n_node = Surface_contact_mesh_pt->nnode();
   std::cout << "n_node = " << n_node << std::endl;
   
    unsigned nel = Surface_contact_mesh_pt->nelement();
    for(unsigned e=0; e<nel;e++){
        AxiSymNonlinearSurfaceContactElement<ELEMENT>* el_pt = 
          dynamic_cast<AxiSymNonlinearSurfaceContactElement<ELEMENT>*>(
            Surface_contact_mesh_pt->element_pt(e));
          
        unsigned nr_node = el_pt->nnode();
        
        for(unsigned inod=0;inod<nr_node;inod++){
          SolidNode* nod_pt=dynamic_cast<SolidNode*>(el_pt->node_pt(inod));
          std::cout << "(" <<nod_pt->x(0) << ", " << nod_pt->xi(1) + nod_pt->x(1) << "): ";
          
          std::cout <<" | ";
          unsigned n_value=nod_pt->nvalue();
          for(unsigned ival=0;ival<n_value;ival++){
            std::cout<<nod_pt->eqn_number(ival)<<" ";
          }
          std::cout <<" | ";
          

            unsigned n_position_type=nod_pt->nposition_type();
            for(unsigned k=0;k<n_position_type;k++)
            {
              for(unsigned i=0;i<2;i++)
              {
                std::cout << nod_pt->position_eqn_number(k,i)<<" ";
              }
            }

            std::cout<<std::endl;
          }
    }
    std::cout << std::endl;
    */
    /*
   ofstream descr_file;
    char filename_dof[100];
    sprintf(filename_dof,"%s/dof%i.dat",
            Doc_info.directory().c_str(),
            Doc_info.number());
   descr_file.open(filename_dof);
   describe_dofs(descr_file);
   descr_file.close();

      
    DenseDoubleMatrix jacobian;
    DoubleVector residuals;          
    this->get_jacobian(residuals,jacobian);

    char filename_1[100];
    sprintf(filename_1,"%s/jacobian%i.dat",
            Doc_info.directory().c_str(),
            Doc_info.number());

    char filename_2[100];
    sprintf(filename_2,"%s/residuals%i.dat",
            Doc_info.directory().c_str(),
            Doc_info.number());
    
    jacobian.output(filename_1);
    residuals.output(filename_2); 
    */
}

 
  /// Output contact element situation
 void actions_after_newton_step() { newton_step++; }
 
  /// Update function (empty)
 void actions_before_adapt();        
 /// Update function (empty)
 void actions_after_adapt();

#ifdef REFINE

 /// Access function for the solid mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 

#else

 /// Access function for the solid mesh
 ElasticRectangularQuadMesh<ELEMENT>*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 

#endif
  
 /// Access function to the mesh of surface traction elements
 SolidMesh*& traction_mesh_pt(){return Traction_mesh_pt;} 

 /// Doc the solution
 void doc_solution();

  /// Return the internal pressure that enforces the volume constraint
 double get_interal_pressure();

 ///Run it with increasing pressure
 void run_function(const unsigned nstep, const double p_increment);
 
 //run function to increase the volume 
 void increaseVolume(const double V_target);
 
 /// \short Function to save solution to diretory
 void save_solution();
 void save_solution(const char * directory);
 
 /// Lowers height adaptivly
 void lower_height(double targetHeight);
 void lower_height_w_revert(double targetHeight);

  /// Change a parameter via adaptive continuation
  void change_parameter(double &parameter, double target);
 
 void set_under_relaxation_factor(double ur){under_relaxation_factor = ur;}

  double calc_inflated_vol(double t, double lambda);

  /// Returns global variables in a string
  std::string get_global_variables_as_string();

private:

 /// \short Create contact elements
 /// This will impose the boundary condition on the top
 void create_contact_elements();
 
 void delete_contact_elements();

 /// \short Pass pointer to contact function to the
 /// elements in the contact mesh
 void set_contact_pt();

 /// Update contact height
 void update_height();


 /// \short Pass pointer to traction function to the
 /// elements in the traction mesh
 void set_traction_pt();

 /// \short Create traction elements
 void create_traction_elements();

 /// Delete traction elements
 void delete_traction_elements();
 
 
  /// \short Create volume constraint enforcer elements
 void create_vol_const_elements();

 /// \short Delete volume constraint enforcer elements
 void delete_vol_const_elements();

 /// Trace file
 ofstream Trace_file;
 
 /// Pointers to node whose position we're tracing
 Node* Trace_node_pt;

#ifdef REFINE

 /// Pointer to solid mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>* Solid_mesh_pt;

#else

 /// Pointer to solid mesh
 ElasticRectangularQuadMesh<ELEMENT>* Solid_mesh_pt;

#endif
 
 /// Pointer to mesh of volume constraint elements and the master element
 SolidMesh* Vol_const_mesh_pt;
 Mesh* Vol_const_master_mesh_pt;

 /// Pointer to mesh of traction elements
 SolidMesh* Traction_mesh_pt;
 
/// Pointer to the "surface" mesh
 Mesh* Surface_contact_mesh_pt;
 
 /// ID of contact boundary
 unsigned Contact_boundary_id;

 /// DocInfo object for output
 DocInfo Doc_info;
 
 unsigned newton_step;
 
};


//===========start_of_constructor======================================= 
/// Constructor: 
//====================================================================== 
template<class ELEMENT>
CantileverProblem<ELEMENT>::CantileverProblem(const bool& incompress) 
{
 //Set max number of newton iterations
  Max_newton_iterations = Global_Physical_Variables::MaxNSteps;
 
  // Because the newton method often seems to oscillate, it will reach
  // a large residual and then converge back to a solution. Hence, there
  // is a command line option of setting it to something bigger than the 
  // default 10, say 1000 should give enough room
  Max_residuals = Global_Physical_Variables::MaxResidual;
    
  Newton_solver_tolerance = Global_Physical_Variables::newton_tol;
    
  // Settings for bifurcation tracking, currently not working with under-
  // relaxed newton method
  Desired_newton_iterations_ds = 35; //number of steps in arc length solver?
  
  Minimum_ds = 1e-5; //Minimum ds, no point in going far below that
   
  Scale_arc_length = false;
  
  Bifurcation_detection = false;
    
  // 0.65 appears to be best value: faster than 0.55 and at 0.7 it never converges
  // If the contact pressure is pre-multiplied by 1e4, then 0.95 is optimal for 1 ele in r
  under_relaxation_factor = 0.5; 
  
  newton_step=0; //counter variable, duplicate, clean up 
  
  //Flag to enable 
  //enable_globally_convergent_newton_method();

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

  if(Global_Physical_Variables::printDebugInfo){
    cout << "Now create the mesh " << endl;
  }
 
 #ifdef REFINE

 //Now create the mesh 
 solid_mesh_pt() = new ElasticRefineableRectangularQuadMesh<ELEMENT>(
  n_r,n_theta,l_r,l_theta,origin);
 
 // Set error estimator
 dynamic_cast<ElasticRefineableRectangularQuadMesh<ELEMENT>* >(
  solid_mesh_pt())->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 
#else
 
 //Now create the mesh 
  solid_mesh_pt() = new ElasticRectangularQuadMesh<ELEMENT>(
  n_r,n_theta,l_r,l_theta,origin);

#endif

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
 
#ifdef REFINE
 
 // Refine the mesh uniformly
 dynamic_cast<ElasticRefineableRectangularQuadMesh<ELEMENT>* >(
  solid_mesh_pt())->refine_uniformly();
  
#endif
 
 if(Global_Physical_Variables::enforce_volume_constraint){
    //Creat meshes for surface traction elements
    // that will impose volume constraint and the
    // master mesh that contains a single element
    // that enforces the volume constraint
    Vol_const_mesh_pt = new SolidMesh;
    Vol_const_master_mesh_pt = new Mesh;
    
    
    create_vol_const_elements();     
 }
 else{
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
 }
    

 
 // Create the surface mesh as an empty mesh
 Surface_contact_mesh_pt=new Mesh;


if(Global_Physical_Variables::printDebugInfo){
 	cout << "Going into create_contact_elements()" << endl;
  }

 //set which boundary will be in contact
 Contact_boundary_id = 1;

 // Build 'em 
 create_contact_elements();

 // Set contact pointer
 set_contact_pt();
 
 // Set boundary condition and complete the build of all elements
 //complete_problem_setup();

if(Global_Physical_Variables::printDebugInfo){
 	cout << "Adding meshes together" << endl;
  }



// Either add traction mesh or volume control traction mesh
  if(Global_Physical_Variables::enforce_volume_constraint){
    //add meshes for volume control
    add_sub_mesh(Vol_const_master_mesh_pt);  
    add_sub_mesh(Vol_const_mesh_pt);
  }
  else{
    // Add traction sub-mesh
    add_sub_mesh(traction_mesh_pt());      
  }

 // Add the sub meshes to the problem
 add_sub_mesh(Surface_contact_mesh_pt);

 // Solid mesh is first sub-mesh
 add_sub_mesh(solid_mesh_pt());
  
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

 // Set output directorytraction_mesh_pt()
 #ifdef REFINE
 Doc_info.set_directory("RESLT_REFINE");
#else
 Doc_info.set_directory("RESLT");
#endif


 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",Doc_info.directory().c_str());
 Trace_file.open(filename);
 
  if(Global_Physical_Variables::printDebugInfo){
 	cout << "constructor finished " << endl;
  }

} //end of constructor



//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of traction elements
//========================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::actions_before_adapt()
{
  
 #ifdef REFINE

 if(Global_Physical_Variables::enforce_volume_constraint){  
  // Kill volume constraint elements
  delete_vol_const_elements
}
 else{
  // Kill the traction elements and wipe surface mesh
  delete_traction_elements();
 }  
 
 //delete the contact mesh
 delete_contact_elements();
 
 
 #endif
}// end of actions_before_adapt



//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Rebuild the mesh of traction elements
//========================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::actions_after_adapt()
{
  
  #ifdef REFINE

 
 if(Global_Physical_Variables::enforce_volume_constraint){  
   //creates the surface elements that enforce the volume constraint
    create_vol_const_elements();      
 }
 else{
    create_traction_elements();
    
    set_traction_pt();
 }
 
 // Build 'em 
 create_contact_elements();

 // Set contact pointer
 set_contact_pt();
 
   
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();
 
  #endif
 
}// end of actions_after_adapt

//==================start_of_set_traction_pt==============================
/// Set pointer to contact function for the relevant
/// elements in the contact mesh
//========================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::set_contact_pt()
{
  if(Global_Physical_Variables::printDebugInfo){
 	cout << "Starting  set_contact_pt()" << endl;
  }

   // How many surface elements are in the surface mesh
   unsigned n_element = Surface_contact_mesh_pt->nelement();

   // Loop over the contact elements to pass pointer to penetrator
   //-------------------------------------------------------------
   n_element=Surface_contact_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement 
     AxiSymNonlinearSurfaceContactElement<ELEMENT> *el_pt = 
      dynamic_cast<AxiSymNonlinearSurfaceContactElement<ELEMENT>*>(
       Surface_contact_mesh_pt->element_pt(e));
     
     // Set pointer to penetrator
     el_pt->set_penetrator_pt(Global_Physical_Variables::Penetrator_pt);
    }
 
}// end of set contact pt


//==================start_of_create_contact_elements==============================
/// Create contact surface elements to impose compression
//========================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::create_contact_elements()
  {

   if(Global_Physical_Variables::printDebugInfo){
 	cout << "start_of_create_contact_elements  with boundary " << Contact_boundary_id << endl;
  }

   // How many bulk elements are adjacent to boundary b?
   unsigned b=Contact_boundary_id; 
   unsigned n_element = solid_mesh_pt()->nboundary_element(b);
   bool outputed = false;
   
   // Loop over the bulk elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      solid_mesh_pt()->boundary_element_pt(b,e));
     
     //What is the face index of element e along boundary b
     int face_index = solid_mesh_pt()->face_index_at_boundary(b,e);
     
     // Build the corresponding contact element
     AxiSymNonlinearSurfaceContactElement<ELEMENT>* contact_element_pt = new 
      AxiSymNonlinearSurfaceContactElement<ELEMENT>(bulk_elem_pt,face_index);
      
     // Set the scaling factor for the resiudal
     if(Global_Physical_Variables::contact_pressure_prefactor != 1.0)
       {
	 contact_element_pt->set_pressure_prefactor(Global_Physical_Variables::contact_pressure_prefactor);
       }

     //make sure there is no sticl
      if(contact_element_pt-> is_stick_enabled()){
          contact_element_pt->disable_stick();
     }
     
    if (CommandLineArgs::command_line_flag_has_been_set("--nintpt"))
     {
        if(!outputed){cout << "setting integration schema with " << Global_Physical_Variables::Nintpt << " points." << endl; outputed=true;}
        contact_element_pt->set_integration_scheme(new MyIntegral(Global_Physical_Variables::Nintpt));
     }

    if(Global_Physical_Variables::printDebugInfo){
      std::cout << "Setting contact options now. " << std::endl;
    }

  if (CommandLineArgs::command_line_flag_has_been_set("--contact_use_isoparametric") 
      && Global_Physical_Variables::Contact_use_isoparametric_flag)
     {
       contact_element_pt->set_isoparametric_flag_pt(Global_Physical_Variables::true_flag);
     }

    if (CommandLineArgs::command_line_flag_has_been_set("--contact_use_collocated_penetration")
	&& Global_Physical_Variables::Contact_use_collocated_penetration_flag)
     {
       contact_element_pt->set_collocated_penetration_flag_pt(Global_Physical_Variables::true_flag);
     }
    if (CommandLineArgs::command_line_flag_has_been_set("--contact_use_collocated_contact_pressure")
	&& Global_Physical_Variables::Contact_use_collocated_contact_pressure_flag)
     {
       contact_element_pt->set_collocated_contact_pressure_flag_pt(Global_Physical_Variables::true_flag);
     }

     //Add the contact element to the surface mesh
     Surface_contact_mesh_pt->add_element_pt(contact_element_pt);
     
    } //end of loop over bulk elements adjacent to boundary b
  }

  
//============start_of_delete_contact_elements==============================
/// Delete contact elements and wipe the contact meshes
//=======================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::delete_contact_elements()
{
 // How many surface elements are in the surface mesh
 unsigned n_element = Surface_contact_mesh_pt->nelement();
 
 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  { 
   // Kill surface element
   delete Surface_contact_mesh_pt->element_pt(e);
  }
 
 // Wipe the mesh
 Surface_contact_mesh_pt->flush_element_and_node_storage();

} // end of delete_traction_elements
  
  
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
 // For the quater circle, boundary 1 is the curved top and boundary 3 the inside
 //Traction elements are located on boundary:
 unsigned b=3;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = solid_mesh_pt()->nboundary_element(b);
 if(Global_Physical_Variables::printDebugInfo){
 	cout << "Going to create " << n_element << " elements." << endl;
  }
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    solid_mesh_pt()->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b
   int face_index = solid_mesh_pt()->face_index_at_boundary(b,e);
   
   //Use axisymmetric traction elements
   Traction_mesh_pt->add_element_pt(new AxisymmetricSolidTractionElement<ELEMENT>
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
  { //if(Global_Physical_Variables::printDebugInfo){
 // cout << "Starting " << e << " ..." ;
  //}
   // Kill surface element
   delete Traction_mesh_pt->element_pt(e);
  }
 
 // Wipe the mesh
 Traction_mesh_pt->flush_element_and_node_storage();

} // end of delete_traction_elements

//============start_of_create_vol_const_elements==============================
/// Create traction elements 
//=======================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::create_vol_const_elements()
{
  
 //first create the master element that will add the equation for the
 // volume constraint
  Vol_const_master_mesh_pt->add_element_pt(
    new AxisSymSolidVolumeConstraintElement<ELEMENT>(&Global_Physical_Variables::Volume));
  
   //Cast to a master element - need it to add the surface emsh pointer at
   // the end
   AxisSymSolidVolumeConstraintElement<ELEMENT> *el_pt = 
      dynamic_cast<AxisSymSolidVolumeConstraintElement<ELEMENT>*>
        (Vol_const_master_mesh_pt->element_pt(0));
  
 // For the quater circle, boundary 1 is the curved top and boundary 3 the inside
 //Volume constraint elements are located on boundary:
 unsigned b=3;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = solid_mesh_pt()->nboundary_element(b);
 if(Global_Physical_Variables::printDebugInfo){
    cout << "Going to create " << n_element << " volume constraint elements."; //<< endl;
  }
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    solid_mesh_pt()->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b
   int face_index = solid_mesh_pt()->face_index_at_boundary(b,e);

   //Use axisymmetric traction elements
   Vol_const_mesh_pt->add_element_pt(new AxisymmetricSolidTractionVolumeConstraintElement<ELEMENT>
                                    (bulk_elem_pt,face_index, Vol_const_master_mesh_pt->element_pt(0)));
  }  

 // ensure that the surface mesh gets added as external data to the volume constraint element
 el_pt->set_mesh_pt_and_external(Vol_const_mesh_pt);
 
 
 if(Global_Physical_Variables::printDebugInfo){
    cout << " Vol_const_mesh_pt->nnode() " << Vol_const_mesh_pt->nnode() << "." << endl; 
 //   cout << " Vol_const_mesh_pt->nelements() " << Vol_const_mesh_pt->n_elements() << "."<< endl;
  }
 
} // end of create_vol_const_elements


//============start_of_delete_traction_elements==============================
/// Delete traction elements and wipe the  traction meshes
//=======================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::delete_vol_const_elements()
{
 // How many surface elements are in the surface mesh
 unsigned n_element = Vol_const_mesh_pt->nelement();

 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete Vol_const_mesh_pt->element_pt(e);
  }
 
 // Wipe the mesh
 Vol_const_mesh_pt->flush_element_and_node_storage();
 
 Vol_const_master_mesh_pt->element_pt(0);
 Vol_const_master_mesh_pt->flush_element_and_node_storage();

} // end of delete_vol_const_elements



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
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
        Doc_info.number());
 some_file.open(filename);
 solid_mesh_pt()->output(some_file,n_plot);
 some_file.close();

sprintf(filename,"%s/soln%i_coarse.dat",Doc_info.directory().c_str(),
        Doc_info.number());
 some_file.open(filename);
 solid_mesh_pt()->output(some_file, 2);
 some_file.close();

 // Write trace file: Load/displacement characteristics
 Trace_file << Global_Physical_Variables::P  << " " 
            << Trace_node_pt->x(0) << " " 
            << Trace_node_pt->x(1) << " " 
            << std::endl;
            
// Output traction elements
//  if(Global_Physical_Variables::enforce_volume_constraint){
//     //First output master element
//     dynamic_cast<AxisSymSolidVolumeConstraintElement<ELEMENT>* >(
//        Vol_const_master_mesh_pt->element_pt(0))->output(some_file,n_plot);
// 
//     unsigned nel= Vol_const_mesh_pt->nelement();
//     for (unsigned e=0;e<nel;e++)
//     {
//     dynamic_cast<AxisymmetricSolidTractionVolumeConstraintElement<ELEMENT>* >(
//         Vol_const_mesh_pt->element_pt(e))->output(some_file,n_plot);
//     } 
//  }
//  else{}

  // Output contact elements

 sprintf(filename,"%s/contact%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);

 some_file << "#Contact mesh details for the compression of a sphere." << std::endl;
 some_file <<  get_global_variables_as_string() << std::endl ;
 
 //document contact options
 some_file << "# Contact options are: " << std::endl;
 some_file <<  dynamic_cast<AxiSymNonlinearSurfaceContactElement<ELEMENT>* >(
		    Surface_contact_mesh_pt->element_pt(0))->get_contact_options_in_string()
           << std::endl;

 unsigned nel=Surface_contact_mesh_pt->nelement();

 double f_normal =0;
 double f_orthoganal = 0;
 Vector<double> cont_f(2);
for (unsigned e=0;e<nel;e++)
{
  dynamic_cast<AxiSymNonlinearSurfaceContactElement<ELEMENT>* >(
      Surface_contact_mesh_pt->element_pt(e))->output(some_file,n_plot);
  /*
  force += dynamic_cast<AxiSymNonlinearSurfaceContactElement<ELEMENT>* >(
									 Surface_contact_mesh_pt->element_pt(e))->get_force();
  area += dynamic_cast<AxiSymNonlinearSurfaceContactElement<ELEMENT>* >(
									Surface_contact_mesh_pt->element_pt(e))->get_area();
  */
  dynamic_cast<AxiSymNonlinearSurfaceContactElement<ELEMENT>* >(
		Surface_contact_mesh_pt->element_pt(e))->resulting_contact_force(cont_f);


  ///Set whether or not to use isoparametric basis function for pressure
  //bool* Use_isoparametric_flag_pt;

  ///Set options for basis/test functions for penetration and pressure
  //bool* Use_collocated_penetration_flag_pt;

  ///Set options for basis/test functions for penetration and pressure
  //bool* Use_collocated_contact_pressure_flag_pt;

  f_normal += cont_f[0];
  f_orthoganal += cont_f[1];
}

 std::cout << " Contact force = (" <<f_normal << ", " 
	   << f_orthoganal << ")."  << std::endl;
some_file.close();
 // Increment label for output files
 Doc_info.number()++;
 
 // reset counter
 newton_step=0;
 
 
 double vol=0;
 cout << endl;
  if(Global_Physical_Variables::enforce_volume_constraint){
    
    unsigned n_element=Vol_const_mesh_pt->nelement();
    for(unsigned i=0;i<n_element;i++)
    {
    //Cast to a solid traction element
    AxisymmetricSolidTractionVolumeConstraintElement<ELEMENT> *el_pt = 
        dynamic_cast<AxisymmetricSolidTractionVolumeConstraintElement<ELEMENT>*>
        (Vol_const_mesh_pt->element_pt(i));

    //Set the traction function
    vol += el_pt->get_volume();
    
    //cout << "Element " << i << " has volume " << el_pt->get_volume() << endl;
    }
  }
 else{
    unsigned n_element=traction_mesh_pt()->nelement();
    for(unsigned i=0;i<n_element;i++)
    {
    //Cast to a solid traction element
    AxisymmetricSolidTractionElement<ELEMENT> *el_pt = 
        dynamic_cast<AxisymmetricSolidTractionElement<ELEMENT>*>
        (traction_mesh_pt()->element_pt(i));

    //Set the traction function
    vol += el_pt->get_volume();
    
    //cout << "Element " << i << " has volume " << el_pt->get_volume() << endl;
    }
 }  
  cout << "\n\nTotal Volume is (raw) " << vol << endl;
  
 sprintf(filename,"%s/Volume.txt",Doc_info.directory().c_str());
 some_file.open(filename, std::ios_base::app);
if(Global_Physical_Variables::enforce_volume_constraint){

  AxisSymSolidVolumeConstraintElement<ELEMENT> *el_pt = 
      dynamic_cast<AxisSymSolidVolumeConstraintElement<ELEMENT>*>
        (Vol_const_master_mesh_pt->element_pt(0));
 
  some_file << "H = \t" <<  Global_Physical_Variables::H 
	   << "\t Volume = \t" << vol 
	    << "\t Pressure = \t" << el_pt->internal_data_pt(0)->value(0) 
	    << "\t Normal Contact Force = \t" << f_normal<< std:: endl;
 }
 some_file.close();
 
//   std::ofstream myfile;
//    myfile.open ("traction_elements_debug.txt", std::ios::app );   
//    myfile << "\n\n\n############################## \nFinished Newton Solve \n############################## " << Doc_info.number() -1 << std::endl;
//    myfile.close();

} //end doc

template<class ELEMENT>
double CantileverProblem<ELEMENT>::get_interal_pressure(){
  AxisSymSolidVolumeConstraintElement<ELEMENT> *el_pt = 
      dynamic_cast<AxisSymSolidVolumeConstraintElement<ELEMENT>*>
        (Vol_const_master_mesh_pt->element_pt(0));
 
  return el_pt->internal_data_pt(0)->value(0);
  
}


template<class ELEMENT>
void CantileverProblem<ELEMENT>::run_function(const unsigned nstep, const double p_increment){

 Global_Physical_Variables::P =0.0;
 //solve once in undeformed configuration, hopefully this should converge immediatly
 #ifdef REFINE
    newton_solve(Global_Physical_Variables::max_adapt);
  #else
    newton_solve();
  #endif
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
   #ifdef REFINE
    newton_solve(Global_Physical_Variables::max_adapt);
  #else
    newton_solve();
  #endif
   // Doc solution
   doc_solution();

  }

}


template<class ELEMENT>
void CantileverProblem<ELEMENT>::increaseVolume(const double V_target){

 double increment = V_target - Global_Physical_Variables::Volume;
 newton_solve(); //Solve once, should return immediatly
 
  //Save current solution
  char filename[100];
  sprintf(filename,"sol_adaptive_increase_volume.dat");
  ofstream Solution_output_file;
  ifstream Solution_input_file;

  Solution_output_file.open(filename);
  dump(Solution_output_file);    
  Solution_output_file.close();


 unsigned steps = 0;

 while(Global_Physical_Variables::Volume < V_target){
   Global_Physical_Variables::Volume += increment;
   
   cout << endl;
   cout << "Step " << steps <<  
        ". Incrementing Volume to " << Global_Physical_Variables::Volume << 
        " from " << Global_Physical_Variables::Volume  - increment << endl;
   try{
   // Solve the problem    
   newton_solve();

   // Doc solution
   doc_solution();
   
   //save solution to file so it can be read in case
   // next step fails
   std::cout << "Saving solution." << std::endl;
   Solution_output_file.open(filename);
   dump(Solution_output_file);    
   Solution_output_file.close();

   }
   catch(...){
     std::cout << "Decreasing increment." << std::endl;

     Global_Physical_Variables::Volume -= increment;
     increment = increment/2.0;
     
     if(abs(increment) < 1e-5){
       std::cout << "Increment too small at " << increment << ". Returning. " << std::endl;
       return;
     }

    //load previous soltution
    Solution_input_file.open(filename);
    read(Solution_input_file);
    Solution_input_file.close();

   }
   steps++;
 }
 /* double change = V_target - Global_Physical_Variables::Volume ;
 unsigned nstep = (int)( change/increment);
 for(unsigned i=0;i<nstep;i++)
  {
   // Increment Volume
   Global_Physical_Variables::Volume += increment;
    
   cout << endl;
   cout << "Step " << i << " out of " << nstep << 
        ". Incrementing Volume to " << Global_Physical_Variables::Volume << 
        " from " << Global_Physical_Variables::Volume  - increment << endl;
   // Solve the problem    
   newton_solve();

   // Doc solution
   doc_solution();

   }*/
  Global_Physical_Variables::Volume = V_target;
  newton_solve();
  doc_solution();
  save_solution(); //save solution so it can be read from file in futur
}


template<class ELEMENT>
void CantileverProblem<ELEMENT>::save_solution(){
  save_solution(Doc_info.directory().c_str());
}

template<class ELEMENT>
void CantileverProblem<ELEMENT>::save_solution(const char * directory){

  std::string globals = get_global_variables_as_string();
    char filename[200];
    sprintf(filename,"%s/sol_%s.dat",
            directory,
	    globals.c_str());
    
    ofstream Solution_output_file;
    Solution_output_file.open(filename);
    Solution_output_file.precision(17);
    Solution_output_file.setf(ios::fixed);
    Solution_output_file.setf(ios::showpoint);
    //write solution to file
    dump(Solution_output_file);
    
    Solution_output_file.close();
}


template<class ELEMENT>
/// Function to lower the parameter adaptivly 
void CantileverProblem<ELEMENT>::change_parameter(double &parameter, double target){
  double ds;
  
  int kk = 0; //number of steps
  int nr_errors=0; //Record the number of times the setp-size had to be halfed
  int steps_since_reset = 0;  
  //Save current solution
  char filename[100];
  sprintf(filename,"sol_adaptive_change_parameter.dat");
  ofstream Solution_output_file;
  ifstream Solution_input_file;

  //Assumes problem is currently solved
  // Let's solve once more just to be certain
  std::cout << " Initial Newton solve, should converge immediatly." << std::endl;
  #ifdef REFINE
     newton_solve(Global_Physical_Variables::max_adapt);
  #else
     newton_solve();
  #endif
  
  Solution_output_file.open(filename);
  Solution_output_file.precision(17);
  Solution_output_file.setf(ios::fixed);
  Solution_output_file.setf(ios::showpoint);
  dump(Solution_output_file);    
  Solution_output_file.close();

  if(parameter > target)
    {
	  ds = -1e-2;
    }
  else
    {
      ds = 1e-2;
    }

  // Make sure initial step isn't too big
  if(std::abs(parameter - target) < ds)
    {
      ds = target - parameter;
    }

  while(abs(parameter -  target) > 1e-4)
    {
      if(ds < 0 && parameter+ds < target){
	break;
      }
      if(ds > 0 && parameter+ds > target){
	break;
      }

      parameter+=ds;
          
      //Update Volume in case the parameter is lambda
      Global_Physical_Variables::Volume = 
	calc_inflated_vol(Global_Physical_Variables::t, Global_Physical_Variables::lambda);
      std::cout << std::endl;
      std::cout << "Step " << kk << 
            ". Decrementing to " << parameter
		<< " with step of " << ds << std::endl;
     std::cout << "H = " << Global_Physical_Variables::H
	       << " lambda = " << Global_Physical_Variables::lambda
	       << " Volume = " << Global_Physical_Variables::Volume
	       << " P = " << get_interal_pressure() << std::endl;
      try{
        // Solve the problem    
        #ifdef REFINE
          newton_solve(Global_Physical_Variables::max_adapt);
        #else
          newton_solve();
        #endif
          
        //increase step size iff the number of newton steps is 
        //smaller than 10 (this assuming we are using an under-relaxed
        //Newton method)
	if(newton_step < 12 && steps_since_reset > 3 && std::abs(2*ds) < std::abs(parameter - target))
	    {
	      std::cout << "Doubeling step size from " 
			<< ds << " to " << 2*ds 
			<< ". Current Parameter = " << parameter
			<< " and target = " << target << "." << std::endl; 
	      ds = ds*2.0;
	    }
        
        // Doc solution
        doc_solution();

        Solution_output_file.open(filename);
	Solution_output_file.precision(17);
	Solution_output_file.setf(ios::fixed);
	Solution_output_file.setf(ios::showpoint);
        dump(Solution_output_file);    
        Solution_output_file.close();
        
        kk++;
        steps_since_reset++;
      }
      catch(...){//Didn't working
        //Revert parameter
        parameter-=ds;
        
        //half step-size
        ds = ds/2.0;
        
        steps_since_reset = 0;
        
        if(abs(ds) < 1e-5){
          std::cout << "Step is too small, aborting!" << std::endl;
          return;
        }
        
        //load previous soltution
        Solution_input_file.open(filename);
        read(Solution_input_file);
        Solution_input_file.close();
        
        nr_errors++;
      }
    } 
  //Save previous height
  double final = parameter;

  //Reach below target height
  parameter = target;
  std::cout << "Setting height to " << target << 
    " and solving one last time." << std::endl; 
  try{
    // Solve the problem    
    #ifdef REFINE
    newton_solve(Global_Physical_Variables::max_adapt) ;
    #else
    newton_solve();
    #endif  // Doc solution
    doc_solution();
  }
  catch(...){
    std::cout << "LAst Newton solve failed, reverting to previous solution and returning" 
	      << std::endl;
    //Revert height
    parameter = final;
        
    //load previous soltution
    Solution_input_file.open(filename);
    read(Solution_input_file);
    Solution_input_file.close();

    try{
      // Solve the problem at previous state, should converge immediatly
      #ifdef REFINE
        newton_solve(Global_Physical_Variables::max_adapt) ;
      #else
        newton_solve();
      #endif  // Doc solution
      doc_solution();
    }
    catch(...){
      nr_errors++;
    }
    nr_errors++;
  }

  std::cout << "The 'lower_parameter' function took " << kk +1 << " steps and" <<
    " had to lower the step-size " << nr_errors << " times. " << std::endl;
}

template<class ELEMENT>
/// Function to lower the height H adaptivly to reach states of high compression. 
void CantileverProblem<ELEMENT>::lower_height(double targetHeight){
  
  //double diff = Global_Physical_Variables::H - targetHeight;
  double ds = -1e-2; 
  
  int kk = 0; //number of steps
  int nr_errors=0; //Record the number of times the setp-size had to be halfed
  int steps_since_reset = 0;  
  //Save current solution
  char filename[100];
  sprintf(filename,"sol_adaptive_lower_height.dat");
  ofstream Solution_output_file;
  ifstream Solution_input_file;

  //Assumes problem is currently solved
  // Let's solve once more just to be certain
  std::cout << " Initial Newton solve, should converge immediatly." << std::endl;
  #ifdef REFINE
     newton_solve(Global_Physical_Variables::max_adapt);
  #else
     newton_solve();
  #endif
  
  Solution_output_file.open(filename);
  dump(Solution_output_file);    
  Solution_output_file.close();
  
  while(Global_Physical_Variables::H > targetHeight)
    {
      //Lower contact height
      if(Global_Physical_Variables::H+ds < targetHeight){
	break;
      }
      Global_Physical_Variables::H+=ds;
          
      cout << endl;
      cout << "Step " << kk << 
            ". Decrementing Height to " << Global_Physical_Variables::H 
            << " with step of " << ds << endl;
      
      try{
        // Solve the problem    
        #ifdef REFINE
          newton_solve(Global_Physical_Variables::max_adapt);
        #else
          newton_solve();
        #endif
          
        //increase step size iff the number of newton steps is 
        //smaller than 10 (this assuming we are using an under-relaxed
        //Newton method)
        if(newton_step < 12 && steps_since_reset > 3){ds = ds*2.0;}
        
        // Doc solution
        doc_solution();

        Solution_output_file.open(filename);
        dump(Solution_output_file);    
        Solution_output_file.close();
        
        kk++;
        steps_since_reset++;
      }
      catch(...){//Didn't working
        //Revert height
        Global_Physical_Variables::H-=ds;
        
        //half step-size
        ds = ds/2.0;
        
        steps_since_reset = 0;
        
        if(abs(ds) < 1e-5){
          std::cout << "Step is too small, aborting!" << std::endl;
          return;
        }
        
        //load previous soltution
        Solution_input_file.open(filename);
        read(Solution_input_file);
        Solution_input_file.close();
        
        nr_errors++;
      }
    } 
  //Save previous height
  double finalH = Global_Physical_Variables::H;

  //Reach below target height
  Global_Physical_Variables::H = targetHeight;
  std::cout << "Setting height to " << targetHeight << 
    " and solving one last time." << std::endl; 
  try{
    // Solve the problem    
    #ifdef REFINE
    newton_solve(Global_Physical_Variables::max_adapt) ;
    #else
    newton_solve();
    #endif  // Doc solution
    doc_solution();
  }
  catch(...){
    std::cout << "LAst Newton solve failed, reverting to previous solution and returning" 
	      << std::endl;
    //Revert height
    Global_Physical_Variables::H = finalH;
        
    //load previous soltution
    Solution_input_file.open(filename);
    read(Solution_input_file);
    Solution_input_file.close();

    try{
      // Solve the problem at previous state, should converge immediatly
      #ifdef REFINE
        newton_solve(Global_Physical_Variables::max_adapt) ;
      #else
        newton_solve();
      #endif  // Doc solution
      doc_solution();
    }
    catch(...){
      nr_errors++;
    }
    nr_errors++;
  }

  std::cout << "The 'lower_height' function took " << kk +1 << " steps and" <<
    " had to lower the step-size " << nr_errors << " times. " << std::endl;
}

template<class ELEMENT>
/// Function to lower the height H adaptivly to reach states of high compression. 
void CantileverProblem<ELEMENT>::lower_height_w_revert(double targetHeight){
  
  //double diff = Global_Physical_Variables::H - targetHeight;
  double ds = -5e-2; 
  
  int kk = 0; //number of steps
  int nr_errors=0; //Record the number of times the setp-size had to be halfed
  
  //Save current solution
  char filename1[100];
  sprintf(filename1,"sol_adaptive_lower_height_1.dat");
  char filename2[100];
  sprintf(filename2,"sol_adaptive_lower_height_2.dat");
  char filename3[100];
  sprintf(filename3,"sol_adaptive_lower_height_3.dat");
  
  ofstream Solution_output_file;
  ifstream Solution_input_file;


  Solution_output_file.open(filename1);
  dump(Solution_output_file);    
  Solution_output_file.close();
  
  while(Global_Physical_Variables::H > targetHeight)
    {
      //Lower contact height
      Global_Physical_Variables::H+=ds;
          
      cout << endl;
      cout << "Step " << kk << 
            ". Decrementing Height to " << Global_Physical_Variables::H 
            << " with step of " << ds << endl;
      
      try{
        // Solve the problem    
        #ifdef REFINE
          newton_solve(Global_Physical_Variables::max_adapt);
        #else
          newton_solve();
        #endif
          
        //increase step size iff the number of newton steps is 
        //smaller than 10 (this assuming we are using an under-relaxed
        //Newton method)
        if(newton_step < 12){ds = ds*2.0;}
        
        // Doc solution
        doc_solution();

        //Copy filed to other files
        if(kk>1){
          std::ifstream  src(filename2, std::ios::binary);
          std::ofstream  dst(filename3,   std::ios::binary);
          dst << src.rdbuf();
        }
        std::ifstream  src(filename1, std::ios::binary);
        std::ofstream  dst(filename2,   std::ios::binary);
        dst << src.rdbuf();
        
        Solution_output_file.open(filename1);
	Solution_output_file.precision(17);
	Solution_output_file.setf(ios::fixed);
	Solution_output_file.setf(ios::showpoint);
        dump(Solution_output_file);    
        Solution_output_file.close();
        
        kk++;
      }
      catch(...){//Didn't working        
        //load previous soltution
        if(kk > 2){
          //Revert height
          Global_Physical_Variables::H-= 3*ds;
          Solution_input_file.open(filename3);
        }
        else{
          if(kk==2){
            //Revert height
            Global_Physical_Variables::H-=2*ds;
            std::cout << "Reading solution one step back." << std::endl;
            Solution_input_file.open(filename2);
          }
          else{
            //Revert height
            Global_Physical_Variables::H-=ds;
            std::cout << "Reading previous solution." << std::endl;
            Solution_input_file.open(filename1);
          }
        }

         //half step-size
        ds = ds/2.0;
        
        if(abs(ds) < 1e-5){
          std::cout << "Step is too small, aborting!" << std::endl;
          return;
        }
        
        read(Solution_input_file);
        Solution_input_file.close();
        
        nr_errors++;
      }
    } 
  
  //Reach below target height
  Global_Physical_Variables::H = targetHeight;
  std::cout << "Setting height to " << targetHeight << 
    " and solving one last time." << std::endl; 
  // Solve the problem    
  #ifdef REFINE
    newton_solve(Global_Physical_Variables::max_adapt);
  #else
    newton_solve();
  #endif  // Doc solution
  doc_solution();
  
  std::cout << "The 'lower_height' function took " << kk +1 << " steps and" <<
    " had to lower the step-size " << nr_errors << " times. " << std::endl;
}

/// Calcuates the relevant volume based on tickness t and desiered 
/// inflantion lambda, assuming that the volume of the elastic material
/// remains constant
template<class ELEMENT>
/// Function to lower the height H adaptivly to reach states of high compression. 
double CantileverProblem<ELEMENT>::calc_inflated_vol(double t, double lambda){
  //initial volume under the elastic solid
  double v0 = pow((1 - Global_Physical_Variables::t), 3) / 3;
  
  //volume of elastic material
  double ve0 = 1/3.0 - v0;
  
  //asuming that ve = ve0
  double inner_radius = pow( pow(lambda, 3) - 3.0*ve0 ,1/3.0);
  
  double v = pow(inner_radius,3) / 3.0;
  
  return v;
}


template<class ELEMENT>
/// Function that returns string with all global variables
std::string CantileverProblem<ELEMENT>::get_global_variables_as_string()
{
    char constitutivelaw[50];
    if(!Global_Physical_Variables::constitutive_law.compare("MR"))
      {
	sprintf(constitutivelaw,"CL=%s_C1=%.4f_C2=%.4f",
	    Global_Physical_Variables::constitutive_law.c_str(),
            Global_Physical_Variables::C1,
            Global_Physical_Variables::C2);
    	
      }
    else{
      if(!Global_Physical_Variables::constitutive_law.compare("GH"))
	{
	  sprintf(constitutivelaw,"CL=%s_E=%.4f_Nu=%.4f",
            Global_Physical_Variables::constitutive_law.c_str(),
            Global_Physical_Variables::E,
            Global_Physical_Variables::Nu);    	
      }
    }


    char buff[200];
    sprintf(buff,"nreler=%d_nreletheta=%d_H=%.4f_lambda=%.2f_vol=%.4f_t=%.2f_%s",
            Global_Physical_Variables::n_ele_r,
            Global_Physical_Variables::n_ele_theta,
            Global_Physical_Variables::H,
            Global_Physical_Variables::lambda,
            Global_Physical_Variables::Volume,
            Global_Physical_Variables::t,
	    constitutivelaw);

  std::string globals = buff;
  return globals;
}



//====================================================================
//####################################################################
//====================================================================


//===================================================================
/// Function to lower a parameter passed by reference
//===================================================================

void change_parameter(double &parameter, double target){
 // Initial values for parameter values
 Global_Physical_Variables::P=0.0; 
 Global_Physical_Variables::Gravity=0.0;

 // Create penetrator
 Global_Physical_Variables::Penetrator_pt =
  new AxiSymPenetrator(&Global_Physical_Variables::H);

 //The number of elements such that after so they stay approximatlly
 // square even after some streching 
 Global_Physical_Variables::n_ele_theta =(int) (1.571/(Global_Physical_Variables::t/   Global_Physical_Variables::n_ele_r)  + 0.5)*3;


 cout << "Thickness = " << Global_Physical_Variables::t << " with " <<
        Global_Physical_Variables::n_ele_r << ", " << 
        Global_Physical_Variables::n_ele_theta << 
        " elements in the r and theta direction" << endl;

bool incompressible = true;
#ifdef REFINE
      //TODO: Write refinable element
#else
  CantileverProblem<AxisymQPVDElementWithPressure> problem2(incompressible);
#endif
 cout << endl;

 //  std::cout << "Press Enter\n";
 //  std::cin.ignore(); 
  
 ifstream Solution_input_file; //for loading solution
 

 if(!Global_Physical_Variables::loadSol.empty()){
  
  ofstream Solution_output_file;
  ifstream Solution_input_file;

  //load previous soltution
  std::cout << "Loading solution: " << Global_Physical_Variables::loadSol << std::endl;
  Solution_input_file.open(Global_Physical_Variables::loadSol.c_str());
  problem2.read(Solution_input_file);
  Solution_input_file.close();

   //Lets check how ti looks
   problem2.doc_solution();

   problem2.set_under_relaxation_factor(0.4);  
   Global_Physical_Variables::Volume = problem2.calc_inflated_vol(Global_Physical_Variables::t, Global_Physical_Variables::lambda);
   std::cout << "Initial Newton Solve, shoudl converge immediatly. "
	     << "H = " << Global_Physical_Variables::H
	     << " Volume = " << Global_Physical_Variables::Volume
	     << " P = " << problem2.get_interal_pressure() << std::endl;
   //solve once in undeformed configuration, hopefully this should converge immediatly
   #ifdef REFINE
     Global_Physical_Variables::max_adapt=1;
     problem2.newton_solve(Global_Physical_Variables::max_adapt);
   #else
     problem2.newton_solve();
   #endif

   //Lets check that displacement is zero
   problem2.doc_solution();

}
 else{
   std::cout << "No loadSol variable provided, cannot load solution. " 
	     << "Will try a newton solve with current parameters. This might fail. " 
	     << std::endl; 
 }

   double stepsize = 0.02; 
   double current_Target = parameter;

   if(parameter > target)
     {
       stepsize = stepsize*-1.0;
       while(parameter > target){
	 std::cout << "Current parameter value = " << parameter
		   << " Target = " << target
		   << " Target for next solution = " << current_Target
		   << std::endl;
	 current_Target += stepsize;
	 problem2.change_parameter(parameter, current_Target);
	 problem2.save_solution();//to have something to restart from in futur
	 problem2.save_solution(".");
	 std::cout << "Current parameter = " << parameter 
		   << " Target parameter = " << target
		   << std::endl;
       }
     }
   else
     {
       while(parameter < target){
	 std::cout << "Current parameter value = " << parameter
		   << " Target = " << target
		   << " Target for next solution = " << current_Target
		   << std::endl;
	 current_Target += stepsize;
	 problem2.change_parameter(parameter, current_Target);
	 problem2.save_solution();//to have something to restart from in futur
	 problem2.save_solution(".");
	 std::cout << "Current parameter = " << parameter 
		   << " Target parameter = " << target
		   << std::endl;
       }

     }
   

     
}


//===================================================================
/// Function to lower the C1 variable of the Moony-Rivling 
/// consitutive law
//===================================================================

void lower_C1(){
 // Set the consitutive law via command line argument

 if(Global_Physical_Variables::constitutive_law == "GH"){
   std::cout << "Trying to set constitutive law to GeneralisedHookean. " 
	     << "'lower_C1' only works with Moony-Rivlin law! "
	     << "Function will return without doing anything"
	     << std::endl;

   return;
 }

 // Initial values for parameter values
 Global_Physical_Variables::P=0.0; 
 Global_Physical_Variables::Gravity=0.0;

 // Create penetrator
 Global_Physical_Variables::Penetrator_pt =
  new AxiSymPenetrator(&Global_Physical_Variables::H);
  
 Global_Physical_Variables::n_ele_theta =(int) (1.571/(Global_Physical_Variables::t/   Global_Physical_Variables::n_ele_r)  + 0.5)*3;


 cout << "Thickness = " << Global_Physical_Variables::t << " with " <<
        Global_Physical_Variables::n_ele_r << ", " << 
        Global_Physical_Variables::n_ele_theta << 
        " elements in the r and theta direction" << endl;

bool incompressible = true;
#ifdef REFINE
      //TODO: Write refinable element
#else
  CantileverProblem<AxisymQPVDElementWithPressure> problem2(incompressible);
#endif
 cout << endl;

 //  std::cout << "Press Enter\n";
 //  std::cin.ignore(); 
  
 ifstream Solution_input_file; //for loading solution
 

 if(!Global_Physical_Variables::loadSol.empty()){
  
  ofstream Solution_output_file;
  ifstream Solution_input_file;

  //load previous soltution
  std::cout << "Loading solution: " << Global_Physical_Variables::loadSol << std::endl;
  Solution_input_file.open(Global_Physical_Variables::loadSol.c_str());
  problem2.read(Solution_input_file);
  Solution_input_file.close();

   //Lets check how ti looks
   problem2.doc_solution();

   problem2.set_under_relaxation_factor(0.4);  
   Global_Physical_Variables::Volume = problem2.calc_inflated_vol(Global_Physical_Variables::t, Global_Physical_Variables::lambda);
   std::cout << "Initial Newton Solve, shoudl converge immediatly. "
	     << "H = " << Global_Physical_Variables::H
	     << " Volume = " << Global_Physical_Variables::Volume
	     << " P = " << problem2.get_interal_pressure() << std::endl;
   //solve once in undeformed configuration, hopefully this should converge immediatly
   #ifdef REFINE
     Global_Physical_Variables::max_adapt=1;
     problem2.newton_solve(Global_Physical_Variables::max_adapt);
   #else
     problem2.newton_solve();
   #endif

   //Lets check that displacement is zero
   problem2.doc_solution();

   double stepsize = 0.02; 
   double target = 0.1;
   double current_Target = Global_Physical_Variables::C1;
   while(Global_Physical_Variables::C1 > target){
     std::cout << "Current C1 = " << Global_Physical_Variables::C1 
	       << " Target height = " << target
	       << " Target for next solution = " << current_Target
	       << std::endl;
     current_Target -= stepsize;
     problem2.change_parameter(Global_Physical_Variables::C1, current_Target);
     problem2.save_solution();//to have something to restart from in futur
     problem2.save_solution(".");
     std::cout << "Current C1 = " << Global_Physical_Variables::C1 
	       << " Target C1 = " << target
	       << std::endl;
   }
   
}
 else{
   std::cout << "No loadSol variable provided, cannot load solution. " 
	     << "'lower_C1' function will now return without doing anything." 
	     << std::endl; 
 }
     
}


//===================================================================
/// Function for standard run
//===================================================================
void standard_run(double step_increment){
 // Set the consitutive law via command line argument

 // Initial values for parameter values
 Global_Physical_Variables::P=0.0; 
 Global_Physical_Variables::Gravity=0.0;

 // Create penetrator
 Global_Physical_Variables::Penetrator_pt =
  new AxiSymPenetrator(&Global_Physical_Variables::H);
  
if(false){
 Global_Physical_Variables::n_ele_theta = 2.0* Global_Physical_Variables::n_ele_r;
}
else{
 Global_Physical_Variables::n_ele_theta =(int) (1.571/(Global_Physical_Variables::t/   Global_Physical_Variables::n_ele_r)  + 0.5)*3;
}

 cout << "Thickness = " << Global_Physical_Variables::t << " with " <<
        Global_Physical_Variables::n_ele_r << ", " << 
        Global_Physical_Variables::n_ele_theta << 
        " elements in the r and theta direction" << endl;

bool incompressible = true;
#ifdef REFINE
      //TODO: Write refinable element
#else
  CantileverProblem<AxisymQPVDElementWithPressure> problem2(incompressible);
#endif
 cout << endl;

 //  std::cout << "Press Enter\n";
 //  std::cin.ignore(); 
  
 ifstream Solution_input_file; //for loading solution
 

 if(!Global_Physical_Variables::loadSol.empty()){
   //char filename[100];
   //sprintf(filename,"RESLT/sol_nreler=1_nreletheta=32_H=0.6000_vol=0.8243_t=0.10.dat");
  
  ofstream Solution_output_file;
  ifstream Solution_input_file;

  //load previous soltution
  std::cout << "Loading solution: " << Global_Physical_Variables::loadSol << std::endl;
  Solution_input_file.open(Global_Physical_Variables::loadSol.c_str());
  problem2.read(Solution_input_file);
  Solution_input_file.close();

   //Lets check how ti looks
   problem2.doc_solution();

   problem2.set_under_relaxation_factor(0.4);  
   Global_Physical_Variables::Volume = problem2.calc_inflated_vol(Global_Physical_Variables::t, Global_Physical_Variables::lambda);
   std::cout << "Initial Newton Solve, shoudl converge immediatly. "
	     << "H = " << Global_Physical_Variables::H
	     << " Volume = " << Global_Physical_Variables::Volume
	     << " P = " << problem2.get_interal_pressure() << std::endl;
   //solve once in undeformed configuration, hopefully this should converge immediatly
   #ifdef REFINE
     Global_Physical_Variables::max_adapt=1;
     problem2.newton_solve(Global_Physical_Variables::max_adapt);
   #else
     problem2.newton_solve();
   #endif

   //Lets check that displacement is zero
   problem2.doc_solution();

   double H_stepsize = step_increment; //0.02 * 0.6 * Global_Physical_Variables::lambda;
   double currentH_Target = Global_Physical_Variables::H;
   while(Global_Physical_Variables::H > 0.39*Global_Physical_Variables::lambda){
     std::cout << "Current height = " << Global_Physical_Variables::H 
	       << " Target height = " << 0.39*Global_Physical_Variables::lambda
	       << " Target for next solution = " << currentH_Target
	       << std::endl;
     currentH_Target -= H_stepsize;
     problem2.change_parameter(Global_Physical_Variables::H, currentH_Target);
     problem2.save_solution();//to have something to restart from in futur
     problem2.save_solution(".");
     std::cout << "Current height = " << Global_Physical_Variables::H 
	       << " Target height = " << 0.39*Global_Physical_Variables::lambda
	       << std::endl;
   }
   
}
 else{

   Global_Physical_Variables::P =0.0;

   //Set height to 0.01 aboe the predicted height of capsule
   Global_Physical_Variables::H = Global_Physical_Variables::lambda; 

   // save lambda
   double  save_lambda = Global_Physical_Variables::lambda; 
 
   Global_Physical_Variables::lambda = 1.0;

 cout << "Starting Newton Solve " << endl;
   //solve once in undeformed configuration, hopefully this should converge immediatly
 #ifdef REFINE
 Global_Physical_Variables::max_adapt=1;
   problem2.newton_solve(Global_Physical_Variables::max_adapt);
#else
   problem2.newton_solve();
#endif
 
   //Lets check that displacement is zero
   problem2.doc_solution();

   //Increase volume first
   //run the function that increases volume
   problem2.set_under_relaxation_factor(1.0); //when there is no contact, can use normal newton solve

   Global_Physical_Variables::lambda = save_lambda;

   std::cout << " Aiming for volume of " << problem2.calc_inflated_vol(Global_Physical_Variables::t, Global_Physical_Variables::lambda) <<
     " corrseponding to lambda of " << Global_Physical_Variables::lambda << " from a volume of " << Global_Physical_Variables::Volume <<
     " . Contact height is " << Global_Physical_Variables::H << "." << std::endl;

   Global_Physical_Variables::lambda  = 1.0;
   problem2.change_parameter(Global_Physical_Variables::lambda, save_lambda);

   problem2.set_under_relaxation_factor(0.4);  //max under-relaxation possible for 2 elemtne in r direction    

   double H_stepsize = step_increment; //0.02 * 0.6 * Global_Physical_Variables::lambda;
   double currentH_Target = Global_Physical_Variables::H -0.1;
   while(Global_Physical_Variables::H > 0.39*Global_Physical_Variables::lambda){
     std::cout << "Current height = " << Global_Physical_Variables::H 
	       << " Target height = " << 0.39*Global_Physical_Variables::lambda
	       << " Target for next solution = " << currentH_Target
	       << std::endl;
     currentH_Target -= H_stepsize;
     problem2.lower_height(currentH_Target);
     problem2.save_solution();//to have something to restart from in futur
     problem2.save_solution(".");
     std::cout << "Current height = " << Global_Physical_Variables::H 
	       << " Target height = " << 0.39*Global_Physical_Variables::lambda
	       << std::endl;
     
   }
}
}


void standard()
{
  standard_run(0.2);
}


//=====================================================================
// Function to demonstrate problem with reloading solution from file
//====================================================================

void reload_solution_test(){


 if(Global_Physical_Variables::constitutive_law == "GH"){
   //Create generalised Hookean constitutive equations
   Global_Physical_Variables::Constitutive_law_pt = 
     new GeneralisedHookean(&Global_Physical_Variables::Nu,
                        &Global_Physical_Variables::E);

   std::cout << "Setting constitutive law to GeneralisedHookean" << std::endl;
 }

 if(Global_Physical_Variables::constitutive_law == "MR"){
  // Create MooneyRivlin constitutive equations
  Global_Physical_Variables::Strain_energy_function_pt = 
  new MooneyRivlin(&Global_Physical_Variables::C1,
                        &Global_Physical_Variables::C2);

  // Define a constitutive law (based on strAxisymmetricSolidTractionElementain energy function)
  Global_Physical_Variables::Constitutive_law_pt = 
  new IsotropicStrainEnergyFunctionConstitutiveLaw(
   Global_Physical_Variables::Strain_energy_function_pt);

   std::cout << "Setting constitutive law to MooneyRivlin" << std::endl;
 }

 // Initial values for parameter values
 Global_Physical_Variables::P=0.0; 
 Global_Physical_Variables::Gravity=0.0;


  if(Global_Physical_Variables::printDebugInfo){
        cout << " " << endl;
 	cout << "Initial Newton solve: " << endl;
 }

 // Create penetrator
 Global_Physical_Variables::Penetrator_pt =
  new AxiSymPenetrator(&Global_Physical_Variables::H);
  
if(false){
 Global_Physical_Variables::n_ele_theta = 2.0* Global_Physical_Variables::n_ele_r;
}
else{
 Global_Physical_Variables::n_ele_theta =(int) (1.571/(Global_Physical_Variables::t/   Global_Physical_Variables::n_ele_r)  + 0.5)*3;
}

 cout << "Thickness = " << Global_Physical_Variables::t << " with " <<
        Global_Physical_Variables::n_ele_r << ", " << 
        Global_Physical_Variables::n_ele_theta << 
        " elements in the r and theta direction" << endl;

 bool incompressible = true;

 // Creating two problems
 CantileverProblem<AxisymQPVDElementWithPressure> problem1(incompressible); 

 cout << endl;

   // Set volume to correspond to the actual volume of the undefomred sphere
   Global_Physical_Variables::Volume = pow((1 - Global_Physical_Variables::t), 3) / 3;

   //Set height to 0.01 aboe the predicted height of capsule
   Global_Physical_Variables::H = Global_Physical_Variables::lambda; 


 cout << "Starting Newton Solve " << endl;
   //solve once in undeformed configuration, hopefully this should converge immediatly
   problem1.newton_solve();

   //Lets check that displacement is zero
   problem1.doc_solution();

   problem1.set_under_relaxation_factor(0.4);  //max under-relaxation possible for 2 elemtne in r direction    

   double finalTarget = 0.8;
   double H_stepsize = 0.1; //0.02 * 0.6 * Global_Physical_Variables::lambda;
   double currentH_Target = Global_Physical_Variables::H;
   while(Global_Physical_Variables::H > finalTarget*Global_Physical_Variables::lambda){
     std::cout << "Current height = " << Global_Physical_Variables::H 
	       << " Target height = " << finalTarget*Global_Physical_Variables::lambda
	       << " Target for next solution = " << currentH_Target
	       << std::endl;
     currentH_Target -= H_stepsize;
     problem1.lower_height(currentH_Target);
     problem1.save_solution(".");
     std::cout << "Current height = " << Global_Physical_Variables::H 
	       << " Target height = " << finalTarget*Global_Physical_Variables::lambda
	       << std::endl;
     
   }
  
  ofstream Solution_output_file;
  ifstream Solution_input_file;
  CantileverProblem<AxisymQPVDElementWithPressure> problem2(incompressible);
 
 //load previous soltution
  std::cout << "Loading solution: " << Global_Physical_Variables::loadSol << std::endl;
  std::string file = "sol_nreler=1_nreletheta=48_H=0.8000_vol=0.2430_t=0.10.dat";
  Solution_input_file.open(file.c_str());
  problem2.read(Solution_input_file);
  Solution_input_file.close();

   //Lets check how it looks
   problem2.doc_solution();

   problem2.set_under_relaxation_factor(0.4);  
   Global_Physical_Variables::Volume = problem2.calc_inflated_vol(Global_Physical_Variables::t, Global_Physical_Variables::lambda);
   std::cout << "Initial Newton Solve, shoudl converge immediatly. "
	     << "H = " << Global_Physical_Variables::H
	     << " Volume = " << Global_Physical_Variables::Volume
	     << " P = " << problem2.get_interal_pressure() << std::endl;
   //solve once in undeformed configuration, hopefully this should converge immediatly
   problem2.newton_solve();

   //Lets check that displacement is zero
   problem2.doc_solution();
   finalTarget = finalTarget - 0.02;
   H_stepsize = 0.02; //0.02 * 0.6 * Global_Physical_Variables::lambda;
   currentH_Target = Global_Physical_Variables::H;
   while(Global_Physical_Variables::H > finalTarget*Global_Physical_Variables::lambda){
     std::cout << "Current height = " << Global_Physical_Variables::H 
	       << " Target height = " << finalTarget*Global_Physical_Variables::lambda
	       << " Target for next solution = " << currentH_Target
	       << std::endl;
     currentH_Target -= H_stepsize;
     problem2.lower_height(currentH_Target);
     problem2.save_solution();//to have something to restart from in futur
     problem2.save_solution(".");
     std::cout << "Current height = " << Global_Physical_Variables::H 
	       << " Target height = " << finalTarget*Global_Physical_Variables::lambda
	       << std::endl;
   }
   


}


void reload_solution_fresh_test(){

 // Initial values for parameter values
 Global_Physical_Variables::P=0.0; 
 Global_Physical_Variables::Gravity=0.0;

 // Create penetrator
 Global_Physical_Variables::Penetrator_pt =
  new AxiSymPenetrator(&Global_Physical_Variables::H);
  
if(false){
 Global_Physical_Variables::n_ele_theta = 2.0* Global_Physical_Variables::n_ele_r;
}
else{
 Global_Physical_Variables::n_ele_theta =(int) (1.571/(Global_Physical_Variables::t/   Global_Physical_Variables::n_ele_r)  + 0.5)*3;
}

 cout << "Thickness = " << Global_Physical_Variables::t << " with " <<
        Global_Physical_Variables::n_ele_r << ", " << 
        Global_Physical_Variables::n_ele_theta << 
        " elements in the r and theta direction" << endl;

 bool incompressible = true;

 cout << endl;

   // Set volume to correspond to the actual volume of the undefomred sphere
   Global_Physical_Variables::Volume = pow((1 - Global_Physical_Variables::t), 3) / 3;
 
  ofstream Solution_output_file;
  ifstream Solution_input_file;
  CantileverProblem<AxisymQPVDElementWithPressure> problem2(incompressible);
 

   problem2.set_under_relaxation_factor(Global_Physical_Variables::under_relaxation);  
   Global_Physical_Variables::Volume = problem2.calc_inflated_vol(Global_Physical_Variables::t, Global_Physical_Variables::lambda);


 //load previous soltution
  std::cout << "Loading solution: " << Global_Physical_Variables::loadSol << std::endl;
  Solution_input_file.open(Global_Physical_Variables::loadSol.c_str());
  problem2.read(Solution_input_file);
  Solution_input_file.close();

   //Lets check how ti looks
   problem2.doc_solution();

   std::cout << "Initial Newton Solve, shoudl converge immediatly. "
	     << "H = " << Global_Physical_Variables::H
	     << " Volume = " << Global_Physical_Variables::Volume
	     << " P = " << problem2.get_interal_pressure() << std::endl;

   //solve once in undeformed configuration, hopefully this should converge immediatly
   problem2.newton_solve();

   //Lets check that displacement is zero
   problem2.doc_solution();

   problem2.save_solution();
   problem2.save_solution(".");


}

//=======start_of_main==================================================
/// Driver for cantilever beam loaded by surface traction and/or
/// gravity
//======================================================================
int main(int argc, char **argv){
  typedef boost::iostreams::tee_device<ostream, ofstream> TeeDevice;
  typedef boost::iostreams::stream<TeeDevice> TeeStream;

  // Set the precision fo output to maximum (i.e. double) to avoid problems
  // due to rounding when reloading
  // cout.precision(17);  

  /// Redirect a copy of cout to file
  remove("spherical_solid_w_contact.log");
  ofstream logFile;
  logFile.open("spherical_solid_w_contact.log");

  ostream tmp(cout.rdbuf()); // <----
  TeeDevice outputDevice(tmp, logFile); // <----
  TeeStream logger(outputDevice);    

  cout.rdbuf(logger.rdbuf());
  cout << "spherical_solid_w_contact.cc log info" << endl;

   // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // Define possible command line arguments and parse the ones that
 // were actually specified
 
   // Number of elements in the r direction. Default is 4
 CommandLineArgs::specify_command_line_flag("--nreler",
                                            &Global_Physical_Variables::n_ele_r);
 
  // Thickness of shell. Default is 0.1
 CommandLineArgs::specify_command_line_flag("--tshell",
                                            &Global_Physical_Variables::t);
 
 // Number of integration points for new integration scheme. 
 // Use normal Gauss rule if not specified.
 CommandLineArgs::specify_command_line_flag("--nintpt",
                                            &Global_Physical_Variables::Nintpt);

  // Set max number of Newton Iterations. Default is 10
 CommandLineArgs::specify_command_line_flag("--maxnsteps",
                                            &Global_Physical_Variables::MaxNSteps);

  // Change bool to switch between applied traction and volume constraint
 CommandLineArgs::specify_command_line_flag("--enfvol",
                                            &Global_Physical_Variables::enforce_volume_constraint);
 
  // Set teh degree of preinflaction, with 1 being the undeformed state
 CommandLineArgs::specify_command_line_flag("--lambda",
                                            &Global_Physical_Variables::lambda);

  // Set the material properties of Mooney-Rivlin Constitutive Law
 CommandLineArgs::specify_command_line_flag("--C1",
                                            &Global_Physical_Variables::C1);
  // Set the material properties of Mooney-Rivlin Constitutive Law
 CommandLineArgs::specify_command_line_flag("--C2",
                                            &Global_Physical_Variables::C2);

 // Set the material properties of Mooney-Rivlin Constitutive Law
 CommandLineArgs::specify_command_line_flag("--load",
                                            &Global_Physical_Variables::loadSol);

 // Set height of penetrator
 CommandLineArgs::specify_command_line_flag("--H",
                                            &Global_Physical_Variables::H);

 // Set consitutive law
 CommandLineArgs::specify_command_line_flag("--claw",
                                            &Global_Physical_Variables::constitutive_law);

 CommandLineArgs::specify_command_line_flag("--program",
                                            &Global_Physical_Variables::program);

 //set the newton tolerance
 CommandLineArgs::specify_command_line_flag("--newtontol",
                                            &Global_Physical_Variables::newton_tol);

 //prefactor for contact residues to improve scaling
 CommandLineArgs::specify_command_line_flag("--contactscaling",
                                            &Global_Physical_Variables::contact_pressure_prefactor);

 CommandLineArgs::specify_command_line_flag("--underrelaxation",
                                            &Global_Physical_Variables::under_relaxation);


 // set the contact options
  ///Set whether or not to use isoparametric basis function for pressure
 CommandLineArgs::specify_command_line_flag("--contact_use_isoparametric",
                                            &Global_Physical_Variables::Contact_use_isoparametric_flag);

  ///Set options for basis/test functions for penetration and pressure
 CommandLineArgs::specify_command_line_flag("--contact_use_collocated_penetration",
                                            &Global_Physical_Variables::Contact_use_collocated_penetration_flag);

  ///Set options for basis/test functions for penetration and pressure
 CommandLineArgs::specify_command_line_flag("--contact_use_collocated_contact_pressure",
                                            &Global_Physical_Variables::Contact_use_collocated_contact_pressure_flag);

 CommandLineArgs::specify_command_line_flag("--debugInfo",
                                            &Global_Physical_Variables::printDebugInfo);



  // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
 
 
  // Doc use of integration scheme
 if (CommandLineArgs::command_line_flag_has_been_set("--nintpt"))
  {
   oomph_info << "Setting new integration scheme with nintpt="
              << Global_Physical_Variables::Nintpt << std::endl;
  }
 else
  {
   oomph_info << "Using Gauss scheme" << std::endl;
  }
 

  if (CommandLineArgs::command_line_flag_has_been_set("--nreler"))
  {
   oomph_info << "Setting number of radial elements ="
              << Global_Physical_Variables::n_ele_r << std::endl;
  }
 else
  {
   oomph_info << "Using default number (4) of radial elemens" << std::endl;
  }
  if (CommandLineArgs::command_line_flag_has_been_set("--lambda"))
  {
   oomph_info << "Preinflation (lambda)  ="
              << Global_Physical_Variables::lambda << std::endl;
  }

 if(Global_Physical_Variables::printDebugInfo){
 	cout << "Starting Main " << endl;
 }


 if(Global_Physical_Variables::constitutive_law == "GH"){
   //Create generalised Hookean constitutive equations
   Global_Physical_Variables::Constitutive_law_pt = 
     new GeneralisedHookean(&Global_Physical_Variables::Nu,
                        &Global_Physical_Variables::E);

   std::cout << "Setting constitutive law to GeneralisedHookean" << std::endl;
 }

 if(Global_Physical_Variables::constitutive_law == "MR"){
  // Create MooneyRivlin constitutive equations
  Global_Physical_Variables::Strain_energy_function_pt = 
  new MooneyRivlin(&Global_Physical_Variables::C1,
                        &Global_Physical_Variables::C2);

  // Define a constitutive law (based on strAxisymmetricSolidTractionElementain energy function)
  Global_Physical_Variables::Constitutive_law_pt = 
  new IsotropicStrainEnergyFunctionConstitutiveLaw(
   Global_Physical_Variables::Strain_energy_function_pt);

   std::cout << "Setting constitutive law to MooneyRivlin" << std::endl;
 }



 // switch which sub-routine to run 
 if(!Global_Physical_Variables::program.compare("standard")){
   std::cout << "Starting standard_run()" <<std::endl;
   standard();
 }
 else{
   if(!Global_Physical_Variables::program.compare("standardSmallSteps")){
     std::cout << "standardSmallSteps" <<std::endl;
     standard_run(0.01);
   }
   if(!Global_Physical_Variables::program.compare("lowerC1")){
     std::cout << "lower_C1()" <<std::endl;
     lower_C1();
   }
   if(!Global_Physical_Variables::program.compare("lowerHto0.8")){
     std::cout << "lowerHto0.8()" <<std::endl;
     change_parameter(Global_Physical_Variables::H, 0.8);
   }
   if(!Global_Physical_Variables::program.compare("increaseVol"))
     {
       std::cout << "Decreasing/increaseing Volume to lambda = 1.1."
		 << " (Current lambda = " 
		 << Global_Physical_Variables::lambda
		 << " )"  <<std::endl;
       //increase volume
       change_parameter(Global_Physical_Variables::lambda, 1.1);
     }
  if(!Global_Physical_Variables::program.compare("reloadSolution"))
     {
       std::cout << "Starting program reload_solution_fresh_test()" << std::endl;
       reload_solution_fresh_test();
     } 
 }
 //reload_solution_test();
 // reload_solution_fresh_test();

    
 logger.close();
 
return 0;
 
} //end of main
