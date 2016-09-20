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

#include "axisym_solid_traction_elements.h"

//Need to place the function in a .cc file as otherwise the compiler complains about an undefined function

//=======================================================================
/// Return the residuals for the AxisSymSolidVolumeConstraintElement
//=======================================================================
template <class ELEMENT>
void AxisSymSolidVolumeConstraintElement<ELEMENT>::
fill_in_contribution_to_residuals(Vector<double> &residuals)
{
  // Note: This element can only be used with the associated 
  // AxisymmetricSolidTractionVolumeConstraintElement elements 
  // which compute the actual enclosed volume; here we only add 
  // the contribution to the residual
  
  //volume under the mesh
  double vol =0.0;
  
  //Get volume from face elements
  unsigned n_element=Traction_vol_const_mesh_pt->nelement();
  for(unsigned i=0;i<n_element;i++)
    {
      //Cast to a solid traction element
      AxisymmetricSolidTractionVolumeConstraintElement<ELEMENT> *el_pt = 
        dynamic_cast<AxisymmetricSolidTractionVolumeConstraintElement<ELEMENT>*>(Traction_vol_const_mesh_pt->element_pt(i));
          
      vol += el_pt->get_volume();
    }
  
//   const int local_eqn = this->ptraded_local_eqn();
//   if(local_eqn >= 0)
//    {
  residuals[0] -= vol - *Prescribed_volume_pt;
//    }
}

  // Function to set pointer to surface element mesh and external data
  // This cannot be done in the constructor as the traction elements
  // have not been created yet.
  template <class ELEMENT>
  void AxisSymSolidVolumeConstraintElement<ELEMENT>::
  void set_mesh_pt_and_external(SolidMesh* traction_vol_const_mesh_pt){
    
    ///Need to loop over all traction elements and make them external data
    if(added_external){
      std::cout << "Already added external data! Function AxisSymSolidVolumeConstraintElement::set_mesh_pt_and_external will return now." << std::endl;
      return;
    }
    else{
      added_external=true;
    }
    
    Traction_vol_const_mesh_pt = traction_vol_const_mesh_pt;
    

    
  //The nodes must be external data
  unsigned n_element =traction_vol_const_mesh_pt->nelement();
  for(unsigned i=0;i<n_element;i++)
    {
      //Cast to a solid traction element
      AxisymmetricSolidTractionVolumeConstraintElement<ELEMENT> *el_pt = 
        dynamic_cast<AxisymmetricSolidTractionVolumeConstraintElement<ELEMENT>*>
          (traction_vol_const_mesh_pt->element_pt(i));
    
      unsigned nnod_el=el_pt->nnode();
      // Loop over all nodes in element
      for (unsigned j=0;j<nnod_el;j++)
      {
        //add each node as external data
          add_external_data(el_pt->node_pt(j),true);
      }
    }
  }
