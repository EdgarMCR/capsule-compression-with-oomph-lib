//=======================================================================
/// Return the residuals for the AxisymmetricSolidTractionElements
//=======================================================================
template<class ELEMENT>
void AxisymmetricSolidTractionElement<ELEMENT>::
fill_in_contribution_to_residuals(Vector<double> &residuals)
{
 //Find out how many nodes there are
 unsigned n_node = nnode();

 //Find out how many positional dofs there are
 //unsigned n_position_type = nnodal_position_type(); //Old version
 unsigned n_position_type = this->nnodal_position_type();

 //Integer to hold the local equation number
 int local_eqn=0;

 //Set up memory for the shape functions
 //The surface is 1D, so we only have one local derivative
 Shape psi(n_node,n_position_type);
 DShape dpsids(n_node,n_position_type,1); 

 //Set the value of n_intpt
 unsigned n_intpt = integral_pt()->nweight();

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   //Only need to call the local derivatives
   dshape_local_at_knot(ipt,psi,dpsids);

   //Calculate the global position and lagrangian coordinate
   Vector<double> interpolated_x(2,0.0); 
   Vector<double> interpolated_xi(2,0.0);

   //Calculate the global and lagrangian derivtives wrt the local coordinates
   Vector<double> interpolated_dxds(2,0.0); 
   Vector<double> interpolated_dxids(2,0.0);
 
   //NEW
   //Also calculate the surface Vectors (derivatives wrt local coordinates)
   DenseMatrix<double> interpolated_A(1,2,0.0); 

   //Calculate displacements and derivatives
   for(unsigned l=0;l<n_node;l++) 
    {
     //Loop over positional dofs
     for(unsigned k=0;k<n_position_type;k++)
      {
       //Loop over the number of lagrangian coordinates (2)
       for(unsigned i=0;i<2;i++)
        {
         //Calculate the global position
         interpolated_x[i] += 
          nodal_position_gen(l,bulk_position_type(k),i)*psi(l,k);

         interpolated_xi[i] += 
          this->lagrangian_position_gen(l,bulk_position_type(k),i)*psi(l,k);

         //Calculate the derivatives of the global and lagrangian coordinates
         interpolated_dxds[i] += 
          nodal_position_gen(l,bulk_position_type(k),i)*dpsids(l,k,0);

         interpolated_dxids[i] += 
          this->lagrangian_position_gen(l,bulk_position_type(k),i)
          *dpsids(l,k,0);

          //Loop over LOCAL derivative directions, to calculate the tangent(s)
         for(unsigned j=0;j<1;j++)
          {
           interpolated_A(j,i) += 
            nodal_position_gen(l,bulk_position_type(k),i)*dpsids(l,k,j);
          }
        }
      }
    }

   //Now calculate the entries of the deformed surface metric tensor
   //Now find the local deformed metric tensor from the tangent Vectors
   DenseMatrix<double> A(2);
   //The off-diagonal terms are Zero 
   A(0,1) = A(1,0) = 0.0;
   //The diagonal terms are a little complicated
   A(0,0) =  
    (interpolated_dxds[0] - interpolated_x[1]*interpolated_dxids[1])*
    (interpolated_dxds[0] - interpolated_x[1]*interpolated_dxids[1]) +
    (interpolated_dxds[1] + interpolated_x[0]*interpolated_dxids[1])*
    (interpolated_dxds[1] + interpolated_x[0]*interpolated_dxids[1]);


   A(1,1) =  (interpolated_x[0]*sin(interpolated_xi[1]) +
               interpolated_x[1]*cos(interpolated_xi[1]))*
    (interpolated_x[0]*sin(interpolated_xi[1]) +
     interpolated_x[1]*cos(interpolated_xi[1]));

   //Premultiply the weights and the square-root of the determinant of 
   //the metric tensor
   double W = w*sqrt(A(0,0)*A(1,1));
   
   //Get the outer unit normal
   Vector<double> interpolated_normal(2);
   
   /*if(b)
   {
        std::cout << "ipt = " << ipt << std::endl;
   }*/
   //New method of finding outer unit normal is to call function
   outer_unit_normal(ipt,interpolated_normal);
   
   
   //std::cout << "Modern normel: " << interpolated_normal[0] << ", " << interpolated_normal[1] << std::endl;
   
   //Old way

   //Also find the normal -- just the cross product of the metric tensors
   //but I want to express it in terms of e_r and e_theta components
   //N.B. There is an issue at theta = 0,pi, where the normal is e_{r},
   //but given that I never assemble it, should be OK!
   //The minus sign is chosen to ensure that the normal is really outward 

   //Component in the e_{r} direction
   interpolated_normal[0] = -1.0*
    (interpolated_x[0]*sin(interpolated_xi[1]) +
     interpolated_x[1]*cos(interpolated_xi[1]))*
    (interpolated_dxds[1] + interpolated_x[0]*interpolated_dxids[1]);
   //Component in the e_{theta} direction
   interpolated_normal[1] =  -1.0*
    (interpolated_x[0]*sin(interpolated_xi[1]) +
     interpolated_x[1]*cos(interpolated_xi[1]))*
    (interpolated_x[1]*interpolated_dxids[1] - interpolated_dxds[0]);
   
   //TODO: Fix normal direction!
   //Huge assumption: we are not going to be on north or south face
   //If we're on the north or south face need to flip normal

   //Now adjust and scale the normal
   double length = 0.0;
   for(unsigned i=0;i<2;i++)
    {
     interpolated_normal[i] *= normal_sign();
     length += interpolated_normal[i]*interpolated_normal[i];
    }
   for(unsigned i=0;i<2;i++)
    {
     interpolated_normal[i] /= sqrt(length);
    }

   //Now calculate the load
   Vector<double> traction(2);

   get_traction(ipt, interpolated_xi, interpolated_x, interpolated_normal,
                traction);

   //Normal is outwards
   //get_traction(time(),interpolated_x,interpolated_normal,traction); //Original


   //=====LOAD TERMS  FROM PRINCIPLE OF VIRTUAL DISPLACEMENTS========
   
   // Set volume to zero
   volume = 0.0;
   Vector<double> cart_pos(2);
//   Vector<double> cart_normal(2);
   
   double actual_angle = interpolated_x[1] + interpolated_xi[1] ;

   //Loop over the test functions, nodes of the element
   for(unsigned l=0;l<n_node;l++)
    {
     //Loop of types of dofs
     for(unsigned k=0;k<n_position_type;k++)
      {
         // added by Edgar 31/5/2016
        // dividing by the number of dimensions (2) and taking the vector produce of 
        // the position with the normal. Multiplying the normal by -1 as the elements are on the
        // inner face and point towards the origin while I want the normal going awat from the origin
        // dS = r dr dtheta
        volume += 0.5 * interpolated_x[0]*interpolated_normal[0]*psi(l,k)*W*interpolated_x[0]*-1;
        volume += 0.5 * actual_angle*interpolated_normal[1]*psi(l,k)*W*-1;
        
       //Loop over the displacement components
       for(unsigned i=0;i<2;i++)
        {
         local_eqn = 
          this->position_local_eqn(l,bulk_position_type(k),i);
         /*IF it's not a boundary condition*/
         if(local_eqn >= 0)
          {
           //Add the loading terms to the residuals
           residuals[local_eqn] -= traction[i]*psi(l,k)*W;
          }
          
        }
      } //End of if not boundary condition
    } //End of loop over shape functions
  } //End of loop over integration points
}