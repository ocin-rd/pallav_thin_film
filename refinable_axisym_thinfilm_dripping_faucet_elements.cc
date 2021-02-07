//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
#include "refinable_axisym_thinfilm_dripping_faucet_elements.h"


namespace oomph
{

//========================================================================
/// Add element's contribution to the elemental 
/// residual vector and/or Jacobian matrix.
/// flag=1: compute both
/// flag=0: compute only residual vector
//========================================================================
void RefineableAxisymmetricThinFilmDrippingFaucetEquations::
fill_in_generic_residual_contribution_axisym_thinfilm_dripping_faucet
(Vector<double> &residuals, 
 DenseMatrix<double> &jacobian, 
 const unsigned& flag)
{
 // Return immediately if there are no dofs
 if (ndof()==0) return;

 //Find out how many nodes there are in the element
 const unsigned n_node = nnode();
 
 // Get continuous time from timestepper of first node
 double time=node_pt(0)->time_stepper_pt()->time_pt()->time();

 //Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,1), dtestdx(n_node,1);
 
 //Set the value of n_intpt
 unsigned n_intpt = integral_pt()->nweight();

 //The local index at which the height variable is stored
 unsigned h_nodal_index = this->h_index_axisym_thinfilm_dripping_faucet();
 unsigned u_nodal_index = this->u_index_axisym_thinfilm_dripping_faucet();
 unsigned omega_nodal_index = this->omega_index_axisym_thinfilm_dripping_faucet();

 // Get the Ohnesorg number
 const double ohnesorg = oh();

 // Curvature term
 double curvature = 0.0;

 //Integers to store the local equation and unknown numbers
 int local_eqn_h=0, local_eqn_u=0, local_eqn_omega=0, local_unknown=0;

 // Local storage for pointers to hang_info objects
 HangInfo *hang_info_pt=0, *hang_info2_pt=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
 
   //Call the derivatives of the shape and test functions
   double J = 
     this->dshape_and_dtest_eulerian_at_knot_axisym_thinfilm_dripping_faucet(ipt,psi,dpsidx,test,dtestdx);
 
   //Premultiply the weights and the Jacobian
   double W = w*J;

   //Calculate local values of unknown
   //Allocate and initialise to zero
   double interpolated_z=0.0;
   double interpolated_h=0.0;
   double interpolated_u=0.0;
   double interpolated_omega=0.0;
   double interpolated_dhdz=0.0;
   double interpolated_dudz=0.0;
   double interpolated_domegadz=0.0;
   double interpolated_dhdt=0.0;
   double interpolated_dudt=0.0;
   double interpolated_mesh_velocity=0.0;
   
   //Calculate function value and derivatives:
   //-----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Get the nodal value of the height unknown
     double h_value = raw_nodal_value(l,h_nodal_index);
     double u_value = raw_nodal_value(l,u_nodal_index);
     double omega_value = raw_nodal_value(l,omega_nodal_index);
     interpolated_z += raw_nodal_position(l,0)*psi(l);
     interpolated_h += h_value*psi(l);
     interpolated_u += u_value*psi(l);
     interpolated_omega += omega_value*psi(l);
     interpolated_dhdz += h_value*dpsidx(l,0);
     interpolated_dudz += u_value*dpsidx(l,0);
     interpolated_domegadz += omega_value*dpsidx(l,0);
     interpolated_dhdt += dh_dt_axisym_thinfilm_dripping_faucet(l)*psi(l);
     interpolated_dudt += du_dt_axisym_thinfilm_dripping_faucet(l)*psi(l);
     interpolated_mesh_velocity += this->dnodal_position_dt(l,0)*psi(l);
    }

   // Subtract mesh velocity from dh_dt and du_dt
   //interpolated_dudt -= interpolated_mesh_velocity*interpolated_dudz;
   //interpolated_dhdt -= interpolated_mesh_velocity*interpolated_omega;

   //Get body force function
   //-------------------
   double body_force;
   this->get_body_force_axisym_thinfilm_dripping_faucet(ipt,time,body_force);
 
   // Assemble residuals and Jacobian
 
   // Loop over the nodes for the test functions 
   for(unsigned l=0;l<n_node;l++)
    {
     //Local variables used to store the number of master nodes and the
     //weight associated with the shape function if the node is hanging
     unsigned n_master=1; double hang_weight=1.0;

     //Local bool (is the node hanging)
     bool is_node_hanging = this->node_pt(l)->is_hanging();

     //If the node is hanging, get the number of master nodes
     if(is_node_hanging)
      {
       hang_info_pt = this->node_pt(l)->hanging_pt();
       n_master = hang_info_pt->nmaster();
      }
    //Otherwise there is just one master node, the node itself
    else
     {
      n_master = 1;
     }
   
    //Loop over the master nodes
    for(unsigned m=0;m<n_master;m++)
    {
     //Get the local equation number and hang_weight
     //If the node is hanging
     if(is_node_hanging)
      {
       //Read out the local equation number from the m-th master node
       local_eqn_h = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
					 h_nodal_index);
       local_eqn_u = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
					 u_nodal_index);
       local_eqn_omega = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
					 omega_nodal_index);

       //Read out the weight from the master node
       hang_weight = hang_info_pt->master_weight(m);
      }
     //If the node is not hanging
     else
      {
       //The local equation number comes from the node itself
       local_eqn_h = this->nodal_local_eqn(l,h_nodal_index);
       local_eqn_u = this->nodal_local_eqn(l,u_nodal_index);
       local_eqn_omega = this->nodal_local_eqn(l,omega_nodal_index);

       //The hang weight is one
       hang_weight = 1.0;
      }
     
     //If the nodal equation is not a boundary condition
     if(local_eqn_h >= 0)
      {
       // Add contributions from the thin film model
       residuals[local_eqn_h] += (interpolated_dhdt + interpolated_u*interpolated_omega +
                                interpolated_h/2.0*interpolated_dudz)*test(l)*W*hang_weight;
       
       // Calculate the Jacobian
       if(flag)
        {
         //Local variables to store the number of master nodes
         //and the weights associated with each hanging node
         unsigned n_master2=1; double hang_weight2=1.0;

         //Loop over the nodes for the variables
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           //Local bool (is the node hanging)
           bool is_node2_hanging = this->node_pt(l2)->is_hanging();

           //If the node is hanging, get the number of master nodes
           if(is_node2_hanging)
            {
             hang_info2_pt = this->node_pt(l2)->hanging_pt();
             n_master2 = hang_info2_pt->nmaster();
            }
           //Otherwise there is one master node, the node itself
           else
            {
             n_master2 = 1;
            }
           
           //Loop over the master nodes
           for(unsigned m2=0;m2<n_master2;m2++)
            {
             //Get the local unknown and weight
             //If the node is hanging
             if(is_node2_hanging)
              {
               //Read out the local unknown from the master node
               local_unknown = 
                this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                     h_nodal_index);

               //Read out the hanging weight from the master node
               hang_weight2 = hang_info2_pt->master_weight(m2);
              }
             //If the node is not hanging
             else
              {
               //The local unknown number comes from the node
               local_unknown = this->nodal_local_eqn(l2,h_nodal_index);

               //The hang weight is one
               hang_weight2 = 1.0;
              }

             //If the unknown is not pinned
             if(local_unknown >= 0)
              {
               //Add contribution to Elemental Matrix
	             // obacht not properly implemented -- finite difference
	             jacobian(local_eqn_h,local_unknown) += 
		            dpsidx(l2,0)*dtestdx(l,0)*W*hang_weight*hang_weight2;
              }
            } //End of loop over master nodes
          } //End of loop over nodes
        } //End of Jacobian calculation
       
      } //End of case when residual equation is not pinned
     
     //If the nodal equation is not a boundary condition
     if(local_eqn_u >= 0)
      {
        // Calculate curvature term
        curvature = 1.0/(interpolated_h*sqrt(1.0 + interpolated_omega*interpolated_omega)) -
          interpolated_domegadz/pow(1.0 + interpolated_omega*interpolated_omega, 3.0/2.0);

       // Add contributions from the thin film model
       residuals[local_eqn_u] += (interpolated_h*interpolated_dudt + interpolated_h*interpolated_u*interpolated_dudz)*test(l)*W*hang_weight;
       residuals[local_eqn_u] += (-interpolated_h*curvature + interpolated_h*3.0*ohnesorg*interpolated_dudz)*dtestdx(l,0)*W*hang_weight;
       residuals[local_eqn_u] += (-6.0*ohnesorg*interpolated_omega*interpolated_dudz - curvature*interpolated_omega)*test(l)*W*hang_weight;
       residuals[local_eqn_u] += 3.0*ohnesorg*interpolated_omega*interpolated_dudz*test(l)*W*hang_weight;
       residuals[local_eqn_u] += -interpolated_h*body_force*test(l)*W*hang_weight;

/*       residuals[local_eqn_u] += (interpolated_dudt + interpolated_u*interpolated_dudz)*test(l)*W*hang_weight;
       residuals[local_eqn_u] += (-curvature + 3.0*ohnesorg*interpolated_dudz)*dtestdx(l,0)*W*hang_weight;
       residuals[local_eqn_u] += (-6.0*ohnesorg*interpolated_omega*interpolated_dudz/interpolated_h - body_force)*test(l)*W*hang_weight;
*/       
       // Calculate the Jacobian
       if(flag)
        {
         //Local variables to store the number of master nodes
         //and the weights associated with each hanging node
         unsigned n_master2=1; double hang_weight2=1.0;

         //Loop over the nodes for the variables
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           //Local bool (is the node hanging)
           bool is_node2_hanging = this->node_pt(l2)->is_hanging();

           //If the node is hanging, get the number of master nodes
           if(is_node2_hanging)
            {
             hang_info2_pt = this->node_pt(l2)->hanging_pt();
             n_master2 = hang_info2_pt->nmaster();
            }
           //Otherwise there is one master node, the node itself
           else
            {
             n_master2 = 1;
            }
           
           //Loop over the master nodes
           for(unsigned m2=0;m2<n_master2;m2++)
            {
             //Get the local unknown and weight
             //If the node is hanging
             if(is_node2_hanging)
              {
               //Read out the local unknown from the master node
               local_unknown = 
                this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                     h_nodal_index);

               //Read out the hanging weight from the master node
               hang_weight2 = hang_info2_pt->master_weight(m2);
              }
             //If the node is not hanging
             else
              {
               //The local unknown number comes from the node
               local_unknown = this->nodal_local_eqn(l2,u_nodal_index);

               //The hang weight is one
               hang_weight2 = 1.0;
              }

             //If the unknown is not pinned
             if(local_unknown >= 0)
              {
               //Add contribution to Elemental Matrix
	             // obacht not properly implemented -- finite difference
	             jacobian(local_eqn_h,local_unknown) += 
		            dpsidx(l2,0)*dtestdx(l,0)*W*hang_weight*hang_weight2;
              }
            } //End of loop over master nodes
          } //End of loop over nodes
        } //End of Jacobian calculation
       
      } //End of case when residual equation is not pinned
     // IF it's not a boundary condition
     if(local_eqn_omega >= 0)
      {
       residuals[local_eqn_omega] += (interpolated_dhdz - interpolated_omega)*test(l)*W*hang_weight;
       
       // Calculate the Jacobian
       if(flag)
        {
         //Local variables to store the number of master nodes
         //and the weights associated with each hanging node
         unsigned n_master2=1; double hang_weight2=1.0;

         //Loop over the nodes for the variables
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           //Local bool (is the node hanging)
           bool is_node2_hanging = this->node_pt(l2)->is_hanging();

           //If the node is hanging, get the number of master nodes
           if(is_node2_hanging)
            {
             hang_info2_pt = this->node_pt(l2)->hanging_pt();
             n_master2 = hang_info2_pt->nmaster();
            }
           //Otherwise there is one master node, the node itself
           else
            {
             n_master2 = 1;
            }
           
           //Loop over the master nodes
           for(unsigned m2=0;m2<n_master2;m2++)
            {
             //Get the local unknown and weight
             //If the node is hanging
             if(is_node2_hanging)
              {
               //Read out the local unknown from the master node
               local_unknown = 
                this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                     h_nodal_index);

               //Read out the hanging weight from the master node
               hang_weight2 = hang_info2_pt->master_weight(m2);
              }
             //If the node is not hanging
             else
              {
               //The local unknown number comes from the node
               local_unknown = this->nodal_local_eqn(l2,u_nodal_index);

               //The hang weight is one
               hang_weight2 = 1.0;
              }

             //If the unknown is not pinned
             if(local_unknown >= 0)
              {
               //Add contribution to Elemental Matrix
	             // obacht not properly implemented -- finite difference
	             jacobian(local_eqn_h,local_unknown) += 
		            dpsidx(l2,0)*dtestdx(l,0)*W*hang_weight*hang_weight2;
              }
            } //End of loop over master nodes
          } //End of loop over nodes
        } //End of Jacobian calculation
      } 
    } //End of loop over master nodes for residual
  } //End of loop over nodes
 
} // End of loop over integration points
}


//====================================================================
// Force build of templates
//====================================================================
template class RefineableAxisymmetricThinFilmDrippingFaucetElement<2>;
template class RefineableAxisymmetricThinFilmDrippingFaucetElement<3>;
template class RefineableAxisymmetricThinFilmDrippingFaucetElement<4>;

}
