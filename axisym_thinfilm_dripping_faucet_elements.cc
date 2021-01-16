//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1307 $
//LIC//
//LIC// $LastChangedDate: 2018-01-18 11:30:14 +0000 (Thu, 18 Jan 2018) $
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
//Non-inline functions for axisymmetric thin film DrippingFaucet elements
#include "axisym_thinfilm_dripping_faucet_elements.h"


namespace oomph
{


//======================================================================
/// Set the data for the number of Variables at each node, always three
/// in every case (h, u, \omega)
//======================================================================
 template<unsigned NNODE_1D>
 const unsigned AxisymmetricThinFilmDrippingFaucetElement<NNODE_1D>::Initial_Nvalue = 2;

  double AxisymmetricThinFilmDrippingFaucetEquations::
  Default_Physical_Constant_Value = 0.0;
  double AxisymmetricThinFilmDrippingFaucetEquations::
  Default_Physical_Ratio_Value = 1.0;

//======================================================================
/// Compute element residual Vector and/or element Jacobian matrix 
/// 
/// flag=1: compute both
/// flag=0: compute only residual Vector
///
/// Pure version without hanging nodes
//======================================================================
void AxisymmetricThinFilmDrippingFaucetEquations::
fill_in_generic_residual_contribution_axisym_thinfilm_dripping_faucet
(Vector<double> &residuals, 
 DenseMatrix<double> &jacobian, 
 const unsigned& flag) 
{
 // Return immediately if there are no dofs
 if (ndof()==0) return;

 //Find out how many nodes there are
 const unsigned n_node = nnode();
 
 // Get continuous time from timestepper of first node
 double time=node_pt(0)->time_stepper_pt()->time_pt()->time();

 //Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,1), dtestdx(n_node,1);

 //Index at which the height unknown is stored
 const unsigned h_nodal_index = h_index_axisym_thinfilm_dripping_faucet();
 const unsigned u_nodal_index = u_index_axisym_thinfilm_dripping_faucet();
 const unsigned omega_nodal_index = omega_index_axisym_thinfilm_dripping_faucet();
 
 //Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();

 // Get the Ohnesorg number
 const double ohnesorg = oh();

 //Integers to store the local equation and unknown numbers
 int local_eqn=0, local_unknown=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   //Call the derivatives of the shape and test functions
   double J = dshape_and_dtest_eulerian_at_knot_axisym_thinfilm_dripping_faucet
     (ipt,psi,dpsidx,test,dtestdx);
       
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
    }


   //Get body force function
   //-------------------
   double body_force;
   get_body_force_axisym_thinfilm_dripping_faucet(ipt,time,body_force);

   // Assemble residuals and Jacobian
   //--------------------------------
       
   // Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {
     //Get the local equation
     local_eqn = nodal_local_eqn(l,h_nodal_index);
     // IF it's not a boundary condition
     if(local_eqn >= 0)
      {
       // Add contributions from the thin film model
       residuals[local_eqn] += (interpolated_dhdt + interpolated_u*interpolated_omega +
                                interpolated_h/2.0*interpolated_dudz)*test(l)*W;

       // Calculate the jacobian (currently not implemented)
       //-----------------------
       if(flag)
        {
         //Loop over the height shape functions again
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           local_unknown = nodal_local_eqn(l2,h_nodal_index);
           //If at a non-zero degree of freedom add in the entry
           if(local_unknown >= 0)
            {
             //Add contribution to Elemental Matrix
	           jacobian(local_eqn, local_unknown) += 0.0*psi(l2)*test(l)*W;
            }
          }
        }
      }

     //Get the local equation
     local_eqn = nodal_local_eqn(l,u_nodal_index);
     // IF it's not a boundary condition
     if(local_eqn >= 0)
      {
       // Add contributions from the thin film model
       residuals[local_eqn] += (interpolated_dudt + interpolated_u*interpolated_dudz)*test(l)*W;

       // Calculate the jacobian (currently not implemented)
       //-----------------------
       if(flag)
        {
         //Loop over the height shape functions again
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           local_unknown = nodal_local_eqn(l2,h_nodal_index);
           //If at a non-zero degree of freedom add in the entry
           if(local_unknown >= 0)
            {
             //Add contribution to Elemental Matrix
	           jacobian(local_eqn, local_unknown) += 0.0*psi(l2)*test(l)*W;
            }
          }
        }
      }

      //Get the local equation
      local_eqn = nodal_local_eqn(l,omega_nodal_index);
      // IF it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] += (interpolated_dhdz - interpolated_omega)*test(l)*W;
       }
    }

  } // End of loop over integration points
}   

//======================================================================
/// Self-test:  Return 0 for OK
//======================================================================
unsigned  AxisymmetricThinFilmDrippingFaucetEquations::self_test()
{

 bool passed=true;

 // Check lower-level stuff
 if (FiniteElement::self_test()!=0)
  {
   passed=false;
  }

 // Return verdict
 if (passed)
  {
   return 0;
  }
 else
  {
   return 1;
  }
   
}

//======================================================================
/// Output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
void AxisymmetricThinFilmDrippingFaucetEquations::output(std::ostream &outfile, 
						  const unsigned &nplot)
{

 //Vector of local coordinates
 Vector<double> s(1);
 
 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   outfile << interpolated_x(s,0) << " ";
    
   outfile << interpolated_h_axisym_thinfilm_dripping_faucet(s) << " ";

   outfile << interpolated_u_axisym_thinfilm_dripping_faucet(s) << " ";

   outfile << interpolated_omega_axisym_thinfilm_dripping_faucet(s) << " ";
   
   outfile << interpolated_dhdt_axisym_thinfilm_dripping_faucet(s) << std::endl;   
   
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(outfile,nplot);

}


//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
void AxisymmetricThinFilmDrippingFaucetEquations::output(FILE* file_pt,
						  const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(1);
 
 // Tecplot header info
 fprintf(file_pt,"%s",tecplot_zone_string(nplot).c_str());

 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   fprintf(file_pt,"%g ",interpolated_x(s,0));
    
   fprintf(file_pt,"%g ",interpolated_h_axisym_thinfilm_dripping_faucet(s));
    
   fprintf(file_pt,"%g ",interpolated_u_axisym_thinfilm_dripping_faucet(s));
    
   fprintf(file_pt,"%g ",interpolated_omega_axisym_thinfilm_dripping_faucet(s));
    
   fprintf(file_pt,"%g \n",interpolated_dhdt_axisym_thinfilm_dripping_faucet(s));
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(file_pt,nplot);
}

//=======================================================================
/// Compute norm of the solution
//=======================================================================
void AxisymmetricThinFilmDrippingFaucetEquations::compute_norm(double& norm)
{
 
 // Initialise
 norm=0.0;
 
 //Vector of local coordinates
 Vector<double> s(1);
 
 // Solution
 double h=0.0;
 
 //Find out how many nodes there are in the element
 unsigned n_node = this->nnode();
 
 Shape psi(n_node);
 
 //Set the value of n_intpt
 unsigned n_intpt = this->integral_pt()->nweight();
 
 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   
    //Assign values of s
    s[0] = this->integral_pt()->knot(ipt,0);
   
   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);
   
   // Get jacobian of mapping
   double J=this->J_eulerian(s);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   // Get FE function value
   h=this->interpolated_h_axisym_thinfilm_dripping_faucet(s);
   
   // Add to  norm
   norm+=h*h*W;
  }
}

//====================================================================
// Force build of templates
//====================================================================
template class AxisymmetricThinFilmDrippingFaucetElement<2>;
template class AxisymmetricThinFilmDrippingFaucetElement<3>;
template class AxisymmetricThinFilmDrippingFaucetElement<4>;

}
