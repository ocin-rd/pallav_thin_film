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
//Header file for refineable thin film DrippingFaucet elements

#ifndef OOMPH_REFINEABLE_AXISYM_THINFILM_DRIPPING_FAUCET_ELEMENTS_HEADER
#define OOMPH_REFINEABLE_AXISYM_THINFILM_DRIPPING_FAUCET_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif


//oomph-lib headers
#include "../../../src/generic/refineable_quad_element.h"
#include "../../../src/generic/error_estimator.h"
#include "axisym_thinfilm_dripping_faucet_elements.cc"

namespace oomph
{

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//======================================================================
/// Refineable version of axisymmetric thin film DrippingFaucet equations
///
///
//======================================================================
class RefineableAxisymmetricThinFilmDrippingFaucetEquations : 
  public virtual AxisymmetricThinFilmDrippingFaucetEquations,
  public virtual RefineableElement,
  public virtual ElementWithZ2ErrorEstimator
{
  public:

 /// \short Constructor, simply call other constructors
 RefineableAxisymmetricThinFilmDrippingFaucetEquations() : AxisymmetricThinFilmDrippingFaucetEquations(),
  RefineableElement(), ElementWithZ2ErrorEstimator() 
  { } 

 /// Broken copy constructor
 RefineableAxisymmetricThinFilmDrippingFaucetEquations
   (const RefineableAxisymmetricThinFilmDrippingFaucetEquations& dummy) 
  { 
   BrokenCopy::broken_copy("RefineableAxisymmetricThinFilmDrippingFaucetEquations");
  } 
 
 /// Broken assignment operator
 void operator=(const RefineableAxisymmetricThinFilmDrippingFaucetEquations&) 
  {
   BrokenCopy::broken_assign("RefineableAxisymmetricThinFilmDrippingFaucetEquations");
  }
 
 /// Number of 'flux' terms for Z2 error estimation 
 unsigned num_Z2_flux_terms() {return 1;}

 /// Get 'flux' for Z2 error recovery:  Standard flux from axisym thin film equations
 void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
  {this->get_flux(s,flux);}
 
/// \short Get the function value h in Vector.
/// Note: Given the generality of the interface (this function
/// is usually called from black-box documentation or interpolation routines),
/// the values Vector sets its own size in here.
void get_interpolated_values(const Vector<double>&s,  Vector<double>& values)
 {
  // Set size of Vector: h, u, omega
  values.resize(3);
  
  //Find number of nodes
  unsigned n_node = nnode();
  
  //Local shape function
  Shape psi(n_node);
  
  //Find values of shape function
  shape(s,psi);
  
  //Initialise value of h, u, omega
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;

  //Find the index at which the height unknown is stored
  unsigned h_nodal_index = this->h_index_axisym_thinfilm_dripping_faucet();
  unsigned u_nodal_index = this->u_index_axisym_thinfilm_dripping_faucet();
  unsigned omega_nodal_index = this->omega_index_axisym_thinfilm_dripping_faucet();
  
  //Loop over the local nodes and sum up the values
  for(unsigned l=0;l<n_node;l++)
   {
    values[0] += this->nodal_value(l,h_nodal_index)*psi[l];
    values[1] += this->nodal_value(l,u_nodal_index)*psi[l];
    values[2] += this->nodal_value(l,omega_nodal_index)*psi[l];
   }
 }


 /// \short Get the function value u in Vector.
 /// Note: Given the generality of the interface (this function
 /// is usually called from black-box documentation or interpolation routines),
 /// the values Vector sets its own size in here.
 void get_interpolated_values(const unsigned& t, const Vector<double>&s, 
                              Vector<double>& values)
  {
   // Get the element's timestepper
   TimeStepper* time_stepper_pt = this->node_pt(0)->time_stepper_pt();

   if (!time_stepper_pt->is_steady())
    {
     values.resize(3);

     //Find number of nodes
     unsigned n_node = nnode();
  
     //Local shape function
     Shape psi(n_node);
  
     //Find values of shape function
     shape(s,psi);
  
     //Initialise value of h
     values[0] = 0.0;
     values[1] = 0.0;
     values[2] = 0.0;

     //Find the index at which the height unknown is stored
     const unsigned h_nodal_index = this->h_index_axisym_thinfilm_dripping_faucet();
     const unsigned u_nodal_index = this->u_index_axisym_thinfilm_dripping_faucet();
     const unsigned omega_nodal_index = this->omega_index_axisym_thinfilm_dripping_faucet();
  
     //Loop over the local nodes and sum up the values
     for(unsigned l=0;l<n_node;l++)
      {
       values[0] += nodal_value(t,l,h_nodal_index)*psi[l];
       values[1] += nodal_value(t,l,u_nodal_index)*psi[l];
       values[2] += nodal_value(t,l,omega_nodal_index)*psi[l];
      }
    }
   else
    {
     //Make sure that we call this particular object's steady 
     //get_interpolated_values (it could get overloaded lower down)
     RefineableAxisymmetricThinFilmDrippingFaucetEquations::get_interpolated_values(s,values);
    }
  }

 
 ///  Further build: Copy body force function pointer and Ohnesorg number from father element
 void further_build()
  {
   this->Body_force_fct_pt=dynamic_cast<RefineableAxisymmetricThinFilmDrippingFaucetEquations*>(
    this->father_element_pt())->body_force_fct_pt();
   this->Oh_pt=dynamic_cast<RefineableAxisymmetricThinFilmDrippingFaucetEquations*>(
    this->father_element_pt())->oh_pt();
  }


  private:


/// \short Add element's contribution to elemental residual vector and/or 
/// Jacobian matrix 
/// flag=1: compute both
/// flag=0: compute only residual vector
 void fill_in_generic_residual_contribution_axisym_thinfilm_dripping_faucet(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  const unsigned& flag); 
};


//======================================================================
/// Refineable version of 1D AxisymmetricThinFilmDrippingFaucetElement elements
///
///
//======================================================================
template <unsigned NNODE_1D>
 class RefineableAxisymmetricThinFilmDrippingFaucetElement : 
 public AxisymmetricThinFilmDrippingFaucetElement<NNODE_1D>,
 public virtual RefineableAxisymmetricThinFilmDrippingFaucetEquations,
 public virtual RefineableQElement<1>
{
  public:

 /// \short Constructor, simply call the other constructors 
 RefineableAxisymmetricThinFilmDrippingFaucetElement() : 
  RefineableElement(),
  RefineableAxisymmetricThinFilmDrippingFaucetEquations(),
  RefineableQElement<1>(),
  AxisymmetricThinFilmDrippingFaucetElement<NNODE_1D>()
   {} 


 /// Broken copy constructor
 RefineableAxisymmetricThinFilmDrippingFaucetElement
   (const RefineableAxisymmetricThinFilmDrippingFaucetElement<NNODE_1D>& dummy) 
  { 
   BrokenCopy::broken_copy("RefineableAxisymmetricThinFilmDrippingFaucetElement");
  } 
 
 /// Broken assignment operator
 void operator=(const RefineableAxisymmetricThinFilmDrippingFaucetElement<NNODE_1D>&) 
  {
   BrokenCopy::broken_assign("RefineableAxisymmetricThinFilmDrippingFaucetElement");
  }
 
 /// Number of continuously interpolated values: 2
 unsigned ncont_interpolated_values() const {return 2;}

 /// \short Number of vertex nodes in the element
 unsigned nvertex_node() const
  {return AxisymmetricThinFilmDrippingFaucetElement<NNODE_1D>::nvertex_node();}

 /// \short Pointer to the j-th vertex node in the element
 Node* vertex_node_pt(const unsigned& j) const
  {return AxisymmetricThinFilmDrippingFaucetElement<NNODE_1D>::vertex_node_pt(j);}

 /// Rebuild from sons: empty
 void rebuild_from_sons(Mesh* &mesh_pt) {}

 /// \short Order of recovery shape functions for Z2 error estimation:
 /// Same order as shape functions.
 unsigned nrecovery_order() {return (NNODE_1D-1);}

 ///  \short Perform additional hanging node procedures for variables
 /// that are not interpolated by all nodes. Empty.
 void further_setup_hanging_nodes(){}

};



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//=======================================================================
/// Face geometry for the 1D AxisymmetricThinFilmElement elements: Point elements
//=======================================================================
template<unsigned NNODE_1D>
class FaceGeometry<RefineableAxisymmetricThinFilmDrippingFaucetElement<NNODE_1D> >: 
 public virtual PointElement
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional QElement
 FaceGeometry() : PointElement() {}

};

}

#endif

