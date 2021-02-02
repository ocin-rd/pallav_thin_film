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
//Driver for 1D axisymmetric thin film DrippingFaucet problem

// Generic oomph-lib routines
#include "generic.h"

// Include axisymmetric thin film DrippingFaucet elements/equations
//#include "axisym_thinfilm_dripping_faucet_elements.h"
//#include "axisym_thinfilm_dripping_faucet_elements.cc"
#include "refinable_axisym_thinfilm_dripping_faucet_elements.cc"

// Include the mesh
#include "meshes/one_d_mesh.h"

using namespace std;

using namespace oomph;

namespace oomph
{

//=================================================================
/// Custom element that allows the handling of external data
//=================================================================
template <unsigned NNODE_1D>
class MyRefinableElement : public 
  RefineableAxisymmetricThinFilmDrippingFaucetElement<NNODE_1D>
{

public:
 
 ///\short Compute the element's residual vector and the Jacobian matrix.
 /// Jacobian is computed by finite-differencing.
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
  {
   //Add the contribution to the residuals
   this->fill_in_contribution_to_residuals(residuals);
   
   //Allocate storage for the full residuals (residuals of entire element)
   unsigned n_dof = this->ndof();
   Vector<double> full_residuals(n_dof);

   //Get the residuals for the entire element
   this->get_residuals(full_residuals);

   //There could be internal data
   //(finite-difference the lot by default)
   this->fill_in_jacobian_from_internal_by_fd(full_residuals,jacobian,true);

   //There could also be external data
   //(finite-difference the lot by default)
   this->fill_in_jacobian_from_external_by_fd(full_residuals,jacobian,true);

   //There could also be nodal data
   //(finite-difference the lot by default)
   this->fill_in_jacobian_from_nodal_by_fd(full_residuals,jacobian);
  }

 inline void update_in_external_fd(const unsigned &i)
  {
   //oomph_info << "Aloha!"<<std::endl;
   if(this->external_data_pt(0) != 0 &&
    this->external_data_pt(1) != 0)
    {
     // Get the current tip position
     double z_tip = this->external_data_pt(0)->value(0);
     // Get the current tip velocity
     double u_tip = this->external_data_pt(1)->value(0);
     // Loop over all nodes
     unsigned n_node = this->nnode();
     for(unsigned inod=0; inod<n_node; inod++)
      {
       Node* nod_pt = this->node_pt(inod);
       // Only update the tip node where the velocity is pinned
       if(nod_pt->is_pinned(1))
        {
         //oomph_info << "Updating position!"<<std::endl;
         nod_pt->x(0) = z_tip;
         //oomph_info << "Updating velocity!"<<std::endl;
         nod_pt->set_value(1,u_tip);
        }
      }
    }
  }

 inline void update_before_external_fd()
  {
   const unsigned i=0;
   update_in_external_fd(i);
  }

};

//====================================================================
/// Element to constrain the tip position and its velocity.
//====================================================================
template<class ELEMENT>
class DrippingFaucetConstraintElement : public GeneralisedElement
{
public:

  DrippingFaucetConstraintElement(const double& z_tip, const double& u_tip,
    const Mesh* mesh_pt, Node* tip_node_pt, const double* prescribed_drop_volume_pt)
  {
    // Create an internal data object, which is the unknown
    // position of the drop tip
    this->add_internal_data(new Data(1),true);
    this->internal_data_pt(0)->set_value(0,z_tip);
    // Create an internal data object, which is the unknown
    // velocity of the drop tip
    this->add_internal_data(new Data(1),true);
    this->internal_data_pt(1)->set_value(0,u_tip);

    // Set pointer to the mesh
    Mesh_pt = mesh_pt;
    // Set pointer to the tip node
    Tip_node_pt = tip_node_pt;
    // Set pointer to prescribed drop volume
    Prescribed_drop_volume_pt = prescribed_drop_volume_pt;

    // Add external data
    const unsigned n_node = Mesh_pt->nnode();
    for (unsigned inod=0; inod < n_node; inod++)
    {
      this->add_external_data(Mesh_pt->node_pt(inod), true);
    }
  }

  ~DrippingFaucetConstraintElement(){}

  double calculate_drop_volume()
  {
    // get the volume from the elements
    double volume = 0.0;
    unsigned n_element=Mesh_pt->nelement();
    for(unsigned el=0; el<n_element; el++)
    {
      // Get pointer to element
      ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Mesh_pt->element_pt(el));
    
      // Add its contribution
      volume += el_pt->compute_physical_size();
    }
    return volume;
  }

  void fill_in_contribution_to_residuals(Vector<double>& residuals)
  {
    // Get the current velocity of the tip
    const double current_u_tip = Tip_node_pt->dposition_dt(0);

    // Get the current volume of the drop
    const double current_drop_volume = calculate_drop_volume();

    // Equation for tip position
    const int z_tip_eqn=internal_local_eqn(0,0);

    // Equation for tip velocity
    const int u_tip_eqn=internal_local_eqn(1,0);

    if(z_tip_eqn >= 0)
    {
      // Get the residuals
      residuals[z_tip_eqn] += current_drop_volume - *Prescribed_drop_volume_pt;
    }

    if(u_tip_eqn >= 0)
    {
      // Get the residuals
      residuals[u_tip_eqn] += current_u_tip - this->internal_data_pt(1)->value(0);
    }
  }

 inline void update_in_internal_fd(const unsigned &i)
  {
    // Get the current tip position
    double z_tip = this->internal_data_pt(0)->value(0);
    // Get the current tip velocity
    double u_tip = this->internal_data_pt(1)->value(0);
    //oomph_info << "Updating position!"<<std::endl;
    Tip_node_pt->x(0) = z_tip;
    //oomph_info << "Updating velocity!"<<std::endl;
    Tip_node_pt->set_value(1,u_tip);
  }

 inline void update_before_internal_fd()
  {
   const unsigned i=0;
   update_in_internal_fd(i);
  }

private:

  const Mesh* Mesh_pt;
  Node* Tip_node_pt;
  const double* Prescribed_drop_volume_pt;
};

//==start_of_namespace================================================
/// Namespace for problem parameters
//====================================================================
namespace Problem_Parameter
{
 /// Output directory
 std::string Directory = "RESLT";

 DocInfo Doc_info;

 /// Faucet Radius
 double R = 1.0;

 /// Initial tip position
 double Z_tip = 1.0;
 Data* Tip_position_pt = 0;

 /// Pointer to the tip node
 Node* Tip_node_pt = 0;

 /// Initial tip velocity
 double U_tip = 0.0;
 Data* Tip_velocity_pt = 0;

 /// Magnitude of the slope at the tip
 double Slope_magn = 10000.0;

 /// Ohnesorg number
 double Oh = 1.0;

 /// Weber number
 double We = 1.0;

 /// Bond number
 double Bo = 1.0;

 /// Body force function
 void body_force_function(const double& time, double& force)
 {
  force=Bo;
 }
 // The initial profile h(z)
 double initial_profile(const double& z)
  {
   return R*sqrt(Slope_magn/(Slope_magn*Z_tip - R))*
    sqrt(Z_tip - z + R*R/(4.0*Slope_magn*(Slope_magn*Z_tip - R))) -
    R*R/(2.0*(Slope_magn*Z_tip - R));
  }

 // The initial slope dhdz(z)
 double initial_slope(const double& z)
  {
   return -R/2.0*sqrt(Slope_magn/(Slope_magn*Z_tip - R))/
    sqrt(Z_tip - z + R*R/(4.0*Slope_magn*(Slope_magn*Z_tip - R)));
  }

  /// Initial drop volume
  double V_0 = 0.0;

  /// Function for the current volume
  double drop_volume(const double& t)
  {
    return V_0 + MathematicalConstants::Pi*R*R*sqrt(We)*t;
  }

  double Prescribed_drop_volume = 0.0;

  /// Trace file
  ofstream Trace_file;

} // end of namespace


} //End of namespace extension

//==start_of_problem_class============================================
/// 1D axisymmetric thin film DrippingFaucet problem
//====================================================================
template<class ELEMENT> 
class AxisymmetricThinFilmDrippingFaucetProblem : public Problem
{

public:

 /// Constructor: Pass number of elements and pointer to body force function
 AxisymmetricThinFilmDrippingFaucetProblem
 (const unsigned& n_element, 
  AxisymmetricThinFilmDrippingFaucetEquations::AxisymmetricThinFilmDrippingFaucetBodyForceFctPt body_force_fct_pt);

 /// Destructor (empty)
 ~AxisymmetricThinFilmDrippingFaucetProblem()
  {
    delete Constraint_element_pt->internal_data_pt(0);
    delete Constraint_element_pt->internal_data_pt(1);
    delete Constraint_element_pt;
    delete Constraint_mesh_pt;

    delete Fluid_mesh_pt->spatial_error_estimator_pt();
    delete Fluid_mesh_pt;

    delete this->time_stepper_pt();
  }

 /// Actions before Newton convergence check
 void actions_before_newton_convergence_check()
 {
   Problem_Parameter::Tip_node_pt->x(0) = Problem_Parameter::Tip_position_pt->value(0);
   Problem_Parameter::Tip_node_pt->set_value(1, Problem_Parameter::Tip_velocity_pt->value(0));
 }

 /// \short Actions before implicit timestep: Update the target volume
 void actions_before_implicit_timestep()
 {
   // Current time
   const double t = time_stepper_pt()->time();

   Problem_Parameter::Prescribed_drop_volume = Problem_Parameter::drop_volume(t);
 }

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 /// \short Global error norm for adaptive timestepping
 double global_temporal_error_norm()
 {
    double global_error = 0.0;
    unsigned count=0;
    unsigned num_nod = Fluid_mesh_pt->nnode();
    for (unsigned inod = 0; inod < num_nod; inod++)
    {
      // Get node
      Node* nod_pt=Fluid_mesh_pt->node_pt(inod);
       
      // Error based on velocity
      double err=nod_pt->time_stepper_pt()->temporal_error_in_value(nod_pt,1);
       
      // Add the square of the individual error to the global error
      count++;
      global_error += err*err;
    }

    // Divide by the number of nodes
    global_error /= double(count);
     
    // Return square root...
    global_error=sqrt(global_error);

    return global_error;
  }

 /// \short Doc the solution
 void doc_solution();

 /// Set the tip node
 void set_tip_node_element_pt()
 {
  Problem_Parameter::Tip_node_pt = dynamic_cast<Node*>(Fluid_mesh_pt->boundary_node_pt(1,0));

  unsigned n_element = Fluid_mesh_pt->nelement();
  for(unsigned i=0;i<n_element;i++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(i));
    // Find the element that contains the tip node
    // Loop over all nodes
    unsigned n_node = elem_pt->nnode();
    for(unsigned inod=0; inod<n_node; inod++)
    {
      Node* nod_pt = elem_pt->node_pt(inod);
      if (nod_pt == Problem_Parameter::Tip_node_pt)
      {
        Tip_element_pt = elem_pt;
      }
    }
  }
 }

 void actions_before_adapt()
 {
   delete_constraint_element();
   rebuild_global_mesh();
 }

 void actions_after_adapt()
 {
   set_tip_node_element_pt();
   create_constraint_element();
   complete_problem_setup();
   rebuild_global_mesh();
 }

 void complete_problem_setup()
 {
  // Set the boundary conditions for this problem: By default, all nodal
  // values are free -- we only need to pin the ones that have 
  // Dirichlet conditions.
  Fluid_mesh_pt->boundary_node_pt(0,0)->pin(0);
  Fluid_mesh_pt->boundary_node_pt(0,0)->set_value(0, Problem_Parameter::R);
  Fluid_mesh_pt->boundary_node_pt(0,0)->pin(1);
  Fluid_mesh_pt->boundary_node_pt(0,0)->set_value(1, sqrt(Problem_Parameter::We));

  Problem_Parameter::Tip_node_pt->x(0) = Problem_Parameter::Tip_position_pt->value(0);
  Problem_Parameter::Tip_node_pt->pin(0);
  Problem_Parameter::Tip_node_pt->set_value(0, 0.0);
  Problem_Parameter::Tip_node_pt->pin(1);
  Problem_Parameter::Tip_node_pt->set_value(1, Problem_Parameter::Tip_velocity_pt->value(0));
  Problem_Parameter::Tip_node_pt->pin(2);
  Problem_Parameter::Tip_node_pt->set_value(2, -Problem_Parameter::Slope_magn);

  // Complete the setup of the 1D axisymmetric thin film DrippingFaucet problem:

  // Loop over elements and set pointers to Ohnesorg number 
  // and body force function
  unsigned n_element = Fluid_mesh_pt->nelement();
  for(unsigned i=0;i<n_element;i++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(i));
   
    // Set the Ohnesorg number
    elem_pt->oh_pt() = &Problem_Parameter::Oh;

    //Set the body force function pointer
    elem_pt->body_force_fct_pt() = Body_force_fct_pt;
  }
  // Add tip position and velocity as external data to the tip element
  Tip_element_pt->add_external_data(Problem_Parameter::Tip_position_pt);
  Tip_element_pt->add_external_data(Problem_Parameter::Tip_velocity_pt);
 }

private:

 /// Compute the volume of the drop
 double calculate_drop_volume()
 {
  double volume=0.0;
  unsigned n_element = Fluid_mesh_pt->nelement();
  for(unsigned i=0;i<n_element;i++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(i));

    volume += elem_pt->compute_physical_size();
  }
  return volume;
 }

 /// Create constraint element
 void create_constraint_element()
 {
    Constraint_element_pt = new DrippingFaucetConstraintElement<ELEMENT>(Problem_Parameter::Z_tip,
      Problem_Parameter::U_tip, Fluid_mesh_pt, Problem_Parameter::Tip_node_pt,
      &Problem_Parameter::Prescribed_drop_volume);
    // Set pointer to position and velocity data to allow external access
    Problem_Parameter::Tip_position_pt = Constraint_element_pt->internal_data_pt(0);
    Problem_Parameter::Tip_velocity_pt = Constraint_element_pt->internal_data_pt(1);
    // Add element to mesh
    Constraint_mesh_pt->add_element_pt(Constraint_element_pt);
 }

 /// Delete constraint element
 void delete_constraint_element()
 {
   // Back up the position and velocity of the tip
   Problem_Parameter::Z_tip = Problem_Parameter::Tip_position_pt->value(0);
   Problem_Parameter::U_tip = Problem_Parameter::Tip_velocity_pt->value(0);
   // Remove the external data from the tip element
   Tip_element_pt->flush_external_data();

   delete Constraint_mesh_pt->element_pt(0);

   Constraint_mesh_pt->flush_element_and_node_storage();
 }

 /// Mesh pointer
 RefineableOneDMesh<ELEMENT>* Fluid_mesh_pt;

 /// Constraint element
 DrippingFaucetConstraintElement<ELEMENT>* Constraint_element_pt;

 /// Constraint mesh
 Mesh* Constraint_mesh_pt;

 /// Pointer to the element at the tip
 ELEMENT* Tip_element_pt;

 /// Pointer to body force function
 AxisymmetricThinFilmDrippingFaucetEquations::AxisymmetricThinFilmDrippingFaucetBodyForceFctPt Body_force_fct_pt;

}; // end of problem class

//=====start_of_constructor===============================================
/// \short Constructor for 1D axisymmetric thin film DrippingFaucet problem.
/// Discretise the 1D domain with n_element elements of type ELEMENT.
/// Specify function pointer to body force function. 
//========================================================================
template<class ELEMENT>
AxisymmetricThinFilmDrippingFaucetProblem<ELEMENT>::AxisymmetricThinFilmDrippingFaucetProblem
(const unsigned& n_element,
 AxisymmetricThinFilmDrippingFaucetEquations::AxisymmetricThinFilmDrippingFaucetBodyForceFctPt body_force_fct_pt) : 
 Body_force_fct_pt(body_force_fct_pt)
{ 
 // Problem::Sparse_assembly_method = Perform_assembly_using_two_arrays;

 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps.
 this->add_time_stepper_pt(new BDF<2>(false));
 oomph_info << "Using BDF2\n";

 //linear_solver_pt()=new FD_LU;
 Max_residuals = 1.0e4;
 Max_newton_iterations=20;
 //Newton_solver_tolerance = 1.0e-6;
 //enable_globally_convergent_newton_method();

 // Build mesh and store pointer in Problem
 Fluid_mesh_pt = new RefineableOneDMesh<ELEMENT>(n_element,Problem_Parameter::Z_tip,
						 this->time_stepper_pt());

 // Create/set error estimator
 Fluid_mesh_pt->spatial_error_estimator_pt() = new Z2ErrorEstimator;

 // Set targets for spatial adaptivity
 Fluid_mesh_pt->max_permitted_error()=5.0e-6;
 Fluid_mesh_pt->min_permitted_error()=1.0e-6;

 // Set the tip node and element
 set_tip_node_element_pt();

 /// Create mesh for contraints
 Constraint_mesh_pt = new Mesh();

 /// Create constraint element
 create_constraint_element();

 // Set the initial condition
 unsigned n_node = Fluid_mesh_pt->nnode();
 for(unsigned inod=0;inod<n_node;inod++)
  {
   // Get pointer to node
   Node* nod_pt = dynamic_cast<Node*>(Fluid_mesh_pt->node_pt(inod));
   double z=nod_pt->x(0);
   double h=Problem_Parameter::initial_profile(z);
   nod_pt->set_value(0,h);
   double dhdz=Problem_Parameter::initial_slope(z);
   nod_pt->set_value(2,dhdz);
  }

 // Complete the setup of the 1D axisymmetric thin film DrippingFaucet problem:
 complete_problem_setup();
  // Set the initial drop volume
  Problem_Parameter::V_0 = calculate_drop_volume();
  Problem_Parameter::Prescribed_drop_volume = Problem_Parameter::V_0;

 // Combine meshes
 add_sub_mesh(Fluid_mesh_pt);
 add_sub_mesh(Constraint_mesh_pt);
 build_global_mesh();

 // Setup equation numbering scheme
 oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;

} // end of constructor


//===start_of_doc=========================================================
/// Doc the solution in tecplot format. Label files with label.
//========================================================================
template<class ELEMENT>
void AxisymmetricThinFilmDrippingFaucetProblem<ELEMENT>::doc_solution()
{ 
 oomph_info << "Docing step: " 
	    << Problem_Parameter::Doc_info.number() << std::endl;

 ofstream some_file;
 char filename[1000];

 // Number of plot points
 unsigned npts=5;

 sprintf(filename,"%s/soln%i.dat",
	 Problem_Parameter::Doc_info.directory().c_str(),
	 Problem_Parameter::Doc_info.number());
 some_file.open(filename);
 Fluid_mesh_pt->output(some_file,npts);
 some_file.close();

 // Write trace file
 Problem_Parameter::Trace_file
   << this->time_pt()->time() << " "
   << Problem_Parameter::Tip_node_pt->x(0) << " "
   << Problem_Parameter::Tip_node_pt->value(1) << " "
   << calculate_drop_volume() << " "
   << Fluid_mesh_pt->nelement() << " "
   << Problem_Parameter::Doc_info.number() << " "
   << std::endl;

 // Increment the doc_info number
 Problem_Parameter::Doc_info.number()++;

} // end of doc

 

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//======start_of_main==================================================
/// Driver for 1D Poisson problem
//=====================================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // Sometimes the Newton method produces inverted elements while
 // it's iterating towards a solution with non-inverted elements:
 // accept these!
 FiniteElement::Accept_negative_jacobian = true;

 // Define possible command line arguments and parse the ones that
 // were actually specified
  
 // Output directory
 CommandLineArgs::specify_command_line_flag(
   "--dir", &Problem_Parameter::Directory);

 // Constant timestep
 double dt=0.01;
 CommandLineArgs::specify_command_line_flag("--dt", &dt);

 // Do timestepping until tmax
 double t_max=10.0;
 CommandLineArgs::specify_command_line_flag("--t_max", &t_max);

 // Number of elements
 unsigned n_element=20;
 CommandLineArgs::specify_command_line_flag("--n_element", &n_element);

 // Adaptation?
 unsigned nadapt_interval=0; 
 CommandLineArgs::specify_command_line_flag("--adapt",&nadapt_interval);

 // Ohnesorg number
 CommandLineArgs::specify_command_line_flag("--oh",&Problem_Parameter::Oh);

 // Weber number
 CommandLineArgs::specify_command_line_flag("--we",&Problem_Parameter::We);

 // Bond number
 CommandLineArgs::specify_command_line_flag("--bo",&Problem_Parameter::Bo);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Set up the problem: 
 // Solve a 1D axisymmetric thin film DrippingFaucet problem using a body force
 AxisymmetricThinFilmDrippingFaucetProblem<MyRefinableElement<4> >
   //Element type as template parameter
   problem(n_element,Problem_Parameter::body_force_function);

 // Check whether the problem can be solved
 cout << "\n\n\nProblem self-test ";
 if (problem.self_test()==0)  
  {
   cout << "passed: Problem can be solved." << std::endl;
  }
 else 
  {
   throw OomphLibError("failed!",
                       OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
  }
 // Open trace file
 Problem_Parameter::Trace_file.open((Problem_Parameter::Directory+
                                     "/trace.dat").c_str());

 // Write trace file
 Problem_Parameter::Trace_file  
  << "VARIABLES="
  << "\"time\", "  // 1
  << "\"z_tip\"," // 2
  << "\"u_tip\"," // 3
  << "\"volume\"," // 4
  << "\"n_element\"," // 5
  << "\"doc number\"" // 6
  << std::endl; 

 // Output directory
 Problem_Parameter::Doc_info.set_directory(Problem_Parameter::Directory);

 // Initialise timestepper
 problem.initialise_dt(dt);
   
 // Perform impulsive start from current state
 problem.assign_initial_values_impulsive();

 //Output initial condition
 problem.doc_solution();

 // Refine problem uniformly 2 times (to check automatic unrefinement)
  problem.unsteady_newton_solve(dt);
  problem.doc_solution();
 for(unsigned i=0;i<3;i++) 
 { 
   problem.refine_uniformly();
   problem.unsteady_newton_solve(dt);
   problem.doc_solution();
  }

 while(problem.time_pt()->time() < t_max)
  {
   // Adapt?
   if( CommandLineArgs::command_line_flag_has_been_set("--adapt") &&
       ((Problem_Parameter::Doc_info.number())%nadapt_interval==0))
    {
     oomph_info << "obacht: time for another adaptation at doc number ="
		<< Problem_Parameter::Doc_info.number() << std::endl;
     problem.adapt();
    }
   problem.unsteady_newton_solve(dt);
   problem.doc_solution();
  }
} // end of main


