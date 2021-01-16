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
#include "num_rec.h"

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
//==start_of_namespace================================================
/// Namespace for problem parameters
//====================================================================
namespace Problem_Parameter
{
 /// Output directory
 std::string Directory = "RESLT";

 DocInfo Doc_info;

 /// Initial length
 double L_start = 1.0;

 /// Ohnesorg number
 double Oh = 1.0;

 /// Weber number
 double We = 1.0;

 /// Body force function
 void body_force_function(const double& time, double& force)
 {
  force=-1.0;
 }
 // The initial profile h(z)
 double initial_profile(const double& z)
  {
   if(z <= 1.0)
    {
     return cos(MathematicalConstants::Pi*z);
    }
   else
    {
     return 0.0;
    }
  }

 // The initial slope dhdz(z)
 double initial_slope(const double& z)
  {
   if(z <= 1.0)
    {
     return -MathematicalConstants::Pi*sin(MathematicalConstants::Pi*z);
    }
   else
    {
     return 0.0;
    }
  }

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
   delete mesh_pt();
  }

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 /// \short Doc the solution, pass the number of the case considered,
 /// so that output files can be distinguished.
 void doc_solution();

private:

 /// Mesh pointer
 RefineableOneDMesh<ELEMENT>* Fluid_mesh_pt;

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
 Fluid_mesh_pt = new RefineableOneDMesh<ELEMENT>(n_element,Problem_Parameter::L_start,
						 this->time_stepper_pt());

 // Create/set error estimator
 Fluid_mesh_pt->spatial_error_estimator_pt() = new Z2ErrorEstimator;

 // Set targets for spatial adaptivity
 Fluid_mesh_pt->max_permitted_error()=5.0e-6;
 Fluid_mesh_pt->min_permitted_error()=1.0e-6;

 Problem::mesh_pt() = Fluid_mesh_pt;

 // Set the boundary conditions for this problem: By default, all nodal
 // values are free -- we only need to pin the ones that have 
 // Dirichlet conditions. 

 // Pin the single nodal value at the single node on mesh 
 // boundary 0 (= the left domain boundary at x=0)
 //mesh_pt()->boundary_node_pt(0,0)->pin(0);
 
 // Pin the single nodal value at the single node on mesh 
 // boundary 1 (= the right domain boundary at x=1)
 //mesh_pt()->boundary_node_pt(1,0)->pin(0);
 Fluid_mesh_pt->boundary_node_pt(0,0)->pin(1);
 Fluid_mesh_pt->boundary_node_pt(0,0)->set_value(1, 0.0);
 Fluid_mesh_pt->boundary_node_pt(1,0)->pin(0);
 Fluid_mesh_pt->boundary_node_pt(1,0)->set_value(0, 0.0);

 // Set the initial condition
 unsigned n_node = mesh_pt()->nnode();
 for(unsigned inod=0;inod<n_node;inod++)
  {
   // Get pointer to node
   Node* nod_pt = dynamic_cast<Node*>(mesh_pt()->node_pt(inod));
   double z=nod_pt->x(0);
   double h=Problem_Parameter::initial_profile(z);
   nod_pt->set_value(0,h);
   double dhdz=Problem_Parameter::initial_slope(z);
   nod_pt->set_value(2,dhdz);
  }

 // Complete the setup of the 1D axisymmetric thin film DrippingFaucet problem:

 // Loop over elements and set pointers to DrippingFaucet number 
 // and body force function
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   // Set the DrippingFaucet number
   elem_pt->oh_pt() = &Problem_Parameter::Oh;

   //Set the body force function pointer
   elem_pt->body_force_fct_pt() = Body_force_fct_pt;
  }

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

 // Compute drop volume
 double current_vol=0.0;
 double nel=Fluid_mesh_pt->nelement();
 for(unsigned e=0;e<nel;e++)
  {
   // Get element
   ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));

   // Add to physical size (actual volume)
   current_vol+=el_pt->compute_physical_size();
  }

 // Write trace file
 Problem_Parameter::Trace_file
   << this->time_pt()->time() << " "
   << Fluid_mesh_pt->boundary_node_pt(0,0)->value(0) << " "
   << current_vol << " "
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

 // Define possible command line arguments and parse the ones that
 // were actually specified
  
 // Output directory
 CommandLineArgs::specify_command_line_flag(
   "--dir", &Problem_Parameter::Directory);

 // Constant timestep
 double dt=1.0/40.0;
 CommandLineArgs::specify_command_line_flag("--dt", &dt);

 // Do timestepping until tmax
 double t_max=500.0;
 CommandLineArgs::specify_command_line_flag("--t_max", &t_max);

 // Number of elements
 unsigned n_element=120;
 CommandLineArgs::specify_command_line_flag("--n_element", &n_element);

 // Adaptation?
 unsigned nadapt_interval=0; 
 CommandLineArgs::specify_command_line_flag("--adapt",&nadapt_interval);

 // Ohnesorg number
 CommandLineArgs::specify_command_line_flag("--oh",&Problem_Parameter::Oh);

 // Weber number
 CommandLineArgs::specify_command_line_flag("--we",&Problem_Parameter::We);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Set up the problem: 
 // Solve a 1D axisymmetric thin film DrippingFaucet problem using a body force
 AxisymmetricThinFilmDrippingFaucetProblem<RefineableAxisymmetricThinFilmDrippingFaucetElement<4> >
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
  << "\"height\"," // 2
  << "\"volume\"," // 3
  << "\"n_element\"," // 4
  << "\"doc number\"" // 5
  << std::endl; 

 // Output directory
 Problem_Parameter::Doc_info.set_directory(Problem_Parameter::Directory);

 // Initialise timestepper
 problem.initialise_dt(dt);
   
 // Perform impulsive start from current state
 problem.assign_initial_values_impulsive();

 // Refine problem uniformly 4 times (to check automatic unrefinement)
 for(unsigned i=0;i<3;i++) { problem.refine_uniformly(); }

 //Output initial condition
 problem.doc_solution();

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


