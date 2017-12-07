// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Mon 27 Nov 2017 21:10:41

/**
 * @file MDM_mass_eigenstates.cpp
 * @brief implementation of the MDM model class
 *
 * Contains the definition of the MDM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Mon 27 Nov 2017 21:10:41 with FlexibleSUSY
 * 1.7.4 (git commit: unknown) and SARAH 4.12.2 .
 */

#include "MDM_mass_eigenstates.hpp"
#include "eigen_utils.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"
#include "numerics2.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "root_finder.hpp"
#include "fixed_point_iterator.hpp"
#include "gsl_utils.hpp"
#include "config.h"
#include "parallel.hpp"
#include "pv.hpp"
#include "functors.hpp"




#include <cmath>
#include <iostream>
#include <memory>
#include <algorithm>

#include <gsl/gsl_multiroots.h>

namespace flexiblesusy {

#define CLASSNAME MDM_mass_eigenstates

#define PHYSICAL(parameter) physical.parameter
#define INPUT(parameter) model->get_input().parameter
#define LOCALINPUT(parameter) input.parameter
#define MODELPARAMETER(parameter) model->get_##parameter()

#define HIGGS_2LOOP_CORRECTION_AT_AS     two_loop_corrections.higgs_at_as
#define HIGGS_2LOOP_CORRECTION_AB_AS     two_loop_corrections.higgs_ab_as
#define HIGGS_2LOOP_CORRECTION_AT_AT     two_loop_corrections.higgs_at_at
#define HIGGS_2LOOP_CORRECTION_ATAU_ATAU two_loop_corrections.higgs_atau_atau
#define TOP_POLE_QCD_CORRECTION          two_loop_corrections.top_qcd
#define HIGGS_3LOOP_CORRECTION_AT_AS_AS  1

CLASSNAME::MDM_mass_eigenstates(const MDM_input_parameters& input_)
   : MDM_soft_parameters(input_)
   , number_of_ewsb_iterations(100)
   , number_of_mass_iterations(20)
   , ewsb_loop_order(2)
   , pole_mass_loop_order(2)
   , calculate_sm_pole_masses(false)
   , calculate_bsm_pole_masses(true)
   , force_output(false)
   , precision(1.0e-3)
   , ewsb_iteration_precision(1.0e-5)
   , physical()
   , problems(MDM_info::particle_names)
   , two_loop_corrections()
   , MVG(0), MHp(0), MFv(Eigen::Array<double,3,1>::Zero()), MFc(0), MFg(0), MFn
      (0), MAh(0), Mhh(0), MFd(Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array
      <double,3,1>::Zero()), MFe(Eigen::Array<double,3,1>::Zero()), MVWp(0), MVP(0
      ), MVZ(0)

   , Vd(Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ud(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), Vu(Eigen::Matrix<std::complex<double>,3,
      3>::Zero()), Uu(Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ve(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ue(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZZ(Eigen::Matrix<double,2,2>::Zero())

   , PhaseFc(1,0), PhaseFn(1,0), PhaseFg(1,0)

{
}

CLASSNAME::~MDM_mass_eigenstates()
{
}

void CLASSNAME::do_calculate_sm_pole_masses(bool flag)
{
   calculate_sm_pole_masses = flag;
}

bool CLASSNAME::do_calculate_sm_pole_masses() const
{
   return calculate_sm_pole_masses;
}

void CLASSNAME::do_calculate_bsm_pole_masses(bool flag)
{
   calculate_bsm_pole_masses = flag;
}

bool CLASSNAME::do_calculate_bsm_pole_masses() const
{
   return calculate_bsm_pole_masses;
}

void CLASSNAME::do_force_output(bool flag)
{
   force_output = flag;
}

bool CLASSNAME::do_force_output() const
{
   return force_output;
}

void CLASSNAME::set_ewsb_loop_order(unsigned loop_order)
{
   ewsb_loop_order = loop_order;
}

void CLASSNAME::set_two_loop_corrections(const Two_loop_corrections& two_loop_corrections_)
{
   two_loop_corrections = two_loop_corrections_;
}

const Two_loop_corrections& CLASSNAME::get_two_loop_corrections() const
{
   return two_loop_corrections;
}

void CLASSNAME::set_number_of_ewsb_iterations(std::size_t iterations)
{
   number_of_ewsb_iterations = iterations;
}

std::size_t CLASSNAME::get_number_of_ewsb_iterations() const
{
   return number_of_ewsb_iterations;
}

void CLASSNAME::set_number_of_mass_iterations(std::size_t iterations)
{
   number_of_mass_iterations = iterations;
}

std::size_t CLASSNAME::get_number_of_mass_iterations() const
{
   return number_of_mass_iterations;
}

void CLASSNAME::set_precision(double precision_)
{
   precision = precision_;
   ewsb_iteration_precision = precision_;
}

void CLASSNAME::set_pole_mass_loop_order(unsigned loop_order)
{
   pole_mass_loop_order = loop_order;
}

unsigned CLASSNAME::get_pole_mass_loop_order() const
{
   return pole_mass_loop_order;
}

void CLASSNAME::set_ewsb_iteration_precision(double precision)
{
   ewsb_iteration_precision = precision;
}

double CLASSNAME::get_ewsb_iteration_precision() const
{
   return ewsb_iteration_precision;
}

double CLASSNAME::get_precision() const
{
   return precision;
}

double CLASSNAME::get_ewsb_loop_order() const
{
   return ewsb_loop_order;
}

const MDM_physical& CLASSNAME::get_physical() const
{
   return physical;
}

MDM_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const MDM_physical& physical_)
{
   physical = physical_;
}

const Problems<MDM_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems() const
{
   return problems;
}

Problems<MDM_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems()
{
   return problems;
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @return array of tadpoles
 */
Eigen::Matrix<double, CLASSNAME::number_of_ewsb_equations, 1> CLASSNAME::tadpole_equations() const
{
   Eigen::Matrix<double, number_of_ewsb_equations, 1> tadpole(Eigen::Matrix<double, number_of_ewsb_equations, 1>::Zero());

   tadpole[0] = get_ewsb_eq_hh_1();

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh());

      if (ewsb_loop_order > 1) {

      }
   }

   return tadpole;
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @param tadpole array of tadpole
 */
void CLASSNAME::tadpole_equations(double tadpole[number_of_ewsb_equations]) const
{
   const auto tadpole_(tadpole_equations());
   std::copy(tadpole_.data(), tadpole_.data() + number_of_ewsb_equations, tadpole);
}

/**
 * Method which calculates the tadpoles at loop order specified in the
 * pointer to the CLASSNAME::EWSB_args struct.
 *
 * @param x GSL vector of EWSB output parameters
 * @param params pointer to CLASSNAME::EWSB_args struct
 * @param f GSL vector with tadpoles
 *
 * @return GSL_EDOM if x contains Nans, GSL_SUCCESS otherwise.
 */
int CLASSNAME::tadpole_equations(const gsl_vector* x, void* params, gsl_vector* f)
{
   if (!is_finite(x)) {
      gsl_vector_set_all(f, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   const CLASSNAME::EWSB_args* ewsb_args
      = static_cast<CLASSNAME::EWSB_args*>(params);
   MDM_mass_eigenstates* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   model->set_LamH(gsl_vector_get(x, 0));


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   const auto tadpole(model->tadpole_equations());

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, tadpole[i]);

   return IsFinite(tadpole) ? GSL_SUCCESS : GSL_EDOM;
}

/**
 * This method solves the EWSB conditions iteratively, trying several
 * root finding methods until a solution is found.
 */
int CLASSNAME::solve_ewsb_iteratively()
{
   EWSB_args params = {this, ewsb_loop_order};

   std::unique_ptr<EWSB_solver> solvers[] = {
      std::unique_ptr<EWSB_solver>(new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_relative>(CLASSNAME::ewsb_step, &params, number_of_ewsb_iterations, fixed_point_iterator::Convergence_tester_relative(ewsb_iteration_precision))),
      std::unique_ptr<EWSB_solver>(new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, gsl_multiroot_fsolver_hybrids)),
      std::unique_ptr<EWSB_solver>(new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, gsl_multiroot_fsolver_broyden))
   };

   const std::size_t number_of_solvers = sizeof(solvers)/sizeof(*solvers);
   const auto x_init(ewsb_initial_guess());

   VERBOSE_MSG("Solving EWSB equations ...");
   VERBOSE_MSG("\tInitial guess: x_init = " << x_init.transpose());

   int status;
   for (std::size_t i = 0; i < number_of_solvers; ++i) {
      VERBOSE_MSG("\tStarting EWSB iteration using solver " << i);
      status = solve_ewsb_iteratively_with(solvers[i].get(), x_init);
      if (status == EWSB_solver::SUCCESS) {
         VERBOSE_MSG("\tSolver " << i << " finished successfully!");
         break;
      }
#ifdef ENABLE_VERBOSE
      else {
         WARNING("\tSolver " << i << " could not find a solution!"
                 " (requested precision: " << ewsb_iteration_precision << ")");
      }
#endif
   }

   if (status == EWSB_solver::SUCCESS) {
      problems.unflag_no_ewsb();
   } else {
      problems.flag_no_ewsb();
#ifdef ENABLE_VERBOSE
      WARNING("\tCould not find a solution to the EWSB equations!"
              " (requested precision: " << ewsb_iteration_precision << ")");
#endif
   }

   return status;
}

/**
 * Solves EWSB equations with given EWSB solver
 *
 * @param solver EWSB solver
 * @param x_init initial values
 *
 * @return status of the EWSB solver
 */
int CLASSNAME::solve_ewsb_iteratively_with(
   EWSB_solver* solver,
   const Eigen::Matrix<double, number_of_ewsb_equations, 1>& x_init
)
{
   const int status = solver->solve(&x_init[0]);

   LamH = solver->get_solution(0);


   return status;
}

int CLASSNAME::solve_ewsb_iteratively(unsigned loop_order)
{
   // temporarily set `ewsb_loop_order' to `loop_order' and do
   // iteration
   const unsigned old_loop_order = ewsb_loop_order;
   ewsb_loop_order = loop_order;
   const int status = solve_ewsb_iteratively();
   ewsb_loop_order = old_loop_order;
   return status;
}


int CLASSNAME::solve_ewsb_tree_level()
{
   int error = 0;

   const double old_LamH = LamH;

   LamH = Re(mu2/Sqr(v));

   const bool is_finite = IsFinite(LamH);

   if (!is_finite) {
      LamH = old_LamH;
      error = 1;
   }


   return error;
}

int CLASSNAME::solve_ewsb_tree_level_custom()
{
   int error = 0;



   return error;
}

int CLASSNAME::solve_ewsb_one_loop()
{
   return solve_ewsb_iteratively(1);
}

int CLASSNAME::solve_ewsb()
{
   VERBOSE_MSG("\tSolving EWSB at " << ewsb_loop_order << "-loop order");

   if (ewsb_loop_order == 0)
      return solve_ewsb_tree_level();

   return solve_ewsb_iteratively(ewsb_loop_order);
}

Eigen::Matrix<double, CLASSNAME::number_of_ewsb_equations, 1> CLASSNAME::ewsb_initial_guess()
{
   Eigen::Matrix<double, number_of_ewsb_equations, 1> x_init(Eigen::Matrix<double, number_of_ewsb_equations, 1>::Zero());

   x_init[0] = LamH;


   return x_init;
}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * Throws exception of type EEWSBStepFailed if new EWSB parameters are
 * inf or nan.
 *
 * @return new set of EWSB output parameters
 */
Eigen::Matrix<double, CLASSNAME::number_of_ewsb_equations, 1> CLASSNAME::ewsb_step() const
{
   double tadpole[number_of_ewsb_equations] = { 0. };
   Eigen::Matrix<double, number_of_ewsb_equations, 1> ewsb_parameters(Eigen::Matrix<double, number_of_ewsb_equations, 1>::Zero());

   if (ewsb_loop_order > 0) {
      tadpole[0] += Re(tadpole_hh());

      if (ewsb_loop_order > 1) {

      }
   }

   double LamH;

   LamH = Re((mu2*v + tadpole[0])/Power(v,3));

   const bool is_finite = IsFinite(LamH);


   if (!is_finite)
      throw EEWSBStepFailed();

   ewsb_parameters[0] = LamH;


   return ewsb_parameters;
}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * @param x old EWSB output parameters
 * @param params further function parameters
 * @param f new EWSB output parameters
 *
 * @return Returns status of CLASSNAME::ewsb_step
 */
int CLASSNAME::ewsb_step(const gsl_vector* x, void* params, gsl_vector* f)
{
   if (!is_finite(x)) {
      gsl_vector_set_all(f, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   const CLASSNAME::EWSB_args* ewsb_args
      = static_cast<CLASSNAME::EWSB_args*>(params);
   MDM_mass_eigenstates* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   const double LamH = gsl_vector_get(x, 0);

   model->set_LamH(LamH);


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   Eigen::Matrix<double, number_of_ewsb_equations, 1> ewsb_parameters;
   ewsb_parameters[0] = LamH;


   int status = GSL_SUCCESS;

   try {
      ewsb_parameters = model->ewsb_step();
      status = GSL_SUCCESS;
   } catch (...) {
      status = GSL_EDOM;
   }

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, ewsb_parameters[i]);

   return status;
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "MDM\n"
           "========================================\n";
   MDM_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MHp = " << MHp << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MFc = " << MFc << '\n';
   ostr << "MFg = " << MFg << '\n';
   ostr << "MFn = " << MFn << '\n';
   ostr << "MAh = " << MAh << '\n';
   ostr << "Mhh = " << Mhh << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MVWp = " << MVWp << '\n';

   ostr << "----------------------------------------\n"
           "tree-level DRbar mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "Vd = " << Vd << '\n';
   ostr << "Ud = " << Ud << '\n';
   ostr << "Vu = " << Vu << '\n';
   ostr << "Uu = " << Uu << '\n';
   ostr << "Ve = " << Ve << '\n';
   ostr << "Ue = " << Ue << '\n';
   ostr << "ZZ = " << ZZ << '\n';

   physical.print(ostr);
}

/**
 * wrapper routines for passarino Veltman functions
 */

double CLASSNAME::A0(double m) const
{
   return passarino_veltman::ReA0(m*m, Sqr(get_scale()));
}

double CLASSNAME::B0(double p, double m1, double m2) const
{
   return passarino_veltman::ReB0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::B1(double p, double m1, double m2) const
{
   return passarino_veltman::ReB1(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::B00(double p, double m1, double m2) const
{
   return passarino_veltman::ReB00(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::B22(double p, double m1, double m2) const
{
   return passarino_veltman::ReB22(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::H0(double p, double m1, double m2) const
{
   return passarino_veltman::ReH0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::F0(double p, double m1, double m2) const
{
   return passarino_veltman::ReF0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::G0(double p, double m1, double m2) const
{
   return passarino_veltman::ReG0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

/**
 * routine which finds the DRbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_masses()
{
   const auto old_LamH = LamH;

   solve_ewsb_tree_level();

   calculate_MVPVZ();
   calculate_MVWp();
   calculate_MFe();
   calculate_MFu();
   calculate_MFd();
   calculate_Mhh();
   calculate_MAh();
   calculate_MFn();
   calculate_MFg();
   calculate_MFc();
   calculate_MFv();
   calculate_MHp();
   calculate_MVG();

   LamH = old_LamH;

}

/**
 * Backward compatibility routine which finds the DRbar mass
 * eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_parameters()
{
   calculate_DRbar_masses();
}

/**
 * routine which finds the pole mass eigenstates and mixings.
 */
void CLASSNAME::calculate_pole_masses()
{
#ifdef ENABLE_THREADS
   auto obj_ptr = this;

   std::future<void> fut_MVG;
   std::future<void> fut_MFv;
   std::future<void> fut_MFc;
   std::future<void> fut_MFg;
   std::future<void> fut_MFn;
   std::future<void> fut_Mhh;
   std::future<void> fut_MVP;
   std::future<void> fut_MVZ;
   std::future<void> fut_MFd;
   std::future<void> fut_MFu;
   std::future<void> fut_MFe;
   std::future<void> fut_MVWp;

   if (calculate_bsm_pole_masses) {
      fut_MFc = run_async([obj_ptr] () { obj_ptr->calculate_MFc_pole(); });
      fut_MFg = run_async([obj_ptr] () { obj_ptr->calculate_MFg_pole(); });
      fut_MFn = run_async([obj_ptr] () { obj_ptr->calculate_MFn_pole(); });
   }

   if (calculate_sm_pole_masses) {
      fut_MVG = run_async([obj_ptr] () { obj_ptr->calculate_MVG_pole(); });
      fut_MFv = run_async([obj_ptr] () { obj_ptr->calculate_MFv_pole(); });
      fut_Mhh = run_async([obj_ptr] () { obj_ptr->calculate_Mhh_pole(); });
      fut_MVP = run_async([obj_ptr] () { obj_ptr->calculate_MVP_pole(); });
      fut_MVZ = run_async([obj_ptr] () { obj_ptr->calculate_MVZ_pole(); });
      fut_MFd = run_async([obj_ptr] () { obj_ptr->calculate_MFd_pole(); });
      fut_MFu = run_async([obj_ptr] () { obj_ptr->calculate_MFu_pole(); });
      fut_MFe = run_async([obj_ptr] () { obj_ptr->calculate_MFe_pole(); });
      fut_MVWp = run_async([obj_ptr] () { obj_ptr->calculate_MVWp_pole(); });
   }

   if (fut_MFc.valid()) fut_MFc.get();
   if (fut_MFg.valid()) fut_MFg.get();
   if (fut_MFn.valid()) fut_MFn.get();
   if (fut_MVG.valid()) fut_MVG.get();
   if (fut_MFv.valid()) fut_MFv.get();
   if (fut_Mhh.valid()) fut_Mhh.get();
   if (fut_MVP.valid()) fut_MVP.get();
   if (fut_MVZ.valid()) fut_MVZ.get();
   if (fut_MFd.valid()) fut_MFd.get();
   if (fut_MFu.valid()) fut_MFu.get();
   if (fut_MFe.valid()) fut_MFe.get();
   if (fut_MVWp.valid()) fut_MVWp.get();

#else
   if (calculate_bsm_pole_masses) {
      calculate_MFc_pole();
      calculate_MFg_pole();
      calculate_MFn_pole();
   }

   if (calculate_sm_pole_masses) {
      calculate_MVG_pole();
      calculate_MFv_pole();
      calculate_Mhh_pole();
      calculate_MVP_pole();
      calculate_MVZ_pole();
      calculate_MFd_pole();
      calculate_MFu_pole();
      calculate_MFe_pole();
      calculate_MVWp_pole();
   }

#endif
}

void CLASSNAME::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(MVG) = MVG;
   PHYSICAL(MHp) = MHp;
   PHYSICAL(MFv) = MFv;
   PHYSICAL(MFc) = MFc;
   PHYSICAL(MFg) = MFg;
   PHYSICAL(MFn) = MFn;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(Mhh) = Mhh;
   PHYSICAL(MFd) = MFd;
   PHYSICAL(Vd) = Vd;
   PHYSICAL(Ud) = Ud;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(Vu) = Vu;
   PHYSICAL(Uu) = Uu;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(Ve) = Ve;
   PHYSICAL(Ue) = Ue;
   PHYSICAL(MVWp) = MVWp;
   PHYSICAL(MVP) = MVP;
   PHYSICAL(MVZ) = MVZ;

}

/**
 * reorders DRbar masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_DRbar_masses()
{

}

/**
 * reorders pole masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_pole_masses()
{

}

/**
 * Checks the pole masses for tachyons
 */
void CLASSNAME::check_pole_masses_for_tachyons()
{
   if (PHYSICAL(Mhh) < 0.) problems.flag_tachyon(MDM_info::hh);

}

/**
 * calculates spectrum for model once the DRbar parameters at
 * at low energies are known
 */
void CLASSNAME::calculate_spectrum()
{
   calculate_DRbar_masses();
   if (pole_mass_loop_order > 0)
      calculate_pole_masses();

   // move goldstone bosons to the front
   reorder_DRbar_masses();
   if (pole_mass_loop_order == 0)
      copy_DRbar_masses_to_pole_masses();
   else
      reorder_pole_masses();

   check_pole_masses_for_tachyons();

   if (problems.have_problem() && !force_output) {
      clear_DRbar_parameters();
      physical.clear();
   }
}

void CLASSNAME::clear_DRbar_parameters()
{
   MVG = 0.;
   MHp = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MFc = 0.;
   MFg = 0.;
   MFn = 0.;
   MAh = 0.;
   Mhh = 0.;
   MFd = Eigen::Matrix<double,3,1>::Zero();
   Vd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ud = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   Vu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Uu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   Ve = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ue = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MVWp = 0.;
   MVP = 0.;
   MVZ = 0.;

   PhaseFc = std::complex<double>(1.,0.);
   PhaseFn = std::complex<double>(1.,0.);
   PhaseFg = std::complex<double>(1.,0.);

}

void CLASSNAME::clear_problems()
{
   problems.unflag_all_tachyons();
}

void CLASSNAME::clear()
{
   MDM_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_DRbar_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MHp = pars(1);
   MFv(0) = pars(2);
   MFv(1) = pars(3);
   MFv(2) = pars(4);
   MFc = pars(5);
   MFg = pars(6);
   MFn = pars(7);
   MAh = pars(8);
   Mhh = pars(9);
   MFd(0) = pars(10);
   MFd(1) = pars(11);
   MFd(2) = pars(12);
   MFu(0) = pars(13);
   MFu(1) = pars(14);
   MFu(2) = pars(15);
   MFe(0) = pars(16);
   MFe(1) = pars(17);
   MFe(2) = pars(18);
   MVWp = pars(19);
   MVP = pars(20);
   MVZ = pars(21);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses() const
{
   Eigen::ArrayXd pars(22);

   pars(0) = MVG;
   pars(1) = MHp;
   pars(2) = MFv(0);
   pars(3) = MFv(1);
   pars(4) = MFv(2);
   pars(5) = MFc;
   pars(6) = MFg;
   pars(7) = MFn;
   pars(8) = MAh;
   pars(9) = Mhh;
   pars(10) = MFd(0);
   pars(11) = MFd(1);
   pars(12) = MFd(2);
   pars(13) = MFu(0);
   pars(14) = MFu(1);
   pars(15) = MFu(2);
   pars(16) = MFe(0);
   pars(17) = MFe(1);
   pars(18) = MFe(2);
   pars(19) = MVWp;
   pars(20) = MVP;
   pars(21) = MVZ;

   return pars;
}

std::string CLASSNAME::name() const
{
   return "MDM";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   MDM_soft_parameters::run_to(scale, eps);
}







double CLASSNAME::get_mass_matrix_VG() const
{
   const double mass_matrix_VG = Re(0);

   return mass_matrix_VG;
}

void CLASSNAME::calculate_MVG()
{
   const auto mass_matrix_VG = get_mass_matrix_VG();
   MVG = mass_matrix_VG;
}

double CLASSNAME::get_mass_matrix_Hp() const
{
   const double mass_matrix_Hp = Re(-mu2 + LamH*Sqr(v) + 0.25*Sqr(g2)*Sqr
      (v));

   return mass_matrix_Hp;
}

void CLASSNAME::calculate_MHp()
{
   const auto mass_matrix_Hp = get_mass_matrix_Hp();
   MHp = mass_matrix_Hp;

   if (MHp < 0.) {
      problems.flag_tachyon(MDM_info::Hp);
   }

   MHp = AbsSqrt(MHp);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fv() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fv;

   mass_matrix_Fv(0,0) = 0;
   mass_matrix_Fv(0,1) = 0;
   mass_matrix_Fv(0,2) = 0;
   mass_matrix_Fv(1,1) = 0;
   mass_matrix_Fv(1,2) = 0;
   mass_matrix_Fv(2,2) = 0;

   Symmetrize(mass_matrix_Fv);

   return mass_matrix_Fv;
}

void CLASSNAME::calculate_MFv()
{
   MFv.setConstant(0);
}

double CLASSNAME::get_mass_matrix_Fc() const
{
   const double mass_matrix_Fc = Re(2*Yc);

   return mass_matrix_Fc;
}

void CLASSNAME::calculate_MFc()
{
   const auto mass_matrix_Fc = get_mass_matrix_Fc();
   MFc = calculate_dirac_singlet_mass(mass_matrix_Fc, PhaseFc);
}

double CLASSNAME::get_mass_matrix_Fg() const
{
   const double mass_matrix_Fg = Re(2*Yc);

   return mass_matrix_Fg;
}

void CLASSNAME::calculate_MFg()
{
   const auto mass_matrix_Fg = get_mass_matrix_Fg();
   MFg = calculate_dirac_singlet_mass(mass_matrix_Fg, PhaseFg);
}

double CLASSNAME::get_mass_matrix_Fn() const
{
   const double mass_matrix_Fn = Re(2*Yc);

   return mass_matrix_Fn;
}

void CLASSNAME::calculate_MFn()
{
   const auto mass_matrix_Fn = get_mass_matrix_Fn();
   MFn = calculate_majorana_singlet_mass(mass_matrix_Fn, PhaseFn);
}

double CLASSNAME::get_mass_matrix_Ah() const
{
   const double mass_matrix_Ah = Re(0.25*(-4*(mu2 - LamH*Sqr(v)) + Sqr(v)
      *Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))));

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah = get_mass_matrix_Ah();
   MAh = mass_matrix_Ah;

   if (MAh < 0.) {
      problems.flag_tachyon(MDM_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

double CLASSNAME::get_mass_matrix_hh() const
{
   const double mass_matrix_hh = Re(-mu2 + 3*LamH*Sqr(v));

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh = get_mass_matrix_hh();
   Mhh = mass_matrix_hh;

   if (Mhh < 0.) {
      problems.flag_tachyon(MDM_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fd() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = 0.7071067811865475*v*Yd(0,0);
   mass_matrix_Fd(0,1) = 0.7071067811865475*v*Yd(1,0);
   mass_matrix_Fd(0,2) = 0.7071067811865475*v*Yd(2,0);
   mass_matrix_Fd(1,0) = 0.7071067811865475*v*Yd(0,1);
   mass_matrix_Fd(1,1) = 0.7071067811865475*v*Yd(1,1);
   mass_matrix_Fd(1,2) = 0.7071067811865475*v*Yd(2,1);
   mass_matrix_Fd(2,0) = 0.7071067811865475*v*Yd(0,2);
   mass_matrix_Fd(2,1) = 0.7071067811865475*v*Yd(1,2);
   mass_matrix_Fd(2,2) = 0.7071067811865475*v*Yd(2,2);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud, eigenvalue_error);
   problems.flag_bad_mass(MDM_info::Fd, eigenvalue_error > precision *
      Abs(MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fu() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = -0.7071067811865475*v*Yu(0,0);
   mass_matrix_Fu(0,1) = -0.7071067811865475*v*Yu(1,0);
   mass_matrix_Fu(0,2) = -0.7071067811865475*v*Yu(2,0);
   mass_matrix_Fu(1,0) = -0.7071067811865475*v*Yu(0,1);
   mass_matrix_Fu(1,1) = -0.7071067811865475*v*Yu(1,1);
   mass_matrix_Fu(1,2) = -0.7071067811865475*v*Yu(2,1);
   mass_matrix_Fu(2,0) = -0.7071067811865475*v*Yu(0,2);
   mass_matrix_Fu(2,1) = -0.7071067811865475*v*Yu(1,2);
   mass_matrix_Fu(2,2) = -0.7071067811865475*v*Yu(2,2);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu, eigenvalue_error);
   problems.flag_bad_mass(MDM_info::Fu, eigenvalue_error > precision *
      Abs(MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fe() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = 0.7071067811865475*v*Ye(0,0);
   mass_matrix_Fe(0,1) = 0.7071067811865475*v*Ye(1,0);
   mass_matrix_Fe(0,2) = 0.7071067811865475*v*Ye(2,0);
   mass_matrix_Fe(1,0) = 0.7071067811865475*v*Ye(0,1);
   mass_matrix_Fe(1,1) = 0.7071067811865475*v*Ye(1,1);
   mass_matrix_Fe(1,2) = 0.7071067811865475*v*Ye(2,1);
   mass_matrix_Fe(2,0) = 0.7071067811865475*v*Ye(0,2);
   mass_matrix_Fe(2,1) = 0.7071067811865475*v*Ye(1,2);
   mass_matrix_Fe(2,2) = 0.7071067811865475*v*Ye(2,2);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue, eigenvalue_error);
   problems.flag_bad_mass(MDM_info::Fe, eigenvalue_error > precision *
      Abs(MFe(0)));
#else
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue);
#endif

}

double CLASSNAME::get_mass_matrix_VWp() const
{
   const double mass_matrix_VWp = Re(0.25*Sqr(g2)*Sqr(v));

   return mass_matrix_VWp;
}

void CLASSNAME::calculate_MVWp()
{
   const auto mass_matrix_VWp = get_mass_matrix_VWp();
   MVWp = mass_matrix_VWp;

   if (MVWp < 0.) {
      problems.flag_tachyon(MDM_info::VWp);
   }

   MVWp = AbsSqrt(MVWp);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_VPVZ() const
{
   Eigen::Matrix<double,2,2> mass_matrix_VPVZ;

   mass_matrix_VPVZ(0,0) = 0.15*Sqr(g1)*Sqr(v);
   mass_matrix_VPVZ(0,1) = -0.19364916731037085*g1*g2*Sqr(v);
   mass_matrix_VPVZ(1,1) = 0.25*Sqr(g2)*Sqr(v);

   Symmetrize(mass_matrix_VPVZ);

   return mass_matrix_VPVZ;
}

void CLASSNAME::calculate_MVPVZ()
{
   const auto mass_matrix_VPVZ(get_mass_matrix_VPVZ());
   Eigen::Array<double,2,1> MVPVZ;


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ, eigenvalue_error
      );
   ZZ.transposeInPlace();
#else
   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ);
   ZZ.transposeInPlace();
#endif


   MVPVZ = AbsSqrt(MVPVZ);

   MVP = 0.;
   MVZ = MVPVZ(1);
}


double CLASSNAME::get_ewsb_eq_hh_1() const
{
   double result = Re(-(mu2*v) + LamH*Power(v,3));

   return result;
}



double CLASSNAME::CpconjHpHphh() const
{
   double result = 0.0;

   result = -2*LamH*v;

   return result;
}

double CLASSNAME::CpconjHpVWpVP() const
{
   double result = 0.0;

   result = 0.3872983346207417*g1*g2*v*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjHpVZVWp() const
{
   double result = 0.0;

   result = -0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpHpgWpCbargZ() const
{
   double result = 0.0;

   result = 0.25*g2*v*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjHpbargWpCgZ() const
{
   double result = 0.0;

   result = 0.05*g2*v*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW())
      );

   return result;
}

double CLASSNAME::CpHpgZbargWp() const
{
   double result = 0.0;

   result = 0.05*g2*v*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW())
      );

   return result;
}

double CLASSNAME::CpconjHpbargZgWp() const
{
   double result = 0.0;

   result = 0.25*g2*v*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpHpconjHpAhAh() const
{
   double result = 0.0;

   result = -2*LamH;

   return result;
}

double CLASSNAME::CpHpconjHphhhh() const
{
   double result = 0.0;

   result = -2*LamH;

   return result;
}

double CLASSNAME::CpHpconjHpconjHpHp() const
{
   double result = 0.0;

   result = -4*LamH;

   return result;
}

std::complex<double> CLASSNAME::CpconjHpVWpAh() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2;

   return result;
}

double CLASSNAME::CpconjHpVWphh() const
{
   double result = 0.0;

   result = -0.5*g2;

   return result;
}

double CLASSNAME::CpconjHpVPHp() const
{
   double result = 0.0;

   result = 0.1*(-3.872983346207417*g1*Cos(ThetaW()) - 5*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjHpVZHp() const
{
   double result = 0.0;

   result = 0.1*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpHpconjHpconjVWpVWp() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpHpconjHpVZVZ() const
{
   std::complex<double> result;

   result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(
      g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpconjHpbarFdFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_0;
   std::complex<double> tmp_1;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2;
      std::complex<double> tmp_3;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_2 += tmp_3;
      tmp_1 += (Vd(gI1,j2)) * tmp_2;
   }
   tmp_0 += tmp_1;
   result += (-1) * tmp_0;

   return result;
}

std::complex<double> CLASSNAME::CpconjHpbarFdFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4;
   std::complex<double> tmp_5;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_6;
      std::complex<double> tmp_7;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_7 += Conj(Ud(gI1,j1))*Yd(j1,j2);
      }
      tmp_6 += tmp_7;
      tmp_5 += (Conj(Vu(gI2,j2))) * tmp_6;
   }
   tmp_4 += tmp_5;
   result += (-1) * tmp_4;

   return result;
}

double CLASSNAME::CpconjHpbarFeFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjHpbarFeFvPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_8;
   std::complex<double> tmp_9;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_9 += Conj(Ue(gI1,j1))*Ye(j1,gI2);
   }
   tmp_8 += tmp_9;
   result += (-1) * tmp_8;

   return result;
}

double CLASSNAME::CpAhhhAh() const
{
   double result = 0.0;

   result = -2*LamH*v;

   return result;
}

std::complex<double> CLASSNAME::CpAhbargWpgWp() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*v*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpAhbargWpCgWpC() const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*v*Sqr(g2);

   return result;
}

double CLASSNAME::CpAhAhAhAh() const
{
   double result = 0.0;

   result = -6*LamH;

   return result;
}

double CLASSNAME::CpAhAhhhhh() const
{
   double result = 0.0;

   result = -2*LamH;

   return result;
}

double CLASSNAME::CpAhAhconjHpHp() const
{
   double result = 0.0;

   result = -2*LamH;

   return result;
}

std::complex<double> CLASSNAME::CpAhVZhh() const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(5*g2*Cos(ThetaW()) + 3.872983346207417
      *g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpAhconjVWpHp() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2;

   return result;
}

double CLASSNAME::CpAhAhconjVWpVWp() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhVZVZ() const
{
   std::complex<double> result;

   result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*
      Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_10;
   std::complex<double> tmp_11;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_12;
      std::complex<double> tmp_13;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_13 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_12 += tmp_13;
      tmp_11 += (Vd(gI1,j2)) * tmp_12;
   }
   tmp_10 += tmp_11;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_10;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_14;
   std::complex<double> tmp_15;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_16;
      std::complex<double> tmp_17;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_17 += Conj(Ud(gI1,j1))*Yd(j1,j2);
      }
      tmp_16 += tmp_17;
      tmp_15 += (Conj(Vd(gI2,j2))) * tmp_16;
   }
   tmp_14 += tmp_15;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_14;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFeFePR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_18;
   std::complex<double> tmp_19;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_20;
      std::complex<double> tmp_21;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_21 += Conj(Ye(j1,j2))*Ue(gI2,j1);
      }
      tmp_20 += tmp_21;
      tmp_19 += (Ve(gI1,j2)) * tmp_20;
   }
   tmp_18 += tmp_19;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_18;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFeFePL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_22;
   std::complex<double> tmp_23;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_24;
      std::complex<double> tmp_25;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_25 += Conj(Ue(gI1,j1))*Ye(j1,j2);
      }
      tmp_24 += tmp_25;
      tmp_23 += (Conj(Ve(gI2,j2))) * tmp_24;
   }
   tmp_22 += tmp_23;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_22;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_26;
   std::complex<double> tmp_27;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_28;
      std::complex<double> tmp_29;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_29 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_28 += tmp_29;
      tmp_27 += (Vu(gI1,j2)) * tmp_28;
   }
   tmp_26 += tmp_27;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_26;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_30;
   std::complex<double> tmp_31;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_32;
      std::complex<double> tmp_33;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_33 += Conj(Uu(gI1,j1))*Yu(j1,j2);
      }
      tmp_32 += tmp_33;
      tmp_31 += (Conj(Vu(gI2,j2))) * tmp_32;
   }
   tmp_30 += tmp_31;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_30;

   return result;
}

double CLASSNAME::CphhAhAh() const
{
   double result = 0.0;

   result = -2*LamH*v;

   return result;
}

double CLASSNAME::Cphhhhhh() const
{
   double result = 0.0;

   result = -6*LamH*v;

   return result;
}

double CLASSNAME::CphhVZVZ() const
{
   double result = 0.0;

   result = 0.5*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CphhbargWpgWp() const
{
   double result = 0.0;

   result = -0.25*v*Sqr(g2);

   return result;
}

double CLASSNAME::CphhbargWpCgWpC() const
{
   double result = 0.0;

   result = -0.25*v*Sqr(g2);

   return result;
}

double CLASSNAME::CphhbargZgZ() const
{
   double result = 0.0;

   result = -0.25*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))
      ;

   return result;
}

double CLASSNAME::CphhconjHpHp() const
{
   double result = 0.0;

   result = -2*LamH*v;

   return result;
}

double CLASSNAME::CphhconjVWpVWp() const
{
   double result = 0.0;

   result = 0.5*v*Sqr(g2);

   return result;
}

double CLASSNAME::CphhhhAhAh() const
{
   double result = 0.0;

   result = -2*LamH;

   return result;
}

double CLASSNAME::Cphhhhhhhh() const
{
   double result = 0.0;

   result = -6*LamH;

   return result;
}

double CLASSNAME::CphhhhconjHpHp() const
{
   double result = 0.0;

   result = -2*LamH;

   return result;
}

std::complex<double> CLASSNAME::CphhVZAh() const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(5*g2*Cos(ThetaW()) + 3.872983346207417
      *g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CphhconjVWpHp() const
{
   double result = 0.0;

   result = 0.5*g2;

   return result;
}

double CLASSNAME::CphhhhconjVWpVWp() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CphhhhVZVZ() const
{
   std::complex<double> result;

   result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*
      Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CphhbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_34;
   std::complex<double> tmp_35;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_36;
      std::complex<double> tmp_37;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_37 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_36 += tmp_37;
      tmp_35 += (Vd(gI1,j2)) * tmp_36;
   }
   tmp_34 += tmp_35;
   result += (-0.7071067811865475) * tmp_34;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_38;
   std::complex<double> tmp_39;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_40;
      std::complex<double> tmp_41;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_41 += Conj(Ud(gI1,j1))*Yd(j1,j2);
      }
      tmp_40 += tmp_41;
      tmp_39 += (Conj(Vd(gI2,j2))) * tmp_40;
   }
   tmp_38 += tmp_39;
   result += (-0.7071067811865475) * tmp_38;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFeFePR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_42;
   std::complex<double> tmp_43;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_44;
      std::complex<double> tmp_45;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_45 += Conj(Ye(j1,j2))*Ue(gI2,j1);
      }
      tmp_44 += tmp_45;
      tmp_43 += (Ve(gI1,j2)) * tmp_44;
   }
   tmp_42 += tmp_43;
   result += (-0.7071067811865475) * tmp_42;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFeFePL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_46;
   std::complex<double> tmp_47;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_48;
      std::complex<double> tmp_49;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_49 += Conj(Ue(gI1,j1))*Ye(j1,j2);
      }
      tmp_48 += tmp_49;
      tmp_47 += (Conj(Ve(gI2,j2))) * tmp_48;
   }
   tmp_46 += tmp_47;
   result += (-0.7071067811865475) * tmp_46;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_50;
   std::complex<double> tmp_51;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_52;
      std::complex<double> tmp_53;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_53 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_52 += tmp_53;
      tmp_51 += (Vu(gI1,j2)) * tmp_52;
   }
   tmp_50 += tmp_51;
   result += (0.7071067811865475) * tmp_50;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_54;
   std::complex<double> tmp_55;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_56;
      std::complex<double> tmp_57;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_57 += Conj(Uu(gI1,j1))*Yu(j1,j2);
      }
      tmp_56 += tmp_57;
      tmp_55 += (Conj(Vu(gI2,j2))) * tmp_56;
   }
   tmp_54 += tmp_55;
   result += (0.7071067811865475) * tmp_54;

   return result;
}

std::complex<double> CLASSNAME::CpVGVGVG() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpVGbargGgG() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3;

   return result;
}

double CLASSNAME::CpVGbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpVGbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpVGbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpVGbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpVGVGVGVG1() const
{
   double result = 0.0;

   result = -16*Sqr(g3);

   return result;
}

double CLASSNAME::CpVGVGVGVG2() const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpVGVGVGVG3() const
{
   double result = 0.0;

   result = 16*Sqr(g3);

   return result;
}

double CLASSNAME::CpVPbargWpgWp() const
{
   double result = 0.0;

   result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVPbargWpCgWpC() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPconjHpHp() const
{
   double result = 0.0;

   result = 0.1*(-3.872983346207417*g1*Cos(ThetaW()) - 5*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPconjVWpHp() const
{
   double result = 0.0;

   result = 0.3872983346207417*g1*g2*v*Cos(ThetaW());

   return result;
}

std::complex<double> CLASSNAME::CpVPVPconjHpHp() const
{
   std::complex<double> result;

   result = 0.1*(g2*Sin(ThetaW())*(7.745966692414834*g1*Cos(ThetaW()) + 5*g2*
      Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos(ThetaW())));

   return result;
}

double CLASSNAME::CpVPconjVWpVWp() const
{
   double result = 0.0;

   result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVPbarFcFcPL() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPbarFcFcPR() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPbarFgFgPL() const
{
   double result = 0.0;

   result = -2*g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVPbarFgFgPR() const
{
   double result = 0.0;

   result = -2*g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVPbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1
      *Cos(ThetaW()) - 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.2581988897471611*g1*Cos(ThetaW())*KroneckerDelta(gI1,gI2);

   return result;
}

double CLASSNAME::CpVPbarFeFePL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(ThetaW()) +
      g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPbarFeFePR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.7745966692414834*g1*Cos(ThetaW())*KroneckerDelta(gI1,gI2);

   return result;
}

double CLASSNAME::CpVPbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1
      *Cos(ThetaW()) + 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5163977794943222*g1*Cos(ThetaW())*KroneckerDelta(gI1,gI2);

   return result;
}

double CLASSNAME::CpVPVPconjVWpVWp1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPVPconjVWpVWp2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPVPconjVWpVWp3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVZhhAh() const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(5*g2*Cos(ThetaW()) + 3.872983346207417
      *g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZVZhh() const
{
   double result = 0.0;

   result = 0.5*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbargWpgWp() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpVZbargWpCgWpC() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZconjHpHp() const
{
   double result = 0.0;

   result = 0.1*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZconjVWpHp() const
{
   double result = 0.0;

   result = -0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

std::complex<double> CLASSNAME::CpVZVZAhAh() const
{
   std::complex<double> result;

   result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*
      Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZhhhh() const
{
   std::complex<double> result;

   result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*
      Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjHpHp() const
{
   std::complex<double> result;

   result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(
      g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())));

   return result;
}

double CLASSNAME::CpVZconjVWpVWp() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFcFcPL() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFcFcPR() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFgFgPL() const
{
   double result = 0.0;

   result = -2*g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFgFgPR() const
{
   double result = 0.0;

   result = -2*g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.16666666666666666*KroneckerDelta(gI1,gI2)*(3*g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.2581988897471611*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFeFePL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) - 0.7745966692414834*
      g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFeFePR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.7745966692414834*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.03333333333333333*KroneckerDelta(gI1,gI2)*(-15*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5163977794943222*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFvFvPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) + 0.7745966692414834
      *g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFvFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpVZVZconjVWpVWp1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZVZconjVWpVWp2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZVZconjVWpVWp3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWpHpAh() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2;

   return result;
}

double CLASSNAME::CpconjVWpHphh() const
{
   double result = 0.0;

   result = 0.5*g2;

   return result;
}

double CLASSNAME::CpconjVWpVPHp() const
{
   double result = 0.0;

   result = 0.3872983346207417*g1*g2*v*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpVWphh() const
{
   double result = 0.0;

   result = 0.5*v*Sqr(g2);

   return result;
}

double CLASSNAME::CpconjVWpVZHp() const
{
   double result = 0.0;

   result = -0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpbargPgWp() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWpbargWpCgP() const
{
   double result = 0.0;

   result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpbargWpCgZ() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpbargZgWp() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpAhAh() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CpVWpconjVWphhhh() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CpVWpconjVWpconjHpHp() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CpconjVWpVWpVP() const
{
   double result = 0.0;

   result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpVZVWp() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWpFnFcPL() const
{
   double result = 0.0;

   result = 1.7320508075688772*g2;

   return result;
}

double CLASSNAME::CpconjVWpFnFcPR() const
{
   double result = 0.0;

   result = 1.7320508075688772*g2;

   return result;
}

double CLASSNAME::CpconjVWpbarFcFgPL() const
{
   double result = 0.0;

   result = 1.4142135623730951*g2;

   return result;
}

double CLASSNAME::CpconjVWpbarFcFgPR() const
{
   double result = 0.0;

   result = 1.4142135623730951*g2;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWpbarFdFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_58;
   std::complex<double> tmp_59;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_59 += Conj(Vu(gI2,j1))*Vd(gI1,j1);
   }
   tmp_58 += tmp_59;
   result += (-0.7071067811865475*g2) * tmp_58;

   return result;
}

double CLASSNAME::CpconjVWpbarFdFuPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWpbarFeFvPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7071067811865475*g2*Ve(gI1,gI2);
   }

   return result;
}

double CLASSNAME::CpconjVWpbarFeFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpVWpconjVWpVPVP1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVPVP2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVPVP3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVZVZ1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVZVZ2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVZVZ3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpconjVWpVWp1() const
{
   double result = 0.0;

   result = -Sqr(g2);

   return result;
}

double CLASSNAME::CpVWpconjVWpconjVWpVWp2() const
{
   double result = 0.0;

   result = -Sqr(g2);

   return result;
}

double CLASSNAME::CpVWpconjVWpconjVWpVWp3() const
{
   double result = 0.0;

   result = 2*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_60;
      std::complex<double> tmp_61;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_61 += Conj(Vd(gI1,j2))*Yd(gO2,j2);
      }
      tmp_60 += tmp_61;
      result += (std::complex<double>(0.,0.7071067811865475)) * tmp_60;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_62;
      std::complex<double> tmp_63;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_63 += Conj(Yd(j1,gO1))*Ud(gI1,j1);
      }
      tmp_62 += tmp_63;
      result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_62;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_64;
      std::complex<double> tmp_65;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_65 += Conj(Vd(gI2,j2))*Yd(gO2,j2);
      }
      tmp_64 += tmp_65;
      result += (-0.7071067811865475) * tmp_64;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_66;
      std::complex<double> tmp_67;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_67 += Conj(Yd(j1,gO1))*Ud(gI2,j1);
      }
      tmp_66 += tmp_67;
      result += (-0.7071067811865475) * tmp_66;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVGFdPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Ud(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVGFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Conj(Vd(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVPFdPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.2581988897471611*g1*Cos(ThetaW())*Ud(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVPFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.12909944487358055*g1*Conj(Vd(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.5*g2*Conj(Vd(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVZFdPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.2581988897471611*g1*Sin(ThetaW())*Ud(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVZFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(Vd(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.12909944487358055*g1*Conj(Vd(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdconjHpFuPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_68;
      std::complex<double> tmp_69;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_69 += Conj(Vu(gI2,j2))*Yd(gO2,j2);
      }
      tmp_68 += tmp_69;
      result += (-1) * tmp_68;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdconjHpFuPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_70;
      std::complex<double> tmp_71;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_71 += Conj(Yu(j1,gO1))*Uu(gI2,j1);
      }
      tmp_70 += tmp_71;
      result += (-1) * tmp_70;
   }

   return result;
}

double CLASSNAME::CpbarUFdconjVWpFuPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdconjVWpFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(Vu(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_72;
      std::complex<double> tmp_73;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_73 += Conj(Vu(gI1,j2))*Yu(gO2,j2);
      }
      tmp_72 += tmp_73;
      result += (std::complex<double>(0.,0.7071067811865475)) * tmp_72;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_74;
      std::complex<double> tmp_75;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_75 += Conj(Yu(j1,gO1))*Uu(gI1,j1);
      }
      tmp_74 += tmp_75;
      result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_74;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_76;
      std::complex<double> tmp_77;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_77 += Conj(Vu(gI2,j2))*Yu(gO2,j2);
      }
      tmp_76 += tmp_77;
      result += (0.7071067811865475) * tmp_76;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_78;
      std::complex<double> tmp_79;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_79 += Conj(Yu(j1,gO1))*Uu(gI2,j1);
      }
      tmp_78 += tmp_79;
      result += (0.7071067811865475) * tmp_78;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuHpFdPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_80;
      std::complex<double> tmp_81;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_81 += Conj(Vd(gI2,j2))*Yu(gO2,j2);
      }
      tmp_80 += tmp_81;
      result += (-1) * tmp_80;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuHpFdPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_82;
      std::complex<double> tmp_83;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_83 += Conj(Yd(j1,gO1))*Ud(gI2,j1);
      }
      tmp_82 += tmp_83;
      result += (-1) * tmp_82;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVGFuPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Uu(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVGFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Conj(Vu(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVPFuPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.5163977794943222*g1*Cos(ThetaW())*Uu(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVPFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.12909944487358055*g1*Conj(Vu(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += -0.5*g2*Conj(Vu(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

double CLASSNAME::CpbarUFuVWpFdPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVWpFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(Vd(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVZFuPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5163977794943222*g1*Sin(ThetaW())*Uu(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVZFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.5*g2*Conj(Vu(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.12909944487358055*g1*Conj(Vu(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_84;
      std::complex<double> tmp_85;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_85 += Conj(Ve(gI1,j2))*Ye(gO2,j2);
      }
      tmp_84 += tmp_85;
      result += (std::complex<double>(0.,0.7071067811865475)) * tmp_84;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_86;
      std::complex<double> tmp_87;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_87 += Conj(Ye(j1,gO1))*Ue(gI1,j1);
      }
      tmp_86 += tmp_87;
      result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_86;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_88;
      std::complex<double> tmp_89;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_89 += Conj(Ve(gI2,j2))*Ye(gO2,j2);
      }
      tmp_88 += tmp_89;
      result += (-0.7071067811865475) * tmp_88;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_90;
      std::complex<double> tmp_91;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_91 += Conj(Ye(j1,gO1))*Ue(gI2,j1);
      }
      tmp_90 += tmp_91;
      result += (-0.7071067811865475) * tmp_90;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVPFePR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.7745966692414834*g1*Cos(ThetaW())*Ue(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVPFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.3872983346207417*g1*Conj(Ve(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.5*g2*Conj(Ve(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVZFePR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7745966692414834*g1*Sin(ThetaW())*Ue(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVZFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(Ve(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += -0.3872983346207417*g1*Conj(Ve(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeconjHpFvPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += -Ye(gO2,gI2);
   }

   return result;
}

double CLASSNAME::CpbarUFeconjHpFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpbarUFeconjVWpFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpbarUFeconjVWpFvPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   if (gI2 < 3) {
      result += -0.7071067811865475*g2*KroneckerDelta(gI2,gO1);
   }

   return result;
}

double CLASSNAME::CpbarFvHpFePL(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvHpFePR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_92;
   std::complex<double> tmp_93;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_93 += Conj(Ye(j1,gO1))*Ue(gI2,j1);
   }
   tmp_92 += tmp_93;
   result += (-1) * tmp_92;

   return result;
}

double CLASSNAME::CpbarFvVWpFePR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvVWpFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(Ve(gI2,gO1));
   }

   return result;
}

double CLASSNAME::CpbarFvVZFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpbarFvVZFvPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5*KroneckerDelta(gI2,gO1)*(g2*Cos(ThetaW()) + 0.7745966692414834
      *g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFcVPFcPR() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFcVPFcPL() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFcVWpFnPR() const
{
   double result = 0.0;

   result = 1.7320508075688772*g2;

   return result;
}

double CLASSNAME::CpbarFcVWpFnPL() const
{
   double result = 0.0;

   result = 1.7320508075688772*g2;

   return result;
}

double CLASSNAME::CpbarFcVZFcPR() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFcVZFcPL() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFcconjVWpFgPR() const
{
   double result = 0.0;

   result = 1.4142135623730951*g2;

   return result;
}

double CLASSNAME::CpbarFcconjVWpFgPL() const
{
   double result = 0.0;

   result = 1.4142135623730951*g2;

   return result;
}

double CLASSNAME::CpbarFgVPFgPR() const
{
   double result = 0.0;

   result = -2*g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarFgVPFgPL() const
{
   double result = 0.0;

   result = -2*g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarFgVWpFcPR() const
{
   double result = 0.0;

   result = 1.4142135623730951*g2;

   return result;
}

double CLASSNAME::CpbarFgVWpFcPL() const
{
   double result = 0.0;

   result = 1.4142135623730951*g2;

   return result;
}

double CLASSNAME::CpbarFgVZFgPR() const
{
   double result = 0.0;

   result = -2*g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpbarFgVZFgPL() const
{
   double result = 0.0;

   result = -2*g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpFnbarFcVWpPR() const
{
   double result = 0.0;

   result = -1.7320508075688772*g2;

   return result;
}

double CLASSNAME::CpFnbarFcVWpPL() const
{
   double result = 0.0;

   result = -1.7320508075688772*g2;

   return result;
}

double CLASSNAME::CpFnconjVWpFcPR() const
{
   double result = 0.0;

   result = 1.7320508075688772*g2;

   return result;
}

double CLASSNAME::CpFnconjVWpFcPL() const
{
   double result = 0.0;

   result = 1.7320508075688772*g2;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_94;
   std::complex<double> tmp_95;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_96;
      std::complex<double> tmp_97;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_97 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_96 += tmp_97;
      tmp_95 += (Conj(Vd(gI1,j2))) * tmp_96;
   }
   tmp_94 += tmp_95;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_94;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_98;
   std::complex<double> tmp_99;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_100;
      std::complex<double> tmp_101;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_101 += Conj(Yd(j1,j2))*Ud(gI1,j1);
      }
      tmp_100 += tmp_101;
      tmp_99 += (Vd(gO1,j2)) * tmp_100;
   }
   tmp_98 += tmp_99;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_98;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_102;
   std::complex<double> tmp_103;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_104;
      std::complex<double> tmp_105;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_105 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_104 += tmp_105;
      tmp_103 += (Conj(Vd(gI2,j2))) * tmp_104;
   }
   tmp_102 += tmp_103;
   result += (-0.7071067811865475) * tmp_102;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_106;
   std::complex<double> tmp_107;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_108;
      std::complex<double> tmp_109;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_109 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_108 += tmp_109;
      tmp_107 += (Vd(gO1,j2)) * tmp_108;
   }
   tmp_106 += tmp_107;
   result += (-0.7071067811865475) * tmp_106;

   return result;
}

double CLASSNAME::CpbarFdVZFdPR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = -0.2581988897471611*g1*KroneckerDelta(gI2,gO2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarFdVZFdPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.16666666666666666*KroneckerDelta(gI2,gO1)*(3*g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdconjHpFuPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_110;
   std::complex<double> tmp_111;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_112;
      std::complex<double> tmp_113;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_113 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_112 += tmp_113;
      tmp_111 += (Conj(Vu(gI2,j2))) * tmp_112;
   }
   tmp_110 += tmp_111;
   result += (-1) * tmp_110;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdconjHpFuPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_114;
   std::complex<double> tmp_115;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_116;
      std::complex<double> tmp_117;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_117 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_116 += tmp_117;
      tmp_115 += (Vd(gO1,j2)) * tmp_116;
   }
   tmp_114 += tmp_115;
   result += (-1) * tmp_114;

   return result;
}

double CLASSNAME::CpbarFdconjVWpFuPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdconjVWpFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_118;
   std::complex<double> tmp_119;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_119 += Conj(Vu(gI2,j1))*Vd(gO1,j1);
   }
   tmp_118 += tmp_119;
   result += (-0.7071067811865475*g2) * tmp_118;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_120;
   std::complex<double> tmp_121;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_122;
      std::complex<double> tmp_123;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_123 += Conj(Ue(gO2,j1))*Ye(j1,j2);
      }
      tmp_122 += tmp_123;
      tmp_121 += (Conj(Ve(gI1,j2))) * tmp_122;
   }
   tmp_120 += tmp_121;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_120;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_124;
   std::complex<double> tmp_125;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_126;
      std::complex<double> tmp_127;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_127 += Conj(Ye(j1,j2))*Ue(gI1,j1);
      }
      tmp_126 += tmp_127;
      tmp_125 += (Ve(gO1,j2)) * tmp_126;
   }
   tmp_124 += tmp_125;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_124;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_128;
   std::complex<double> tmp_129;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_130;
      std::complex<double> tmp_131;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_131 += Conj(Ue(gO2,j1))*Ye(j1,j2);
      }
      tmp_130 += tmp_131;
      tmp_129 += (Conj(Ve(gI2,j2))) * tmp_130;
   }
   tmp_128 += tmp_129;
   result += (-0.7071067811865475) * tmp_128;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_132;
   std::complex<double> tmp_133;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_134;
      std::complex<double> tmp_135;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_135 += Conj(Ye(j1,j2))*Ue(gI2,j1);
      }
      tmp_134 += tmp_135;
      tmp_133 += (Ve(gO1,j2)) * tmp_134;
   }
   tmp_132 += tmp_133;
   result += (-0.7071067811865475) * tmp_132;

   return result;
}

double CLASSNAME::CpbarFeVZFePR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = -0.7745966692414834*g1*KroneckerDelta(gI2,gO2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarFeVZFePL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI2,gO1)*(g2*Cos(ThetaW()) - 0.7745966692414834*
      g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeconjHpFvPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_136;
   std::complex<double> tmp_137;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_137 += Conj(Ue(gO2,j1))*Ye(j1,gI2);
   }
   tmp_136 += tmp_137;
   result += (-1) * tmp_136;

   return result;
}

double CLASSNAME::CpbarFeconjHpFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpbarFeconjVWpFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeconjVWpFvPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7071067811865475*g2*Ve(gO1,gI2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_138;
   std::complex<double> tmp_139;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_140;
      std::complex<double> tmp_141;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_141 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_140 += tmp_141;
      tmp_139 += (Conj(Vu(gI1,j2))) * tmp_140;
   }
   tmp_138 += tmp_139;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_138;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_142;
   std::complex<double> tmp_143;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_144;
      std::complex<double> tmp_145;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_145 += Conj(Yu(j1,j2))*Uu(gI1,j1);
      }
      tmp_144 += tmp_145;
      tmp_143 += (Vu(gO1,j2)) * tmp_144;
   }
   tmp_142 += tmp_143;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_142;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_146;
   std::complex<double> tmp_147;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_148;
      std::complex<double> tmp_149;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_149 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_148 += tmp_149;
      tmp_147 += (Conj(Vu(gI2,j2))) * tmp_148;
   }
   tmp_146 += tmp_147;
   result += (0.7071067811865475) * tmp_146;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_150;
   std::complex<double> tmp_151;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_152;
      std::complex<double> tmp_153;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_153 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_152 += tmp_153;
      tmp_151 += (Vu(gO1,j2)) * tmp_152;
   }
   tmp_150 += tmp_151;
   result += (0.7071067811865475) * tmp_150;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuHpFdPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_154;
   std::complex<double> tmp_155;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_156;
      std::complex<double> tmp_157;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_157 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_156 += tmp_157;
      tmp_155 += (Conj(Vd(gI2,j2))) * tmp_156;
   }
   tmp_154 += tmp_155;
   result += (-1) * tmp_154;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuHpFdPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_158;
   std::complex<double> tmp_159;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_160;
      std::complex<double> tmp_161;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_161 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_160 += tmp_161;
      tmp_159 += (Vu(gO1,j2)) * tmp_160;
   }
   tmp_158 += tmp_159;
   result += (-1) * tmp_158;

   return result;
}

double CLASSNAME::CpbarFuVPFuPR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5163977794943222*g1*Cos(ThetaW())*KroneckerDelta(gI2,gO2);

   return result;
}

double CLASSNAME::CpbarFuVPFuPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI2,gO1)*(0.7745966692414834*g1
      *Cos(ThetaW()) + 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFuVWpFdPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuVWpFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_162;
   std::complex<double> tmp_163;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_163 += Conj(Vd(gI2,j1))*Vu(gO1,j1);
   }
   tmp_162 += tmp_163;
   result += (-0.7071067811865475*g2) * tmp_162;

   return result;
}

double CLASSNAME::CpbarFuVZFuPR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5163977794943222*g1*KroneckerDelta(gI2,gO2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarFuVZFuPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.03333333333333333*KroneckerDelta(gI2,gO1)*(-15*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}


std::complex<double> CLASSNAME::self_energy_Hp(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpconjHpHphh())*B0(p,MHp,Mhh);
   result += 4*AbsSqr(CpconjHpVWpVP())*(-0.5 + B0(p,0,MVWp));
   result += 4*AbsSqr(CpconjHpVZVWp())*(-0.5 + B0(p,MVWp,MVZ));
   result += -0.5*A0(MAh)*CpHpconjHpAhAh();
   result += -(A0(MHp)*CpHpconjHpconjHpHp());
   result += -0.5*A0(Mhh)*CpHpconjHphhhh();
   result += -(B0(p,MVZ,MVWp)*CpconjHpbargWpCgZ()*CpHpgWpCbargZ());
   result += -(B0(p,MVWp,MVZ)*CpconjHpbargZgWp()*CpHpgZbargWp());
   result += AbsSqr(CpconjHpVWpAh())*F0(p,MAh,MVWp);
   result += AbsSqr(CpconjHpVWphh())*F0(p,Mhh,MVWp);
   result += AbsSqr(CpconjHpVPHp())*F0(p,MHp,0);
   result += AbsSqr(CpconjHpVZHp())*F0(p,MHp,MVZ);
   std::complex<double> tmp_164;
   std::complex<double> tmp_165;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_166;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_166 += (AbsSqr(CpconjHpbarFdFuPL(gI1,gI2)) + AbsSqr(
            CpconjHpbarFdFuPR(gI1,gI2)))*G0(p,MFd(gI1),MFu(gI2));
      }
      tmp_165 += tmp_166;
   }
   tmp_164 += tmp_165;
   result += (3) * tmp_164;
   std::complex<double> tmp_167;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_168;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_168 += (AbsSqr(CpconjHpbarFeFvPL(gI1,gI2)) + AbsSqr(
            CpconjHpbarFeFvPR(gI1,gI2)))*G0(p,MFe(gI1),MFv(gI2));
      }
      tmp_167 += tmp_168;
   }
   result += tmp_167;
   std::complex<double> tmp_169;
   std::complex<double> tmp_170;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_171;
      std::complex<double> tmp_172;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_172 += B0(p,MFd(gI1),MFu(gI2))*(Conj(CpconjHpbarFdFuPR(gI1,
            gI2))*CpconjHpbarFdFuPL(gI1,gI2) + Conj(CpconjHpbarFdFuPL(gI1,gI2))*
            CpconjHpbarFdFuPR(gI1,gI2))*MFu(gI2);
      }
      tmp_171 += tmp_172;
      tmp_170 += (MFd(gI1)) * tmp_171;
   }
   tmp_169 += tmp_170;
   result += (-6) * tmp_169;
   std::complex<double> tmp_173;
   std::complex<double> tmp_174;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_175;
      std::complex<double> tmp_176;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_176 += B0(p,MFe(gI1),MFv(gI2))*(Conj(CpconjHpbarFeFvPR(gI1,
            gI2))*CpconjHpbarFeFvPL(gI1,gI2) + Conj(CpconjHpbarFeFvPL(gI1,gI2))*
            CpconjHpbarFeFvPR(gI1,gI2))*MFv(gI2);
      }
      tmp_175 += tmp_176;
      tmp_174 += (MFe(gI1)) * tmp_175;
   }
   tmp_173 += tmp_174;
   result += (-2) * tmp_173;
   result += 4*CpHpconjHpconjVWpVWp()*(A0(MVWp) - 0.5*Sqr(MVWp));
   result += 2*CpHpconjHpVZVZ()*(A0(MVZ) - 0.5*Sqr(MVZ));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Ah(double p ) const
{
   std::complex<double> result;

   result += -0.5*A0(MAh)*CpAhAhAhAh();
   result += -(A0(MHp)*CpAhAhconjHpHp());
   result += -0.5*A0(Mhh)*CpAhAhhhhh();
   result += -(B0(p,MVWp,MVWp)*Sqr(CpAhbargWpCgWpC()));
   result += -(B0(p,MVWp,MVWp)*Sqr(CpAhbargWpgWp()));
   result += AbsSqr(CpAhhhAh())*B0(p,Mhh,MAh);
   result += AbsSqr(CpAhVZhh())*F0(p,Mhh,MVZ);
   result += 2*AbsSqr(CpAhconjVWpHp())*F0(p,MHp,MVWp);
   std::complex<double> tmp_177;
   std::complex<double> tmp_178;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_179;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_179 += (AbsSqr(CpAhbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpAhbarFdFdPR(gI1,gI2)))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_178 += tmp_179;
   }
   tmp_177 += tmp_178;
   result += (3) * tmp_177;
   std::complex<double> tmp_180;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_181;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_181 += (AbsSqr(CpAhbarFeFePL(gI1,gI2)) + AbsSqr(
            CpAhbarFeFePR(gI1,gI2)))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_180 += tmp_181;
   }
   result += tmp_180;
   std::complex<double> tmp_182;
   std::complex<double> tmp_183;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_184;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_184 += (AbsSqr(CpAhbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpAhbarFuFuPR(gI1,gI2)))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_183 += tmp_184;
   }
   tmp_182 += tmp_183;
   result += (3) * tmp_182;
   std::complex<double> tmp_185;
   std::complex<double> tmp_186;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_187;
      std::complex<double> tmp_188;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_188 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpAhbarFdFdPR(gI1,gI2))
            *CpAhbarFdFdPL(gI1,gI2) + Conj(CpAhbarFdFdPL(gI1,gI2))*CpAhbarFdFdPR(
            gI1,gI2))*MFd(gI2);
      }
      tmp_187 += tmp_188;
      tmp_186 += (MFd(gI1)) * tmp_187;
   }
   tmp_185 += tmp_186;
   result += (-6) * tmp_185;
   std::complex<double> tmp_189;
   std::complex<double> tmp_190;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_191;
      std::complex<double> tmp_192;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_192 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpAhbarFeFePR(gI1,gI2))
            *CpAhbarFeFePL(gI1,gI2) + Conj(CpAhbarFeFePL(gI1,gI2))*CpAhbarFeFePR(
            gI1,gI2))*MFe(gI2);
      }
      tmp_191 += tmp_192;
      tmp_190 += (MFe(gI1)) * tmp_191;
   }
   tmp_189 += tmp_190;
   result += (-2) * tmp_189;
   std::complex<double> tmp_193;
   std::complex<double> tmp_194;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_195;
      std::complex<double> tmp_196;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_196 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpAhbarFuFuPR(gI1,gI2))
            *CpAhbarFuFuPL(gI1,gI2) + Conj(CpAhbarFuFuPL(gI1,gI2))*CpAhbarFuFuPR(
            gI1,gI2))*MFu(gI2);
      }
      tmp_195 += tmp_196;
      tmp_194 += (MFu(gI1)) * tmp_195;
   }
   tmp_193 += tmp_194;
   result += (-6) * tmp_193;
   result += 4*CpAhAhconjVWpVWp()*(A0(MVWp) - 0.5*Sqr(MVWp));
   result += 2*CpAhAhVZVZ()*(A0(MVZ) - 0.5*Sqr(MVZ));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_hh(double p ) const
{
   std::complex<double> result;

   result += 0.5*AbsSqr(CphhAhAh())*B0(p,MAh,MAh);
   result += -(B0(p,MVWp,MVWp)*Sqr(CphhbargWpCgWpC()));
   result += -(B0(p,MVWp,MVWp)*Sqr(CphhbargWpgWp()));
   result += -(B0(p,MVZ,MVZ)*Sqr(CphhbargZgZ()));
   result += AbsSqr(CphhconjHpHp())*B0(p,MHp,MHp);
   result += 4*AbsSqr(CphhconjVWpVWp())*(-0.5 + B0(p,MVWp,MVWp));
   result += -0.5*A0(MAh)*CphhhhAhAh();
   result += -(A0(MHp)*CphhhhconjHpHp());
   result += 0.5*AbsSqr(Cphhhhhh())*B0(p,Mhh,Mhh);
   result += -0.5*A0(Mhh)*Cphhhhhhhh();
   result += 2*AbsSqr(CphhVZVZ())*(-0.5 + B0(p,MVZ,MVZ));
   result += AbsSqr(CphhVZAh())*F0(p,MAh,MVZ);
   result += 2*AbsSqr(CphhconjVWpHp())*F0(p,MHp,MVWp);
   std::complex<double> tmp_197;
   std::complex<double> tmp_198;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_199;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_199 += (AbsSqr(CphhbarFdFdPL(gI1,gI2)) + AbsSqr(
            CphhbarFdFdPR(gI1,gI2)))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_198 += tmp_199;
   }
   tmp_197 += tmp_198;
   result += (3) * tmp_197;
   std::complex<double> tmp_200;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_201;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_201 += (AbsSqr(CphhbarFeFePL(gI1,gI2)) + AbsSqr(
            CphhbarFeFePR(gI1,gI2)))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_200 += tmp_201;
   }
   result += tmp_200;
   std::complex<double> tmp_202;
   std::complex<double> tmp_203;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_204;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_204 += (AbsSqr(CphhbarFuFuPL(gI1,gI2)) + AbsSqr(
            CphhbarFuFuPR(gI1,gI2)))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_203 += tmp_204;
   }
   tmp_202 += tmp_203;
   result += (3) * tmp_202;
   std::complex<double> tmp_205;
   std::complex<double> tmp_206;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_207;
      std::complex<double> tmp_208;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_208 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CphhbarFdFdPR(gI1,gI2))
            *CphhbarFdFdPL(gI1,gI2) + Conj(CphhbarFdFdPL(gI1,gI2))*CphhbarFdFdPR(
            gI1,gI2))*MFd(gI2);
      }
      tmp_207 += tmp_208;
      tmp_206 += (MFd(gI1)) * tmp_207;
   }
   tmp_205 += tmp_206;
   result += (-6) * tmp_205;
   std::complex<double> tmp_209;
   std::complex<double> tmp_210;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_211;
      std::complex<double> tmp_212;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_212 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CphhbarFeFePR(gI1,gI2))
            *CphhbarFeFePL(gI1,gI2) + Conj(CphhbarFeFePL(gI1,gI2))*CphhbarFeFePR(
            gI1,gI2))*MFe(gI2);
      }
      tmp_211 += tmp_212;
      tmp_210 += (MFe(gI1)) * tmp_211;
   }
   tmp_209 += tmp_210;
   result += (-2) * tmp_209;
   std::complex<double> tmp_213;
   std::complex<double> tmp_214;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_215;
      std::complex<double> tmp_216;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_216 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CphhbarFuFuPR(gI1,gI2))
            *CphhbarFuFuPL(gI1,gI2) + Conj(CphhbarFuFuPL(gI1,gI2))*CphhbarFuFuPR(
            gI1,gI2))*MFu(gI2);
      }
      tmp_215 += tmp_216;
      tmp_214 += (MFu(gI1)) * tmp_215;
   }
   tmp_213 += tmp_214;
   result += (-6) * tmp_213;
   result += 4*CphhhhconjVWpVWp()*(A0(MVWp) - 0.5*Sqr(MVWp));
   result += 2*CphhhhVZVZ()*(A0(MVZ) - 0.5*Sqr(MVZ));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VG(double p ) const
{
   std::complex<double> result;

   result += 3*AbsSqr(CpVGbargGgG())*B00(p,MVG,MVG);
   result += -1.5*AbsSqr(CpVGVGVG())*(10*B00(p,0,0) + 0.6666666666666666*Sqr(p)
      + 4*B0(p,0,0)*Sqr(p));
   result += 0;
   std::complex<double> tmp_217;
   std::complex<double> tmp_218;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_219;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_219 += (AbsSqr(CpVGbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVGbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_219 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVGbarFdFdPL(gI1,gI2))*CpVGbarFdFdPR(gI1,gI2));
      }
      tmp_218 += tmp_219;
   }
   tmp_217 += tmp_218;
   result += (0.5) * tmp_217;
   std::complex<double> tmp_220;
   std::complex<double> tmp_221;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_222;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_222 += (AbsSqr(CpVGbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVGbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_222 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVGbarFuFuPL(gI1,gI2))*CpVGbarFuFuPR(gI1,gI2));
      }
      tmp_221 += tmp_222;
   }
   tmp_220 += tmp_221;
   result += (0.5) * tmp_220;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VP(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVPbargWpCgWpC())*B00(p,MVWp,MVWp);
   result += AbsSqr(CpVPbargWpgWp())*B00(p,MVWp,MVWp);
   result += -4*AbsSqr(CpVPconjHpHp())*B00(p,MHp,MHp);
   result += 2*AbsSqr(CpVPconjVWpHp())*B0(p,MVWp,MHp);
   result += A0(MHp)*CpVPVPconjHpHp();
   result += -(A0(MVWp)*(4*CpVPVPconjVWpVWp1() + CpVPVPconjVWpVWp2() +
      CpVPVPconjVWpVWp3()));
   std::complex<double> tmp_223;
   std::complex<double> tmp_224;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_225;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_225 += (AbsSqr(CpVPbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVPbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_225 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVPbarFdFdPL(gI1,gI2))*CpVPbarFdFdPR(gI1,gI2));
      }
      tmp_224 += tmp_225;
   }
   tmp_223 += tmp_224;
   result += (3) * tmp_223;
   std::complex<double> tmp_226;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_227;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_227 += (AbsSqr(CpVPbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVPbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_227 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVPbarFeFePL(gI1,gI2))*CpVPbarFeFePR(gI1,gI2));
      }
      tmp_226 += tmp_227;
   }
   result += tmp_226;
   std::complex<double> tmp_228;
   std::complex<double> tmp_229;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_230;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_230 += (AbsSqr(CpVPbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVPbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_230 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVPbarFuFuPL(gI1,gI2))*CpVPbarFuFuPR(gI1,gI2));
      }
      tmp_229 += tmp_230;
   }
   tmp_228 += tmp_229;
   result += (3) * tmp_228;
   result += (AbsSqr(CpVPbarFcFcPL()) + AbsSqr(CpVPbarFcFcPR()))*H0(p,MFc,MFc);
   result += (AbsSqr(CpVPbarFgFgPL()) + AbsSqr(CpVPbarFgFgPR()))*H0(p,MFg,MFg);
   result += 2*CpVPVPconjVWpVWp1()*Sqr(MVWp);
   result += -(AbsSqr(CpVPconjVWpVWp())*(2*A0(MVWp) + 10*B00(p,MVWp,MVWp) - 2*(
      2*Sqr(MVWp) - 0.3333333333333333*Sqr(p)) + B0(p,MVWp,MVWp)*(2*Sqr(MVWp) + 4*
      Sqr(p))));
   result += 4*B0(p,MFc,MFc)*Re(Conj(CpVPbarFcFcPL())*CpVPbarFcFcPR())*Sqr(MFc)
      ;
   result += 4*B0(p,MFg,MFg)*Re(Conj(CpVPbarFgFgPL())*CpVPbarFgFgPR())*Sqr(MFg)
      ;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVZbargWpCgWpC())*B00(p,MVWp,MVWp);
   result += AbsSqr(CpVZbargWpgWp())*B00(p,MVWp,MVWp);
   result += -4*AbsSqr(CpVZconjHpHp())*B00(p,MHp,MHp);
   result += 2*AbsSqr(CpVZconjVWpHp())*B0(p,MVWp,MHp);
   result += -4*AbsSqr(CpVZhhAh())*B00(p,MAh,Mhh);
   result += 0.5*A0(MAh)*CpVZVZAhAh();
   result += A0(MHp)*CpVZVZconjHpHp();
   result += -(A0(MVWp)*(4*CpVZVZconjVWpVWp1() + CpVZVZconjVWpVWp2() +
      CpVZVZconjVWpVWp3()));
   result += AbsSqr(CpVZVZhh())*B0(p,MVZ,Mhh);
   result += 0.5*A0(Mhh)*CpVZVZhhhh();
   std::complex<double> tmp_231;
   std::complex<double> tmp_232;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_233;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_233 += (AbsSqr(CpVZbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVZbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_233 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVZbarFdFdPL(gI1,gI2))*CpVZbarFdFdPR(gI1,gI2));
      }
      tmp_232 += tmp_233;
   }
   tmp_231 += tmp_232;
   result += (3) * tmp_231;
   std::complex<double> tmp_234;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_235;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_235 += (AbsSqr(CpVZbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVZbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_235 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVZbarFeFePL(gI1,gI2))*CpVZbarFeFePR(gI1,gI2));
      }
      tmp_234 += tmp_235;
   }
   result += tmp_234;
   std::complex<double> tmp_236;
   std::complex<double> tmp_237;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_238;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_238 += (AbsSqr(CpVZbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVZbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_238 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVZbarFuFuPL(gI1,gI2))*CpVZbarFuFuPR(gI1,gI2));
      }
      tmp_237 += tmp_238;
   }
   tmp_236 += tmp_237;
   result += (3) * tmp_236;
   std::complex<double> tmp_239;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_240;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_240 += (AbsSqr(CpVZbarFvFvPL(gI1,gI2)) + AbsSqr(
            CpVZbarFvFvPR(gI1,gI2)))*H0(p,MFv(gI1),MFv(gI2));
         tmp_240 += 4*B0(p,MFv(gI1),MFv(gI2))*MFv(gI1)*MFv(gI2)*Re(Conj(
            CpVZbarFvFvPL(gI1,gI2))*CpVZbarFvFvPR(gI1,gI2));
      }
      tmp_239 += tmp_240;
   }
   result += tmp_239;
   result += (AbsSqr(CpVZbarFcFcPL()) + AbsSqr(CpVZbarFcFcPR()))*H0(p,MFc,MFc);
   result += (AbsSqr(CpVZbarFgFgPL()) + AbsSqr(CpVZbarFgFgPR()))*H0(p,MFg,MFg);
   result += 2*CpVZVZconjVWpVWp1()*Sqr(MVWp);
   result += -(AbsSqr(CpVZconjVWpVWp())*(2*A0(MVWp) + 10*B00(p,MVWp,MVWp) - 2*(
      2*Sqr(MVWp) - 0.3333333333333333*Sqr(p)) + B0(p,MVWp,MVWp)*(2*Sqr(MVWp) + 4*
      Sqr(p))));
   result += 4*B0(p,MFc,MFc)*Re(Conj(CpVZbarFcFcPL())*CpVZbarFcFcPR())*Sqr(MFc)
      ;
   result += 4*B0(p,MFg,MFg)*Re(Conj(CpVZbarFgFgPL())*CpVZbarFgFgPR())*Sqr(MFg)
      ;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWp(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpconjVWpbargPgWp())*B00(p,MVWp,MVP);
   result += AbsSqr(CpconjVWpbargWpCgP())*B00(p,MVP,MVWp);
   result += AbsSqr(CpconjVWpbargWpCgZ())*B00(p,MVZ,MVWp);
   result += AbsSqr(CpconjVWpbargZgWp())*B00(p,MVWp,MVZ);
   result += -4*AbsSqr(CpconjVWpHpAh())*B00(p,MAh,MHp);
   result += -4*AbsSqr(CpconjVWpHphh())*B00(p,Mhh,MHp);
   result += AbsSqr(CpconjVWpVPHp())*B0(p,0,MHp);
   result += AbsSqr(CpconjVWpVWphh())*B0(p,MVWp,Mhh);
   result += AbsSqr(CpconjVWpVZHp())*B0(p,MVZ,MHp);
   result += 0.5*A0(MAh)*CpVWpconjVWpAhAh();
   result += A0(MHp)*CpVWpconjVWpconjHpHp();
   result += -(A0(MVWp)*(4*CpVWpconjVWpconjVWpVWp1() + CpVWpconjVWpconjVWpVWp2(
      ) + CpVWpconjVWpconjVWpVWp3()));
   result += 0.5*A0(Mhh)*CpVWpconjVWphhhh();
   result += 0;
   std::complex<double> tmp_241;
   std::complex<double> tmp_242;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_243;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_243 += (AbsSqr(CpconjVWpbarFdFuPL(gI1,gI2)) + AbsSqr(
            CpconjVWpbarFdFuPR(gI1,gI2)))*H0(p,MFd(gI1),MFu(gI2));
         tmp_243 += 4*B0(p,MFd(gI1),MFu(gI2))*MFd(gI1)*MFu(gI2)*Re(Conj(
            CpconjVWpbarFdFuPL(gI1,gI2))*CpconjVWpbarFdFuPR(gI1,gI2));
      }
      tmp_242 += tmp_243;
   }
   tmp_241 += tmp_242;
   result += (3) * tmp_241;
   std::complex<double> tmp_244;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_245;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_245 += (AbsSqr(CpconjVWpbarFeFvPL(gI1,gI2)) + AbsSqr(
            CpconjVWpbarFeFvPR(gI1,gI2)))*H0(p,MFe(gI1),MFv(gI2));
         tmp_245 += 4*B0(p,MFe(gI1),MFv(gI2))*MFe(gI1)*MFv(gI2)*Re(Conj(
            CpconjVWpbarFeFvPL(gI1,gI2))*CpconjVWpbarFeFvPR(gI1,gI2));
      }
      tmp_244 += tmp_245;
   }
   result += tmp_244;
   result += (AbsSqr(CpconjVWpbarFcFgPL()) + AbsSqr(CpconjVWpbarFcFgPR()))*H0(p
      ,MFc,MFg);
   result += (AbsSqr(CpconjVWpFnFcPL()) + AbsSqr(CpconjVWpFnFcPR()))*H0(p,MFn,
      MFc);
   result += 2*CpVWpconjVWpconjVWpVWp1()*Sqr(MVWp);
   result += -(AbsSqr(CpconjVWpVWpVP())*(A0(MVWp) + 10*B00(p,MVWp,0) - 2*(Sqr(
      MVWp) - 0.3333333333333333*Sqr(p)) + B0(p,MVWp,0)*(Sqr(MVWp) + 4*Sqr(p))));
   result += 0.5*(-(A0(MVZ)*(4*CpVWpconjVWpVZVZ1() + CpVWpconjVWpVZVZ2() +
      CpVWpconjVWpVZVZ3())) + 2*CpVWpconjVWpVZVZ1()*Sqr(MVZ));
   result += -(AbsSqr(CpconjVWpVZVWp())*(A0(MVWp) + A0(MVZ) + 10*B00(p,MVZ,MVWp
      ) - 2*(Sqr(MVWp) + Sqr(MVZ) - 0.3333333333333333*Sqr(p)) + B0(p,MVZ,MVWp)*(
      Sqr(MVWp) + Sqr(MVZ) + 4*Sqr(p))));
   result += 4*MFc*MFg*B0(p,MFc,MFg)*Re(Conj(CpconjVWpbarFcFgPL())*
      CpconjVWpbarFcFgPR());
   result += 4*MFc*MFn*B0(p,MFn,MFc)*Re(Conj(CpconjVWpFnFcPL())*CpconjVWpFnFcPR
      ());

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_246;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_246 += B0(p,MFd(gI1),MAh)*Conj(CpbarUFdFdAhPL(gO2,gI1))*
         CpbarUFdFdAhPR(gO1,gI1)*MFd(gI1);
   }
   result += tmp_246;
   std::complex<double> tmp_247;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_247 += B0(p,MFd(gI2),Mhh)*Conj(CpbarUFdhhFdPL(gO2,gI2))*
         CpbarUFdhhFdPR(gO1,gI2)*MFd(gI2);
   }
   result += tmp_247;
   std::complex<double> tmp_248;
   std::complex<double> tmp_249;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_249 += (-0.5 + B0(p,MFd(gI2),0))*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_248 += tmp_249;
   result += (-5.333333333333333) * tmp_248;
   std::complex<double> tmp_250;
   std::complex<double> tmp_251;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_251 += (-0.5 + B0(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_250 += tmp_251;
   result += (-4) * tmp_250;
   std::complex<double> tmp_252;
   std::complex<double> tmp_253;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_253 += (-0.5 + B0(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_252 += tmp_253;
   result += (-4) * tmp_252;
   std::complex<double> tmp_254;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_254 += B0(p,MFu(gI2),MHp)*Conj(CpbarUFdconjHpFuPL(gO2,gI2))*
         CpbarUFdconjHpFuPR(gO1,gI2)*MFu(gI2);
   }
   result += tmp_254;
   std::complex<double> tmp_255;
   std::complex<double> tmp_256;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_256 += (-0.5 + B0(p,MFu(gI2),MVWp))*Conj(CpbarUFdconjVWpFuPR(gO2,
         gI2))*CpbarUFdconjVWpFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_255 += tmp_256;
   result += (-4) * tmp_255;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_257;
   std::complex<double> tmp_258;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_258 += B1(p,MFd(gI1),MAh)*Conj(CpbarUFdFdAhPR(gO2,gI1))*
         CpbarUFdFdAhPR(gO1,gI1);
   }
   tmp_257 += tmp_258;
   result += (-0.5) * tmp_257;
   std::complex<double> tmp_259;
   std::complex<double> tmp_260;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_260 += B1(p,MFu(gI2),MHp)*Conj(CpbarUFdconjHpFuPR(gO2,gI2))*
         CpbarUFdconjHpFuPR(gO1,gI2);
   }
   tmp_259 += tmp_260;
   result += (-0.5) * tmp_259;
   std::complex<double> tmp_261;
   std::complex<double> tmp_262;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_262 += (0.5 + B1(p,MFu(gI2),MVWp))*Conj(CpbarUFdconjVWpFuPL(gO2,
         gI2))*CpbarUFdconjVWpFuPL(gO1,gI2);
   }
   tmp_261 += tmp_262;
   result += (-1) * tmp_261;
   std::complex<double> tmp_263;
   std::complex<double> tmp_264;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_264 += B1(p,MFd(gI2),Mhh)*Conj(CpbarUFdhhFdPR(gO2,gI2))*
         CpbarUFdhhFdPR(gO1,gI2);
   }
   tmp_263 += tmp_264;
   result += (-0.5) * tmp_263;
   std::complex<double> tmp_265;
   std::complex<double> tmp_266;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_266 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVGFdPL(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2);
   }
   tmp_265 += tmp_266;
   result += (-1.3333333333333333) * tmp_265;
   std::complex<double> tmp_267;
   std::complex<double> tmp_268;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_268 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPL(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2);
   }
   tmp_267 += tmp_268;
   result += (-1) * tmp_267;
   std::complex<double> tmp_269;
   std::complex<double> tmp_270;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_270 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPL(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2);
   }
   tmp_269 += tmp_270;
   result += (-1) * tmp_269;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_271;
   std::complex<double> tmp_272;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_272 += B1(p,MFd(gI1),MAh)*Conj(CpbarUFdFdAhPL(gO2,gI1))*
         CpbarUFdFdAhPL(gO1,gI1);
   }
   tmp_271 += tmp_272;
   result += (-0.5) * tmp_271;
   std::complex<double> tmp_273;
   std::complex<double> tmp_274;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_274 += B1(p,MFu(gI2),MHp)*Conj(CpbarUFdconjHpFuPL(gO2,gI2))*
         CpbarUFdconjHpFuPL(gO1,gI2);
   }
   tmp_273 += tmp_274;
   result += (-0.5) * tmp_273;
   std::complex<double> tmp_275;
   std::complex<double> tmp_276;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_276 += (0.5 + B1(p,MFu(gI2),MVWp))*Conj(CpbarUFdconjVWpFuPR(gO2,
         gI2))*CpbarUFdconjVWpFuPR(gO1,gI2);
   }
   tmp_275 += tmp_276;
   result += (-1) * tmp_275;
   std::complex<double> tmp_277;
   std::complex<double> tmp_278;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_278 += B1(p,MFd(gI2),Mhh)*Conj(CpbarUFdhhFdPL(gO2,gI2))*
         CpbarUFdhhFdPL(gO1,gI2);
   }
   tmp_277 += tmp_278;
   result += (-0.5) * tmp_277;
   std::complex<double> tmp_279;
   std::complex<double> tmp_280;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_280 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPR(gO1,gI2);
   }
   tmp_279 += tmp_280;
   result += (-1.3333333333333333) * tmp_279;
   std::complex<double> tmp_281;
   std::complex<double> tmp_282;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_282 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPR(gO1,gI2);
   }
   tmp_281 += tmp_282;
   result += (-1) * tmp_281;
   std::complex<double> tmp_283;
   std::complex<double> tmp_284;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_284 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPR(gO1,gI2);
   }
   tmp_283 += tmp_284;
   result += (-1) * tmp_283;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_285;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_285 += B0(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPL(gO2,gI1))*
         CpbarUFuFuAhPR(gO1,gI1)*MFu(gI1);
   }
   result += tmp_285;
   std::complex<double> tmp_286;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_286 += B0(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPL(gO2,gI2))*
         CpbarUFuHpFdPR(gO1,gI2)*MFd(gI2);
   }
   result += tmp_286;
   std::complex<double> tmp_287;
   std::complex<double> tmp_288;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_288 += (-0.5 + B0(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPR(gO2,gI2))
         *CpbarUFuVWpFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_287 += tmp_288;
   result += (-4) * tmp_287;
   std::complex<double> tmp_289;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_289 += B0(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPL(gO2,gI2))*
         CpbarUFuhhFuPR(gO1,gI2)*MFu(gI2);
   }
   result += tmp_289;
   std::complex<double> tmp_290;
   std::complex<double> tmp_291;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_291 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_290 += tmp_291;
   result += (-5.333333333333333) * tmp_290;
   std::complex<double> tmp_292;
   std::complex<double> tmp_293;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_293 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_292 += tmp_293;
   result += (-4) * tmp_292;
   std::complex<double> tmp_294;
   std::complex<double> tmp_295;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_295 += (-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_294 += tmp_295;
   result += (-4) * tmp_294;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_296;
   std::complex<double> tmp_297;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_297 += B1(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPR(gO2,gI1))*
         CpbarUFuFuAhPR(gO1,gI1);
   }
   tmp_296 += tmp_297;
   result += (-0.5) * tmp_296;
   std::complex<double> tmp_298;
   std::complex<double> tmp_299;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_299 += B1(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPR(gO2,gI2))*
         CpbarUFuhhFuPR(gO1,gI2);
   }
   tmp_298 += tmp_299;
   result += (-0.5) * tmp_298;
   std::complex<double> tmp_300;
   std::complex<double> tmp_301;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_301 += B1(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPR(gO2,gI2))*
         CpbarUFuHpFdPR(gO1,gI2);
   }
   tmp_300 += tmp_301;
   result += (-0.5) * tmp_300;
   std::complex<double> tmp_302;
   std::complex<double> tmp_303;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_303 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVGFuPL(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2);
   }
   tmp_302 += tmp_303;
   result += (-1.3333333333333333) * tmp_302;
   std::complex<double> tmp_304;
   std::complex<double> tmp_305;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_305 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_304 += tmp_305;
   result += (-1) * tmp_304;
   std::complex<double> tmp_306;
   std::complex<double> tmp_307;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_307 += (0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPL(gO2,gI2))*
         CpbarUFuVWpFdPL(gO1,gI2);
   }
   tmp_306 += tmp_307;
   result += (-1) * tmp_306;
   std::complex<double> tmp_308;
   std::complex<double> tmp_309;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_309 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_308 += tmp_309;
   result += (-1) * tmp_308;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_310;
   std::complex<double> tmp_311;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_311 += B1(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPL(gO2,gI1))*
         CpbarUFuFuAhPL(gO1,gI1);
   }
   tmp_310 += tmp_311;
   result += (-0.5) * tmp_310;
   std::complex<double> tmp_312;
   std::complex<double> tmp_313;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_313 += B1(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPL(gO2,gI2))*
         CpbarUFuhhFuPL(gO1,gI2);
   }
   tmp_312 += tmp_313;
   result += (-0.5) * tmp_312;
   std::complex<double> tmp_314;
   std::complex<double> tmp_315;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_315 += B1(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPL(gO2,gI2))*
         CpbarUFuHpFdPL(gO1,gI2);
   }
   tmp_314 += tmp_315;
   result += (-0.5) * tmp_314;
   std::complex<double> tmp_316;
   std::complex<double> tmp_317;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_317 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPR(gO1,gI2);
   }
   tmp_316 += tmp_317;
   result += (-1.3333333333333333) * tmp_316;
   std::complex<double> tmp_318;
   std::complex<double> tmp_319;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_319 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_318 += tmp_319;
   result += (-1) * tmp_318;
   std::complex<double> tmp_320;
   std::complex<double> tmp_321;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_321 += (0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPR(gO2,gI2))*
         CpbarUFuVWpFdPR(gO1,gI2);
   }
   tmp_320 += tmp_321;
   result += (-1) * tmp_320;
   std::complex<double> tmp_322;
   std::complex<double> tmp_323;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_323 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_322 += tmp_323;
   result += (-1) * tmp_322;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_324;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_324 += B0(p,MFe(gI1),MAh)*Conj(CpbarUFeFeAhPL(gO2,gI1))*
         CpbarUFeFeAhPR(gO1,gI1)*MFe(gI1);
   }
   result += tmp_324;
   std::complex<double> tmp_325;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_325 += B0(p,MFe(gI2),Mhh)*Conj(CpbarUFehhFePL(gO2,gI2))*
         CpbarUFehhFePR(gO1,gI2)*MFe(gI2);
   }
   result += tmp_325;
   std::complex<double> tmp_326;
   std::complex<double> tmp_327;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_327 += (-0.5 + B0(p,MFe(gI2),0))*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_326 += tmp_327;
   result += (-4) * tmp_326;
   std::complex<double> tmp_328;
   std::complex<double> tmp_329;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_329 += (-0.5 + B0(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_328 += tmp_329;
   result += (-4) * tmp_328;
   std::complex<double> tmp_330;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_330 += B0(p,MFv(gI2),MHp)*Conj(CpbarUFeconjHpFvPL(gO2,gI2))*
         CpbarUFeconjHpFvPR(gO1,gI2)*MFv(gI2);
   }
   result += tmp_330;
   std::complex<double> tmp_331;
   std::complex<double> tmp_332;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_332 += (-0.5 + B0(p,MFv(gI2),MVWp))*Conj(CpbarUFeconjVWpFvPR(gO2,
         gI2))*CpbarUFeconjVWpFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_331 += tmp_332;
   result += (-4) * tmp_331;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_333;
   std::complex<double> tmp_334;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_334 += B1(p,MFe(gI1),MAh)*Conj(CpbarUFeFeAhPR(gO2,gI1))*
         CpbarUFeFeAhPR(gO1,gI1);
   }
   tmp_333 += tmp_334;
   result += (-0.5) * tmp_333;
   std::complex<double> tmp_335;
   std::complex<double> tmp_336;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_336 += B1(p,MFv(gI2),MHp)*Conj(CpbarUFeconjHpFvPR(gO2,gI2))*
         CpbarUFeconjHpFvPR(gO1,gI2);
   }
   tmp_335 += tmp_336;
   result += (-0.5) * tmp_335;
   std::complex<double> tmp_337;
   std::complex<double> tmp_338;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_338 += (0.5 + B1(p,MFv(gI2),MVWp))*Conj(CpbarUFeconjVWpFvPL(gO2,
         gI2))*CpbarUFeconjVWpFvPL(gO1,gI2);
   }
   tmp_337 += tmp_338;
   result += (-1) * tmp_337;
   std::complex<double> tmp_339;
   std::complex<double> tmp_340;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_340 += B1(p,MFe(gI2),Mhh)*Conj(CpbarUFehhFePR(gO2,gI2))*
         CpbarUFehhFePR(gO1,gI2);
   }
   tmp_339 += tmp_340;
   result += (-0.5) * tmp_339;
   std::complex<double> tmp_341;
   std::complex<double> tmp_342;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_342 += (0.5 + B1(p,MFe(gI2),0))*Conj(CpbarUFeVPFePL(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2);
   }
   tmp_341 += tmp_342;
   result += (-1) * tmp_341;
   std::complex<double> tmp_343;
   std::complex<double> tmp_344;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_344 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePL(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2);
   }
   tmp_343 += tmp_344;
   result += (-1) * tmp_343;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_345;
   std::complex<double> tmp_346;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_346 += B1(p,MFe(gI1),MAh)*Conj(CpbarUFeFeAhPL(gO2,gI1))*
         CpbarUFeFeAhPL(gO1,gI1);
   }
   tmp_345 += tmp_346;
   result += (-0.5) * tmp_345;
   std::complex<double> tmp_347;
   std::complex<double> tmp_348;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_348 += B1(p,MFv(gI2),MHp)*Conj(CpbarUFeconjHpFvPL(gO2,gI2))*
         CpbarUFeconjHpFvPL(gO1,gI2);
   }
   tmp_347 += tmp_348;
   result += (-0.5) * tmp_347;
   std::complex<double> tmp_349;
   std::complex<double> tmp_350;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_350 += (0.5 + B1(p,MFv(gI2),MVWp))*Conj(CpbarUFeconjVWpFvPR(gO2,
         gI2))*CpbarUFeconjVWpFvPR(gO1,gI2);
   }
   tmp_349 += tmp_350;
   result += (-1) * tmp_349;
   std::complex<double> tmp_351;
   std::complex<double> tmp_352;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_352 += B1(p,MFe(gI2),Mhh)*Conj(CpbarUFehhFePL(gO2,gI2))*
         CpbarUFehhFePL(gO1,gI2);
   }
   tmp_351 += tmp_352;
   result += (-0.5) * tmp_351;
   std::complex<double> tmp_353;
   std::complex<double> tmp_354;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_354 += (0.5 + B1(p,MFe(gI2),0))*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePR(gO1,gI2);
   }
   tmp_353 += tmp_354;
   result += (-1) * tmp_353;
   std::complex<double> tmp_355;
   std::complex<double> tmp_356;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_356 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePR(gO1,gI2);
   }
   tmp_355 += tmp_356;
   result += (-1) * tmp_355;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_357;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_357 += B0(p,MFe(gI2),MHp)*Conj(CpbarFvHpFePL(gO2,gI2))*
         CpbarFvHpFePR(gO1,gI2)*MFe(gI2);
   }
   result += tmp_357;
   std::complex<double> tmp_358;
   std::complex<double> tmp_359;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_359 += (-0.5 + B0(p,MFe(gI2),MVWp))*Conj(CpbarFvVWpFePR(gO2,gI2))*
         CpbarFvVWpFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_358 += tmp_359;
   result += (-4) * tmp_358;
   std::complex<double> tmp_360;
   std::complex<double> tmp_361;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_361 += (-0.5 + B0(p,MFv(gI2),MVZ))*Conj(CpbarFvVZFvPR(gO2,gI2))*
         CpbarFvVZFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_360 += tmp_361;
   result += (-4) * tmp_360;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_362;
   std::complex<double> tmp_363;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_363 += B1(p,MFe(gI2),MHp)*Conj(CpbarFvHpFePR(gO2,gI2))*
         CpbarFvHpFePR(gO1,gI2);
   }
   tmp_362 += tmp_363;
   result += (-0.5) * tmp_362;
   std::complex<double> tmp_364;
   std::complex<double> tmp_365;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_365 += (0.5 + B1(p,MFe(gI2),MVWp))*Conj(CpbarFvVWpFePL(gO2,gI2))*
         CpbarFvVWpFePL(gO1,gI2);
   }
   tmp_364 += tmp_365;
   result += (-1) * tmp_364;
   std::complex<double> tmp_366;
   std::complex<double> tmp_367;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_367 += (0.5 + B1(p,MFv(gI2),MVZ))*Conj(CpbarFvVZFvPL(gO2,gI2))*
         CpbarFvVZFvPL(gO1,gI2);
   }
   tmp_366 += tmp_367;
   result += (-1) * tmp_366;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_368;
   std::complex<double> tmp_369;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_369 += B1(p,MFe(gI2),MHp)*Conj(CpbarFvHpFePL(gO2,gI2))*
         CpbarFvHpFePL(gO1,gI2);
   }
   tmp_368 += tmp_369;
   result += (-0.5) * tmp_368;
   std::complex<double> tmp_370;
   std::complex<double> tmp_371;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_371 += (0.5 + B1(p,MFe(gI2),MVWp))*Conj(CpbarFvVWpFePR(gO2,gI2))*
         CpbarFvVWpFePR(gO1,gI2);
   }
   tmp_370 += tmp_371;
   result += (-1) * tmp_370;
   std::complex<double> tmp_372;
   std::complex<double> tmp_373;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_373 += (0.5 + B1(p,MFv(gI2),MVZ))*Conj(CpbarFvVZFvPR(gO2,gI2))*
         CpbarFvVZFvPR(gO1,gI2);
   }
   tmp_372 += tmp_373;
   result += (-1) * tmp_372;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fc_1(double p ) const
{
   std::complex<double> result;

   result += -4*MFc*(-0.5 + B0(p,MFc,0))*Conj(CpbarFcVPFcPR())*CpbarFcVPFcPL();
   result += -4*MFc*(-0.5 + B0(p,MFc,MVZ))*Conj(CpbarFcVZFcPR())*CpbarFcVZFcPL(
      );
   result += -4*MFg*(-0.5 + B0(p,MFg,MVWp))*Conj(CpbarFcconjVWpFgPR())*
      CpbarFcconjVWpFgPL();
   result += -4*MFn*(-0.5 + B0(p,MFn,MVWp))*Conj(CpbarFcVWpFnPR())*
      CpbarFcVWpFnPL();

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fc_PR(double p ) const
{
   std::complex<double> result;

   result += -(AbsSqr(CpbarFcconjVWpFgPL())*(0.5 + B1(p,MFg,MVWp)));
   result += -(AbsSqr(CpbarFcVPFcPL())*(0.5 + B1(p,MFc,0)));
   result += -(AbsSqr(CpbarFcVWpFnPL())*(0.5 + B1(p,MFn,MVWp)));
   result += -(AbsSqr(CpbarFcVZFcPL())*(0.5 + B1(p,MFc,MVZ)));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fc_PL(double p ) const
{
   std::complex<double> result;

   result += -(AbsSqr(CpbarFcconjVWpFgPR())*(0.5 + B1(p,MFg,MVWp)));
   result += -(AbsSqr(CpbarFcVPFcPR())*(0.5 + B1(p,MFc,0)));
   result += -(AbsSqr(CpbarFcVWpFnPR())*(0.5 + B1(p,MFn,MVWp)));
   result += -(AbsSqr(CpbarFcVZFcPR())*(0.5 + B1(p,MFc,MVZ)));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fg_1(double p ) const
{
   std::complex<double> result;

   result += -4*MFc*(-0.5 + B0(p,MFc,MVWp))*Conj(CpbarFgVWpFcPR())*
      CpbarFgVWpFcPL();
   result += -4*MFg*(-0.5 + B0(p,MFg,0))*Conj(CpbarFgVPFgPR())*CpbarFgVPFgPL();
   result += -4*MFg*(-0.5 + B0(p,MFg,MVZ))*Conj(CpbarFgVZFgPR())*CpbarFgVZFgPL(
      );

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fg_PR(double p ) const
{
   std::complex<double> result;

   result += -(AbsSqr(CpbarFgVPFgPL())*(0.5 + B1(p,MFg,0)));
   result += -(AbsSqr(CpbarFgVWpFcPL())*(0.5 + B1(p,MFc,MVWp)));
   result += -(AbsSqr(CpbarFgVZFgPL())*(0.5 + B1(p,MFg,MVZ)));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fg_PL(double p ) const
{
   std::complex<double> result;

   result += -(AbsSqr(CpbarFgVPFgPR())*(0.5 + B1(p,MFg,0)));
   result += -(AbsSqr(CpbarFgVWpFcPR())*(0.5 + B1(p,MFc,MVWp)));
   result += -(AbsSqr(CpbarFgVZFgPR())*(0.5 + B1(p,MFg,MVZ)));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fn_1(double p ) const
{
   std::complex<double> result;

   result += -4*MFc*(-0.5 + B0(p,MFc,MVWp))*Conj(CpFnbarFcVWpPR())*
      CpFnbarFcVWpPL();
   result += -4*MFc*(-0.5 + B0(p,MFc,MVWp))*Conj(CpFnconjVWpFcPR())*
      CpFnconjVWpFcPL();

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fn_PR(double p ) const
{
   std::complex<double> result;

   result += -(AbsSqr(CpFnbarFcVWpPL())*(0.5 + B1(p,MFc,MVWp)));
   result += -(AbsSqr(CpFnconjVWpFcPL())*(0.5 + B1(p,MFc,MVWp)));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fn_PL(double p ) const
{
   std::complex<double> result;

   result += -(AbsSqr(CpFnbarFcVWpPR())*(0.5 + B1(p,MFc,MVWp)));
   result += -(AbsSqr(CpFnconjVWpFcPR())*(0.5 + B1(p,MFc,MVWp)));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ_heavy(double p ) const
{
   std::complex<double> result;

   result += (AbsSqr(CpVZbarFcFcPL()) + AbsSqr(CpVZbarFcFcPR()))*H0(p,MFc,MFc);
   result += (AbsSqr(CpVZbarFgFgPL()) + AbsSqr(CpVZbarFgFgPR()))*H0(p,MFg,MFg);
   result += 4*B0(p,MFc,MFc)*Re(Conj(CpVZbarFcFcPL())*CpVZbarFcFcPR())*Sqr(MFc)
      ;
   result += 4*B0(p,MFg,MFg)*Re(Conj(CpVZbarFgFgPL())*CpVZbarFgFgPR())*Sqr(MFg)
      ;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWp_heavy(double p ) const
{
   std::complex<double> result;

   result += (AbsSqr(CpconjVWpbarFcFgPL()) + AbsSqr(CpconjVWpbarFcFgPR()))*H0(p
      ,MFc,MFg);
   result += (AbsSqr(CpconjVWpFnFcPL()) + AbsSqr(CpconjVWpFnFcPR()))*H0(p,MFn,
      MFc);
   result += 4*MFc*MFg*B0(p,MFc,MFg)*Re(Conj(CpconjVWpbarFcFgPL())*
      CpconjVWpbarFcFgPR());
   result += 4*MFc*MFn*B0(p,MFn,MFc)*Re(Conj(CpconjVWpFnFcPL())*CpconjVWpFnFcPR
      ());

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_374;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_374 += B0(p,MFd(gI1),MAh)*Conj(CpbarFdFdAhPL(gO2,gI1))*
         CpbarFdFdAhPR(gO1,gI1)*MFd(gI1);
   }
   result += tmp_374;
   std::complex<double> tmp_375;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_375 += B0(p,MFd(gI2),Mhh)*Conj(CpbarFdhhFdPL(gO2,gI2))*
         CpbarFdhhFdPR(gO1,gI2)*MFd(gI2);
   }
   result += tmp_375;
   std::complex<double> tmp_376;
   std::complex<double> tmp_377;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_377 += (-0.5 + B0(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_376 += tmp_377;
   result += (-4) * tmp_376;
   std::complex<double> tmp_378;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_378 += B0(p,MFu(gI2),MHp)*Conj(CpbarFdconjHpFuPL(gO2,gI2))*
         CpbarFdconjHpFuPR(gO1,gI2)*MFu(gI2);
   }
   result += tmp_378;
   std::complex<double> tmp_379;
   std::complex<double> tmp_380;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_380 += (-0.5 + B0(p,MFu(gI2),MVWp))*Conj(CpbarFdconjVWpFuPR(gO2,
         gI2))*CpbarFdconjVWpFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_379 += tmp_380;
   result += (-4) * tmp_379;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_381;
   std::complex<double> tmp_382;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_382 += B1(p,MFd(gI1),MAh)*Conj(CpbarFdFdAhPR(gO2,gI1))*
         CpbarFdFdAhPR(gO1,gI1);
   }
   tmp_381 += tmp_382;
   result += (-0.5) * tmp_381;
   std::complex<double> tmp_383;
   std::complex<double> tmp_384;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_384 += B1(p,MFu(gI2),MHp)*Conj(CpbarFdconjHpFuPR(gO2,gI2))*
         CpbarFdconjHpFuPR(gO1,gI2);
   }
   tmp_383 += tmp_384;
   result += (-0.5) * tmp_383;
   std::complex<double> tmp_385;
   std::complex<double> tmp_386;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_386 += (0.5 + B1(p,MFu(gI2),MVWp))*Conj(CpbarFdconjVWpFuPL(gO2,gI2
         ))*CpbarFdconjVWpFuPL(gO1,gI2);
   }
   tmp_385 += tmp_386;
   result += (-1) * tmp_385;
   std::complex<double> tmp_387;
   std::complex<double> tmp_388;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_388 += B1(p,MFd(gI2),Mhh)*Conj(CpbarFdhhFdPR(gO2,gI2))*
         CpbarFdhhFdPR(gO1,gI2);
   }
   tmp_387 += tmp_388;
   result += (-0.5) * tmp_387;
   std::complex<double> tmp_389;
   std::complex<double> tmp_390;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_390 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPL(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2);
   }
   tmp_389 += tmp_390;
   result += (-1) * tmp_389;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_391;
   std::complex<double> tmp_392;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_392 += B1(p,MFd(gI1),MAh)*Conj(CpbarFdFdAhPL(gO2,gI1))*
         CpbarFdFdAhPL(gO1,gI1);
   }
   tmp_391 += tmp_392;
   result += (-0.5) * tmp_391;
   std::complex<double> tmp_393;
   std::complex<double> tmp_394;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_394 += B1(p,MFu(gI2),MHp)*Conj(CpbarFdconjHpFuPL(gO2,gI2))*
         CpbarFdconjHpFuPL(gO1,gI2);
   }
   tmp_393 += tmp_394;
   result += (-0.5) * tmp_393;
   std::complex<double> tmp_395;
   std::complex<double> tmp_396;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_396 += (0.5 + B1(p,MFu(gI2),MVWp))*Conj(CpbarFdconjVWpFuPR(gO2,gI2
         ))*CpbarFdconjVWpFuPR(gO1,gI2);
   }
   tmp_395 += tmp_396;
   result += (-1) * tmp_395;
   std::complex<double> tmp_397;
   std::complex<double> tmp_398;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_398 += B1(p,MFd(gI2),Mhh)*Conj(CpbarFdhhFdPL(gO2,gI2))*
         CpbarFdhhFdPL(gO1,gI2);
   }
   tmp_397 += tmp_398;
   result += (-0.5) * tmp_397;
   std::complex<double> tmp_399;
   std::complex<double> tmp_400;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_400 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPR(gO1,gI2);
   }
   tmp_399 += tmp_400;
   result += (-1) * tmp_399;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_401;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_401 += B0(p,MFe(gI1),MAh)*Conj(CpbarFeFeAhPL(gO2,gI1))*
         CpbarFeFeAhPR(gO1,gI1)*MFe(gI1);
   }
   result += tmp_401;
   std::complex<double> tmp_402;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_402 += B0(p,MFe(gI2),Mhh)*Conj(CpbarFehhFePL(gO2,gI2))*
         CpbarFehhFePR(gO1,gI2)*MFe(gI2);
   }
   result += tmp_402;
   std::complex<double> tmp_403;
   std::complex<double> tmp_404;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_404 += (-0.5 + B0(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_403 += tmp_404;
   result += (-4) * tmp_403;
   std::complex<double> tmp_405;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_405 += B0(p,MFv(gI2),MHp)*Conj(CpbarFeconjHpFvPL(gO2,gI2))*
         CpbarFeconjHpFvPR(gO1,gI2)*MFv(gI2);
   }
   result += tmp_405;
   std::complex<double> tmp_406;
   std::complex<double> tmp_407;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_407 += (-0.5 + B0(p,MFv(gI2),MVWp))*Conj(CpbarFeconjVWpFvPR(gO2,
         gI2))*CpbarFeconjVWpFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_406 += tmp_407;
   result += (-4) * tmp_406;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_408;
   std::complex<double> tmp_409;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_409 += B1(p,MFe(gI1),MAh)*Conj(CpbarFeFeAhPR(gO2,gI1))*
         CpbarFeFeAhPR(gO1,gI1);
   }
   tmp_408 += tmp_409;
   result += (-0.5) * tmp_408;
   std::complex<double> tmp_410;
   std::complex<double> tmp_411;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_411 += B1(p,MFv(gI2),MHp)*Conj(CpbarFeconjHpFvPR(gO2,gI2))*
         CpbarFeconjHpFvPR(gO1,gI2);
   }
   tmp_410 += tmp_411;
   result += (-0.5) * tmp_410;
   std::complex<double> tmp_412;
   std::complex<double> tmp_413;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_413 += (0.5 + B1(p,MFv(gI2),MVWp))*Conj(CpbarFeconjVWpFvPL(gO2,gI2
         ))*CpbarFeconjVWpFvPL(gO1,gI2);
   }
   tmp_412 += tmp_413;
   result += (-1) * tmp_412;
   std::complex<double> tmp_414;
   std::complex<double> tmp_415;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_415 += B1(p,MFe(gI2),Mhh)*Conj(CpbarFehhFePR(gO2,gI2))*
         CpbarFehhFePR(gO1,gI2);
   }
   tmp_414 += tmp_415;
   result += (-0.5) * tmp_414;
   std::complex<double> tmp_416;
   std::complex<double> tmp_417;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_417 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePL(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2);
   }
   tmp_416 += tmp_417;
   result += (-1) * tmp_416;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_418;
   std::complex<double> tmp_419;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_419 += B1(p,MFe(gI1),MAh)*Conj(CpbarFeFeAhPL(gO2,gI1))*
         CpbarFeFeAhPL(gO1,gI1);
   }
   tmp_418 += tmp_419;
   result += (-0.5) * tmp_418;
   std::complex<double> tmp_420;
   std::complex<double> tmp_421;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_421 += B1(p,MFv(gI2),MHp)*Conj(CpbarFeconjHpFvPL(gO2,gI2))*
         CpbarFeconjHpFvPL(gO1,gI2);
   }
   tmp_420 += tmp_421;
   result += (-0.5) * tmp_420;
   std::complex<double> tmp_422;
   std::complex<double> tmp_423;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_423 += (0.5 + B1(p,MFv(gI2),MVWp))*Conj(CpbarFeconjVWpFvPR(gO2,gI2
         ))*CpbarFeconjVWpFvPR(gO1,gI2);
   }
   tmp_422 += tmp_423;
   result += (-1) * tmp_422;
   std::complex<double> tmp_424;
   std::complex<double> tmp_425;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_425 += B1(p,MFe(gI2),Mhh)*Conj(CpbarFehhFePL(gO2,gI2))*
         CpbarFehhFePL(gO1,gI2);
   }
   tmp_424 += tmp_425;
   result += (-0.5) * tmp_424;
   std::complex<double> tmp_426;
   std::complex<double> tmp_427;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_427 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePR(gO1,gI2);
   }
   tmp_426 += tmp_427;
   result += (-1) * tmp_426;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_428;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_428 += B0(p,MFu(gI1),MAh)*Conj(CpbarFuFuAhPL(gO2,gI1))*
         CpbarFuFuAhPR(gO1,gI1)*MFu(gI1);
   }
   result += tmp_428;
   std::complex<double> tmp_429;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_429 += B0(p,MFd(gI2),MHp)*Conj(CpbarFuHpFdPL(gO2,gI2))*
         CpbarFuHpFdPR(gO1,gI2)*MFd(gI2);
   }
   result += tmp_429;
   std::complex<double> tmp_430;
   std::complex<double> tmp_431;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_431 += (-0.5 + B0(p,MFd(gI2),MVWp))*Conj(CpbarFuVWpFdPR(gO2,gI2))*
         CpbarFuVWpFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_430 += tmp_431;
   result += (-4) * tmp_430;
   std::complex<double> tmp_432;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_432 += B0(p,MFu(gI2),Mhh)*Conj(CpbarFuhhFuPL(gO2,gI2))*
         CpbarFuhhFuPR(gO1,gI2)*MFu(gI2);
   }
   result += tmp_432;
   std::complex<double> tmp_433;
   std::complex<double> tmp_434;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_434 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_433 += tmp_434;
   result += (-4) * tmp_433;
   std::complex<double> tmp_435;
   std::complex<double> tmp_436;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_436 += (-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_435 += tmp_436;
   result += (-4) * tmp_435;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_437;
   std::complex<double> tmp_438;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_438 += B1(p,MFu(gI1),MAh)*Conj(CpbarFuFuAhPR(gO2,gI1))*
         CpbarFuFuAhPR(gO1,gI1);
   }
   tmp_437 += tmp_438;
   result += (-0.5) * tmp_437;
   std::complex<double> tmp_439;
   std::complex<double> tmp_440;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_440 += B1(p,MFu(gI2),Mhh)*Conj(CpbarFuhhFuPR(gO2,gI2))*
         CpbarFuhhFuPR(gO1,gI2);
   }
   tmp_439 += tmp_440;
   result += (-0.5) * tmp_439;
   std::complex<double> tmp_441;
   std::complex<double> tmp_442;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_442 += B1(p,MFd(gI2),MHp)*Conj(CpbarFuHpFdPR(gO2,gI2))*
         CpbarFuHpFdPR(gO1,gI2);
   }
   tmp_441 += tmp_442;
   result += (-0.5) * tmp_441;
   std::complex<double> tmp_443;
   std::complex<double> tmp_444;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_444 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarFuVPFuPL(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2);
   }
   tmp_443 += tmp_444;
   result += (-1) * tmp_443;
   std::complex<double> tmp_445;
   std::complex<double> tmp_446;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_446 += (0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarFuVWpFdPL(gO2,gI2))*
         CpbarFuVWpFdPL(gO1,gI2);
   }
   tmp_445 += tmp_446;
   result += (-1) * tmp_445;
   std::complex<double> tmp_447;
   std::complex<double> tmp_448;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_448 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPL(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2);
   }
   tmp_447 += tmp_448;
   result += (-1) * tmp_447;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_449;
   std::complex<double> tmp_450;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_450 += B1(p,MFu(gI1),MAh)*Conj(CpbarFuFuAhPL(gO2,gI1))*
         CpbarFuFuAhPL(gO1,gI1);
   }
   tmp_449 += tmp_450;
   result += (-0.5) * tmp_449;
   std::complex<double> tmp_451;
   std::complex<double> tmp_452;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_452 += B1(p,MFu(gI2),Mhh)*Conj(CpbarFuhhFuPL(gO2,gI2))*
         CpbarFuhhFuPL(gO1,gI2);
   }
   tmp_451 += tmp_452;
   result += (-0.5) * tmp_451;
   std::complex<double> tmp_453;
   std::complex<double> tmp_454;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_454 += B1(p,MFd(gI2),MHp)*Conj(CpbarFuHpFdPL(gO2,gI2))*
         CpbarFuHpFdPL(gO1,gI2);
   }
   tmp_453 += tmp_454;
   result += (-0.5) * tmp_453;
   std::complex<double> tmp_455;
   std::complex<double> tmp_456;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_456 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPR(gO1,gI2);
   }
   tmp_455 += tmp_456;
   result += (-1) * tmp_455;
   std::complex<double> tmp_457;
   std::complex<double> tmp_458;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_458 += (0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarFuVWpFdPR(gO2,gI2))*
         CpbarFuVWpFdPR(gO1,gI2);
   }
   tmp_457 += tmp_458;
   result += (-1) * tmp_457;
   std::complex<double> tmp_459;
   std::complex<double> tmp_460;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_460 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPR(gO1,gI2);
   }
   tmp_459 += tmp_460;
   result += (-1) * tmp_459;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_461;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_461 += B0(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPL(gO2,gI1))*
         CpbarUFuFuAhPR(gO1,gI1)*MFu(gI1);
   }
   result += tmp_461;
   std::complex<double> tmp_462;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_462 += B0(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPL(gO2,gI2))*
         CpbarUFuHpFdPR(gO1,gI2)*MFd(gI2);
   }
   result += tmp_462;
   std::complex<double> tmp_463;
   std::complex<double> tmp_464;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_464 += (-0.5 + B0(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPR(gO2,gI2))
         *CpbarUFuVWpFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_463 += tmp_464;
   result += (-4) * tmp_463;
   std::complex<double> tmp_465;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_465 += B0(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPL(gO2,gI2))*
         CpbarUFuhhFuPR(gO1,gI2)*MFu(gI2);
   }
   result += tmp_465;
   std::complex<double> tmp_466;
   std::complex<double> tmp_467;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_467 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_466 += tmp_467;
   result += (-4) * tmp_466;
   std::complex<double> tmp_468;
   std::complex<double> tmp_469;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_469 += (-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_468 += tmp_469;
   result += (-4) * tmp_468;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_470;
   std::complex<double> tmp_471;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_471 += B1(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPR(gO2,gI1))*
         CpbarUFuFuAhPR(gO1,gI1);
   }
   tmp_470 += tmp_471;
   result += (-0.5) * tmp_470;
   std::complex<double> tmp_472;
   std::complex<double> tmp_473;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_473 += B1(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPR(gO2,gI2))*
         CpbarUFuhhFuPR(gO1,gI2);
   }
   tmp_472 += tmp_473;
   result += (-0.5) * tmp_472;
   std::complex<double> tmp_474;
   std::complex<double> tmp_475;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_475 += B1(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPR(gO2,gI2))*
         CpbarUFuHpFdPR(gO1,gI2);
   }
   tmp_474 += tmp_475;
   result += (-0.5) * tmp_474;
   std::complex<double> tmp_476;
   std::complex<double> tmp_477;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_477 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_476 += tmp_477;
   result += (-1) * tmp_476;
   std::complex<double> tmp_478;
   std::complex<double> tmp_479;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_479 += (0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPL(gO2,gI2))*
         CpbarUFuVWpFdPL(gO1,gI2);
   }
   tmp_478 += tmp_479;
   result += (-1) * tmp_478;
   std::complex<double> tmp_480;
   std::complex<double> tmp_481;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_481 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_480 += tmp_481;
   result += (-1) * tmp_480;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_482;
   std::complex<double> tmp_483;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_483 += B1(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPL(gO2,gI1))*
         CpbarUFuFuAhPL(gO1,gI1);
   }
   tmp_482 += tmp_483;
   result += (-0.5) * tmp_482;
   std::complex<double> tmp_484;
   std::complex<double> tmp_485;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_485 += B1(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPL(gO2,gI2))*
         CpbarUFuhhFuPL(gO1,gI2);
   }
   tmp_484 += tmp_485;
   result += (-0.5) * tmp_484;
   std::complex<double> tmp_486;
   std::complex<double> tmp_487;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_487 += B1(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPL(gO2,gI2))*
         CpbarUFuHpFdPL(gO1,gI2);
   }
   tmp_486 += tmp_487;
   result += (-0.5) * tmp_486;
   std::complex<double> tmp_488;
   std::complex<double> tmp_489;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_489 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_488 += tmp_489;
   result += (-1) * tmp_488;
   std::complex<double> tmp_490;
   std::complex<double> tmp_491;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_491 += (0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPR(gO2,gI2))*
         CpbarUFuVWpFdPR(gO1,gI2);
   }
   tmp_490 += tmp_491;
   result += (-1) * tmp_490;
   std::complex<double> tmp_492;
   std::complex<double> tmp_493;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_493 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_492 += tmp_493;
   result += (-1) * tmp_492;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::tadpole_hh() const
{
   std::complex<double> result;

   result += -0.5*A0(MAh)*CphhAhAh();
   result += A0(MVWp)*CphhbargWpCgWpC();
   result += A0(MVWp)*CphhbargWpgWp();
   result += A0(MVZ)*CphhbargZgZ();
   result += -(A0(MHp)*CphhconjHpHp());
   result += -0.5*A0(Mhh)*Cphhhhhh();
   std::complex<double> tmp_494;
   std::complex<double> tmp_495;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_495 += A0(MFd(gI1))*(CphhbarFdFdPL(gI1,gI1) + CphhbarFdFdPR(gI1,
         gI1))*MFd(gI1);
   }
   tmp_494 += tmp_495;
   result += (6) * tmp_494;
   std::complex<double> tmp_496;
   std::complex<double> tmp_497;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_497 += A0(MFe(gI1))*(CphhbarFeFePL(gI1,gI1) + CphhbarFeFePR(gI1,
         gI1))*MFe(gI1);
   }
   tmp_496 += tmp_497;
   result += (2) * tmp_496;
   std::complex<double> tmp_498;
   std::complex<double> tmp_499;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_499 += A0(MFu(gI1))*(CphhbarFuFuPL(gI1,gI1) + CphhbarFuFuPR(gI1,
         gI1))*MFu(gI1);
   }
   tmp_498 += tmp_499;
   result += (6) * tmp_498;
   result += 4*CphhconjVWpVWp()*(A0(MVWp) - 0.5*Sqr(MVWp));
   result += 2*CphhVZVZ()*(A0(MVZ) - 0.5*Sqr(MVZ));

   return result * oneOver16PiSqr;

}










void CLASSNAME::calculate_MVG_pole()
{
   // diagonalization with low precision
   PHYSICAL(MVG) = 0.;
}

void CLASSNAME::calculate_MFv_pole()
{
   // diagonalization with low precision
   PHYSICAL(MFv).setConstant(0.);
}

void CLASSNAME::calculate_MFc_pole()
{
   // diagonalization with low precision
   const double M_tree(MFc);
   const double p = MFc;
   const double self_energy_1  = Re(self_energy_Fc_1(p));
   const double self_energy_PL = Re(self_energy_Fc_PL(p));
   const double self_energy_PR = Re(self_energy_Fc_PR(p));
   const auto M_loop = M_tree - self_energy_1 - M_tree * (self_energy_PL
      + self_energy_PR);
   PHYSICAL(MFc) = calculate_singlet_mass(M_loop);
}

void CLASSNAME::calculate_MFg_pole()
{
   // diagonalization with low precision
   const double M_tree(MFg);
   const double p = MFg;
   const double self_energy_1  = Re(self_energy_Fg_1(p));
   const double self_energy_PL = Re(self_energy_Fg_PL(p));
   const double self_energy_PR = Re(self_energy_Fg_PR(p));
   const auto M_loop = M_tree - self_energy_1 - M_tree * (self_energy_PL
      + self_energy_PR);
   PHYSICAL(MFg) = calculate_singlet_mass(M_loop);
}

void CLASSNAME::calculate_MFn_pole()
{
   // diagonalization with low precision
   const double M_tree(MFn);
   const double p = MFn;
   const double self_energy_1  = Re(self_energy_Fn_1(p));
   const double self_energy_PL = Re(self_energy_Fn_PL(p));
   const double self_energy_PR = Re(self_energy_Fn_PR(p));
   const auto M_loop = M_tree - self_energy_1 - M_tree * (self_energy_PL
      + self_energy_PR);
   PHYSICAL(MFn) = calculate_singlet_mass(M_loop);
}

void CLASSNAME::calculate_Mhh_pole()
{
   if (!force_output && problems.is_tachyon(MDM_info::hh))
      return;

   // diagonalization with low precision
   const double M_tree(get_mass_matrix_hh());
   const double p = Mhh;
   double self_energy = Re(self_energy_hh(p));
   const double mass_sqr = M_tree - self_energy;

   PHYSICAL(Mhh) = SignedAbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MVP_pole()
{
   // diagonalization with low precision
   PHYSICAL(MVP) = 0.;
}

void CLASSNAME::calculate_MVZ_pole()
{
   if (!force_output && problems.is_tachyon(MDM_info::VZ))
      return;

   // diagonalization with low precision
   const double M_tree(Sqr(MVZ));
   const double p = MVZ;
   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(MDM_info::VZ);

   PHYSICAL(MVZ) = AbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MFd_pole()
{
   // diagonalization with low precision
   Eigen::Array<double,3,1> MFd_reordered(MFd);
   reorder_vector(MFd_reordered, get_mass_matrix_Fd());

   Eigen::Matrix<double,3,3> self_energy_1;
   Eigen::Matrix<double,3,3> self_energy_PL;
   Eigen::Matrix<double,3,3> self_energy_PR;
   for (unsigned i1 = 0; i1 < 3; ++i1) {
      for (unsigned i2 = 0; i2 < 3; ++i2) {
         const double p = AbsSqrt(MFd_reordered(i1) * MFd_reordered
            (i2));
         self_energy_1(i1,i2)  = Re(self_energy_Fd_1(p,i1,i2));
         self_energy_PL(i1,i2) = Re(self_energy_Fd_PL(p,i1,i2));
         self_energy_PR(i1,i2) = Re(self_energy_Fd_PR(p,i1,i2));
      }
   }
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fd());
   const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
      M_tree * self_energy_PL - self_energy_1);
   const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, PHYSICAL(MFd), PHYSICAL(Vd), PHYSICAL(Ud),
         eigenvalue_error);
      problems.flag_bad_mass(MDM_info::Fd, eigenvalue_error >
         precision * Abs(PHYSICAL(MFd)(0)));
   #else
      fs_svd(M_loop, PHYSICAL(MFd), PHYSICAL(Vd), PHYSICAL(Ud));
   #endif
}

void CLASSNAME::calculate_MFu_pole()
{
   // diagonalization with low precision
   Eigen::Array<double,3,1> MFu_reordered(MFu);
   reorder_vector(MFu_reordered, get_mass_matrix_Fu());

   Eigen::Matrix<double,3,3> self_energy_1;
   Eigen::Matrix<double,3,3> self_energy_PL;
   Eigen::Matrix<double,3,3> self_energy_PR;
   for (unsigned i1 = 0; i1 < 3; ++i1) {
      for (unsigned i2 = 0; i2 < 3; ++i2) {
         const double p = AbsSqrt(MFu_reordered(i1) * MFu_reordered
            (i2));
         self_energy_1(i1,i2)  = Re(self_energy_Fu_1(p,i1,i2));
         self_energy_PL(i1,i2) = Re(self_energy_Fu_PL(p,i1,i2));
         self_energy_PR(i1,i2) = Re(self_energy_Fu_PR(p,i1,i2));
      }
   }
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fu());
   const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
      M_tree * self_energy_PL - self_energy_1);
   const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, PHYSICAL(MFu), PHYSICAL(Vu), PHYSICAL(Uu),
         eigenvalue_error);
      problems.flag_bad_mass(MDM_info::Fu, eigenvalue_error >
         precision * Abs(PHYSICAL(MFu)(0)));
   #else
      fs_svd(M_loop, PHYSICAL(MFu), PHYSICAL(Vu), PHYSICAL(Uu));
   #endif
}

void CLASSNAME::calculate_MFe_pole()
{
   // diagonalization with low precision
   Eigen::Array<double,3,1> MFe_reordered(MFe);
   reorder_vector(MFe_reordered, get_mass_matrix_Fe());

   Eigen::Matrix<double,3,3> self_energy_1;
   Eigen::Matrix<double,3,3> self_energy_PL;
   Eigen::Matrix<double,3,3> self_energy_PR;
   for (unsigned i1 = 0; i1 < 3; ++i1) {
      for (unsigned i2 = 0; i2 < 3; ++i2) {
         const double p = AbsSqrt(MFe_reordered(i1) * MFe_reordered
            (i2));
         self_energy_1(i1,i2)  = Re(self_energy_Fe_1(p,i1,i2));
         self_energy_PL(i1,i2) = Re(self_energy_Fe_PL(p,i1,i2));
         self_energy_PR(i1,i2) = Re(self_energy_Fe_PR(p,i1,i2));
      }
   }
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fe());
   const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
      M_tree * self_energy_PL - self_energy_1);
   const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, PHYSICAL(MFe), PHYSICAL(Ve), PHYSICAL(Ue),
         eigenvalue_error);
      problems.flag_bad_mass(MDM_info::Fe, eigenvalue_error >
         precision * Abs(PHYSICAL(MFe)(0)));
   #else
      fs_svd(M_loop, PHYSICAL(MFe), PHYSICAL(Ve), PHYSICAL(Ue));
   #endif
}

void CLASSNAME::calculate_MVWp_pole()
{
   if (!force_output && problems.is_tachyon(MDM_info::VWp))
      return;

   // diagonalization with low precision
   const double M_tree(Sqr(MVWp));
   const double p = MVWp;
   const double self_energy = Re(self_energy_VWp(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(MDM_info::VWp);

   PHYSICAL(MVWp) = AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWp_pole(double p)
{
   if (!force_output && problems.is_tachyon(MDM_info::VWp))
      return 0.;

   const double self_energy = Re(self_energy_VWp(p));
   const double mass_sqr = Sqr(MVWp) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(MDM_info::VWp);

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVZ_pole(double p)
{
   if (!force_output && problems.is_tachyon(MDM_info::VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = Sqr(MVZ) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(MDM_info::VZ);

   return AbsSqrt(mass_sqr);
}


double CLASSNAME::calculate_MFv_DRbar(double, int) const
{
   return 0.0;
}

double CLASSNAME::calculate_MFe_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fe_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fe_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fe_PR_heavy_rotated(p,
      idx, idx));
   const double drbar_conversion = 1;
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;

   const double m_susy_drbar = m_sm_drbar + self_energy_1 + m_sm_drbar *
      (self_energy_PL + self_energy_PR);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFu_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fu_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fu_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fu_PR_heavy_rotated(p,
      idx, idx));

   const double currentScale = get_scale();
   const double qcd_1l = -0.008443431970194815*(4. - 3.*Log(Sqr(MFu(idx))
      /Sqr(currentScale)))*Sqr(g3);
   double qcd_2l = 0., qcd_3l = 0.;

   if (get_thresholds() > 1) {
      qcd_2l = -0.0041441100714622115*Power(g3,4) -
         0.0015238567409297061*Power(g3,4)*Log(Sqr(currentScale)/Sqr(MFu(idx)))
         - 0.00024060895909416413*Power(g3,4)*Sqr(Log(Power(currentScale,2)
         /Sqr(MFu(idx))));
   }

   if (get_thresholds() > 2) {
      qcd_3l = -0.0008783313853540776*Power(g3,6) -
         0.0004114970933517977*Power(g3,6)*Log(Sqr(currentScale)/Sqr(MFu(idx)))
         - 5.078913443827405e-6*Power(g3,6)*Power(Log(Sqr(currentScale)/Sqr(
         MFu(idx))),3) - 0.0002952541682011665*Power(g3,6)*Log(Sqr(MFu(idx))
         /Sqr(currentScale)) + 0.00005282069981580501*Power(g3,6)*Sqr(Log(Power
         (MFu(idx),2)/Sqr(currentScale))) - 0.00007466002762426286*Power(g3,6)*
         Sqr(Log(Power(currentScale,2)/Sqr(MFu(idx))));
   }

   const double m_susy_drbar = m_pole + self_energy_1 + m_pole * (
      self_energy_PL + self_energy_PR + qcd_1l + qcd_2l + qcd_3l);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFd_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fd_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fd_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fd_PR_heavy_rotated(p,
      idx, idx));
   const double m_tree = MFd(idx);
   const double drbar_conversion = 1;
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;

   const double m_susy_drbar = m_sm_drbar / (1.0 - self_energy_1/m_tree -
      self_energy_PL - self_energy_PR);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MVP_DRbar(double)
{
   return 0.0;
}

double CLASSNAME::calculate_MVZ_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_tachyon(MDM_info::VZ);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWp_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VWp(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_tachyon(MDM_info::VWp);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}


double CLASSNAME::ThetaW() const
{
   return ArcCos(Abs(ZZ(0,0)));
}


std::ostream& operator<<(std::ostream& ostr, const MDM_mass_eigenstates& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
