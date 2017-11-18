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

// File generated at Wed 15 Nov 2017 19:56:29

#include "EW_triplet_standard_model_matching.hpp"
#include "wrappers.hpp"
#include "two_scale_matching.hpp"
#include "two_loop_corrections.hpp"
#include "standard_model.hpp"
#include "EW_triplet_mass_eigenstates.hpp"
#include "EW_triplet_info.hpp"
#include "config.h"
#include "parallel.hpp"
#include <cmath>

using namespace flexiblesusy::standard_model;

namespace flexiblesusy {

#define MODELPARAMETER(p) model.get_##p()
#define SMPARAMETER(p) sm.get_##p()
#define INPUTPARAMETER(p) model.get_input().p

namespace {

template <class T>
class Loop_order_setter {
public:
   Loop_order_setter(T& model_, unsigned loop_order_)
      : model(model_), loop_order(loop_order_)
   {
      model_pole_mass_order = model.get_pole_mass_loop_order();
      model_ewsb_order = model.get_ewsb_loop_order();
      tlc = model.get_two_loop_corrections();

      // set top quark QCD corrections to given loop order
      Two_loop_corrections tlc_new = tlc;
      tlc_new.top_qcd = loop_order_ ? loop_order_ - 1 : 0;

      model.set_pole_mass_loop_order(loop_order);
      model.set_ewsb_loop_order(loop_order);
      model.set_two_loop_corrections(tlc_new);
   }

   ~Loop_order_setter() {
      model.set_pole_mass_loop_order(model_pole_mass_order);
      model.set_ewsb_loop_order(model_ewsb_order);
      model.set_two_loop_corrections(tlc);
      model.calculate_DRbar_masses();
      model.solve_ewsb();
   }
private:
   T& model;
   Two_loop_corrections tlc;       ///< 2-loop corrections
   unsigned loop_order;            ///< temporary loop order
   unsigned model_pole_mass_order; ///< old model pole mass loop order
   unsigned model_ewsb_order;      ///< old model EWSB loop order
};

} // anonymous namespace

/**
 * Calculates \f$\lambda(Q)\f$ at the tree level from the lightest
 * CP-even Higgs boson mass of the EW_triplet.
 *
 * @param sm Standard Model
 * @param model BSM model
 * @param idx Higgs index (in mass ordered Higgs multiplet)
 */
void EW_triplet_standard_model_matching::match_high_to_low_scale_model_tree_level(
   Standard_model& sm, EW_triplet_mass_eigenstates& model, unsigned idx)
{
   model.calculate_DRbar_masses();
   sm.set_Lambdax(Sqr(model.get_Mhh()/sm.get_v()));
   sm.calculate_DRbar_masses();
}


/**
 * Calculates \f$\lambda(Q)\f$ at the 1-loop level from the lightest
 * CP-even Higgs boson mass of the EW_triplet by requiring that the
 * 1-loop Higgs pole masses are equal in both models.
 *
 * @param sm Standard Model
 * @param model BSM model
 * @param loop_order downwards matching loop order
 * @param idx Higgs index (in mass ordered Higgs multiplet)
 */
void EW_triplet_standard_model_matching::match_high_to_low_scale_model(
   Standard_model& sm, EW_triplet_mass_eigenstates& model, unsigned loop_order, unsigned idx)
{
   if (loop_order == 0) {
      match_high_to_low_scale_model_tree_level(sm, model, idx);
      return;
   }

   // temporarily set loop order to `loop_order'
   Loop_order_setter<Standard_model> los_sm(sm, loop_order);
   Loop_order_setter<EW_triplet_mass_eigenstates> los_model(model, loop_order);

   match_high_to_low_scale_model(sm, model, idx);

   model.get_physical().clear();
   sm.get_physical().clear();
}


/**
 * Calculates \f$\lambda(Q)\f$ at the current loop level from the
 * lightest CP-even Higgs boson mass of the EW_triplet by requiring
 * that the Higgs pole masses are equal in both models.
 *
 * @param sm Standard Model
 * @param model BSM model
 * @param idx Higgs index (in mass ordered Higgs multiplet)
 */
void EW_triplet_standard_model_matching::match_high_to_low_scale_model(
   Standard_model& sm, EW_triplet_mass_eigenstates& model, unsigned idx)
{
   model.calculate_DRbar_masses();
   model.solve_ewsb();
   model.calculate_Mhh_pole();

   sm.calculate_DRbar_masses();
   sm.solve_ewsb();
   sm.calculate_Mhh_pole();

   sm.set_Lambdax((Sqr(model.get_physical().Mhh)
                   - Sqr(sm.get_physical().Mhh) + Sqr(sm.get_Mhh()))/Sqr(sm.get_v()));
}

/**
 * Calculates the gauge and Yukawa couplings and the SM-like VEV in
 * the EW_triplet at the tree level from the known Standard Model
 * couplings and the SM vev.
 */
void EW_triplet_standard_model_matching::match_low_to_high_scale_model_tree_level(
   EW_triplet_mass_eigenstates& model, Standard_model& sm)
{
   sm.calculate_DRbar_masses();
   model.set_g1(sm.get_g1()*standard_model_info::normalization_g1/EW_triplet_info::normalization_g1);
   model.set_g2(sm.get_g2()*standard_model_info::normalization_g2/EW_triplet_info::normalization_g2);
   model.set_g3(sm.get_g3()*standard_model_info::normalization_g3/EW_triplet_info::normalization_g3);

   {
      EW_triplet_mass_eigenstates* MODEL = &model;
      const double VEV = sm.get_v();


   }

   Eigen::Matrix<double, 3, 3> upQuarksDRbar    = ZEROMATRIX(3,3);
   Eigen::Matrix<double, 3, 3> downQuarksDRbar  = ZEROMATRIX(3,3);
   Eigen::Matrix<double, 3, 3> downLeptonsDRbar = ZEROMATRIX(3,3);

   for (unsigned gen = 0; gen < 3; gen++) {
      upQuarksDRbar(gen, gen)    = sm.get_MFu(gen);
      downQuarksDRbar(gen, gen)  = sm.get_MFd(gen);
      downLeptonsDRbar(gen, gen) = sm.get_MFe(gen);
   }



   model.calculate_DRbar_masses();
}

static void calculate_SM_pole_masses(EW_triplet_mass_eigenstates& model, Standard_model& sm)
{
#ifdef ENABLE_THREADS
   EW_triplet_mass_eigenstates* obj_ptr = &model;

   auto fut_Standard_model_MFu = run_async([&sm] () { sm.calculate_MFu_pole(); });
   auto fut_Standard_model_MFd = run_async([&sm] () { sm.calculate_MFd_pole(); });
   auto fut_Standard_model_MFe = run_async([&sm] () { sm.calculate_MFe_pole(); });
   auto fut_Standard_model_MVW = run_async([&sm] () { sm.calculate_MVWp_pole(); });
   auto fut_Standard_model_MVZ = run_async([&sm] () { sm.calculate_MVZ_pole(); });

   fut_Standard_model_MFu.get();
   fut_Standard_model_MFd.get();
   fut_Standard_model_MFe.get();
   fut_Standard_model_MVW.get();
   fut_Standard_model_MVZ.get();
#else

   sm.calculate_MFu_pole();
   sm.calculate_MFd_pole();
   sm.calculate_MFe_pole();
   sm.calculate_MVWp_pole();
   sm.calculate_MVZ_pole();
#endif
}

/**
 * Calculates the gauge and Yukawa couplings and the SM-like VEV in
 * the EW_triplet at the 1-loop level from the known Standard Model
 * couplings and the SM vev.
 */
void EW_triplet_standard_model_matching::match_low_to_high_scale_model(
   EW_triplet_mass_eigenstates& model, Standard_model& sm, unsigned loop_order)
{
   if (loop_order == 0) {
      match_low_to_high_scale_model_tree_level(model, sm);
      return;
   }

   // temporarily set loop order to `loop_order'
   Loop_order_setter<Standard_model> los_sm(sm, loop_order);
   Loop_order_setter<EW_triplet_mass_eigenstates> los_model(model, loop_order);

   match_low_to_high_scale_model(model, sm);

   model.get_physical().clear();
   sm.get_physical().clear();
}

/**
 * Calculates the gauge and Yukawa couplings and the SM-like VEV in
 * the EW_triplet at the current loop level from the known Standard
 * Model couplings and the SM vev.
 */
void EW_triplet_standard_model_matching::match_low_to_high_scale_model(
   EW_triplet_mass_eigenstates& model, Standard_model& sm)
{
   Eigen::Matrix<double, 3, 3> upQuarksDRbar    = ZEROMATRIX(3,3);
   Eigen::Matrix<double, 3, 3> downQuarksDRbar  = ZEROMATRIX(3,3);
   Eigen::Matrix<double, 3, 3> downLeptonsDRbar = ZEROMATRIX(3,3);

   model.calculate_DRbar_masses();
   model.solve_ewsb();

   sm.calculate_DRbar_masses();
   sm.solve_ewsb();

   const double alpha_em = Sqr(sm.get_g1() * sm.get_g2() * standard_model_info::normalization_g1 * standard_model_info::normalization_g2)
            /(4. * Pi * (Sqr(sm.get_g1()*standard_model_info::normalization_g1) + Sqr(sm.get_g2()*standard_model_info::normalization_g2)));
   const double alpha_s  = Sqr(sm.get_g3() * standard_model_info::normalization_g3)/(4. * Pi);
   const double currentScale = sm.get_scale();
   double delta_alpha_em = 0., delta_alpha_s = 0.;





   calculate_SM_pole_masses(model, sm);

   // running W, Z masses (via 1L matching)
   const double mW2_1L = Sqr(sm.get_physical().MVWp) - Sqr(model.get_physical().MVWp) + Sqr(model.get_MVWp());
   const double mZ2_1L = Sqr(sm.get_physical().MVZ) - Sqr(model.get_physical().MVZ) + Sqr(model.get_MVZ());

   if (mZ2_1L < 0.)
      model.get_problems().flag_tachyon(EW_triplet_info::VZ);
   else
      model.get_problems().unflag_tachyon(EW_triplet_info::VZ);

   if (mW2_1L < 0.)
      model.get_problems().flag_tachyon(EW_triplet_info::VWp);
   else
      model.get_problems().unflag_tachyon(EW_triplet_info::VWp);



   // define and apply 1L-matched gauge parameters
   const double g1_1L = AbsSqrt(4. * Pi * alpha_em * (1. + delta_alpha_em) * mZ2_1L / mW2_1L) / EW_triplet_info::normalization_g1;
   const double g2_1L = AbsSqrt(4. * Pi * alpha_em * (1. + delta_alpha_em) / (1. - mW2_1L/mZ2_1L)) / EW_triplet_info::normalization_g2;
   const double g3_1L = AbsSqrt(4. * Pi * alpha_s * (1. + delta_alpha_s)) / EW_triplet_info::normalization_g3;

   model.set_g1(g1_1L);
   model.set_g2(g2_1L);
   model.set_g3(g3_1L);

   {
      EW_triplet_mass_eigenstates* MODEL = &model;
      const double VEV = 2. * AbsSqrt(mZ2_1L/(Sqr(g1_1L*EW_triplet_info::normalization_g1) + Sqr(g2_1L*EW_triplet_info::normalization_g2)));


   }


}

} // namespace flexiblesusy
