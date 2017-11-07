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

// File generated at Tue 7 Nov 2017 22:04:07

#ifndef EW_triplet_SLHA_IO_H
#define EW_triplet_SLHA_IO_H

#include "EW_triplet_two_scale_model_slha.hpp"
#include "EW_triplet_info.hpp"
#include "EW_triplet_observables.hpp"
#include "EW_triplet_physical.hpp"
#include "slha_io.hpp"
#include "ckm.hpp"
#include "ew_input.hpp"
#include "lowe.h"

#include <Eigen/Core>
#include <string>
#include <utility>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define PHYSICAL_SLHA(p) model.get_physical_slha().p
#define LOCALPHYSICAL(p) physical.p
#define MODEL model
#define MODELPARAMETER(p) model.get_##p()
#define OBSERVABLES observables
#define LowEnergyConstant(p) Electroweak_constants::p
#define SCALES(p) scales.p

namespace flexiblesusy {

struct EW_triplet_input_parameters;
class Spectrum_generator_settings;

struct EW_triplet_scales {
   EW_triplet_scales() : HighScale(0.), SUSYScale(0.), LowScale(0.) {}
   double HighScale, SUSYScale, LowScale;
};

class EW_triplet_slha_io {
public:
   EW_triplet_slha_io();

   void clear();

   void fill(softsusy::QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(EW_triplet_input_parameters&) const;
   void fill(EW_triplet_mass_eigenstates&) const;
   template <class T> void fill(EW_triplet_slha<T>&) const;
   void fill(Physical_input&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void read_from_source(const std::string&);
   void read_from_stream(std::istream&);
   void set_block(const std::string& str, SLHA_io::Position position = SLHA_io::back) { slha_io.set_block(str, position); }
   void set_blocks(const std::vector<std::string>& vec, SLHA_io::Position position = SLHA_io::back) { slha_io.set_blocks(vec, position); }
   void set_extpar(const EW_triplet_input_parameters&);
   template <class T> void set_extra(const EW_triplet_slha<T>&, const EW_triplet_scales&, const EW_triplet_observables&);
   void set_minpar(const EW_triplet_input_parameters&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class T> void set_spectrum(const EW_triplet_slha<T>&);
   template <class T> void set_spectrum(const EW_triplet<T>&);
   void set_spinfo(const Problems<EW_triplet_info::NUMBER_OF_PARTICLES>&);
   void set_spinfo(const std::vector<std::string>&, const std::vector<std::string>&);
   void set_print_imaginary_parts_of_majorana_mixings(bool);
   void write_to(const std::string&);
   void write_to_file(const std::string& file_name) { slha_io.write_to_file(file_name); }
   void write_to_stream(std::ostream& ostr = std::cout) { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(EW_triplet_input_parameters&, int, double);
   static void fill_extpar_tuple(EW_triplet_input_parameters&, int, double);

   template <class T>
   static void fill_slhaea(SLHAea::Coll&, const EW_triplet_slha<T>&, const softsusy::QedQcd&, const EW_triplet_scales&, const EW_triplet_observables&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const EW_triplet_slha<T>&, const softsusy::QedQcd&, const EW_triplet_scales&, const EW_triplet_observables&);

private:
   SLHA_io slha_io; ///< SLHA io class
   bool print_imaginary_parts_of_majorana_mixings;
   static unsigned const NUMBER_OF_DRBAR_BLOCKS = 6;
   static char const * const drbar_blocks[NUMBER_OF_DRBAR_BLOCKS];

   void set_mass(const EW_triplet_physical&, bool);
   void set_mixing_matrices(const EW_triplet_physical&, bool);
   template <class T> void set_model_parameters(const EW_triplet_slha<T>&);
   void set_ckm(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   void set_pmns(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   double read_scale() const;
   void fill_drbar_parameters(EW_triplet_mass_eigenstates&) const;
   void fill_physical(EW_triplet_physical&) const;
};

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class T>
void EW_triplet_slha_io::fill(EW_triplet_slha<T>& model) const
{
   fill(static_cast<EW_triplet_mass_eigenstates&>(model));
   fill_physical(model.get_physical_slha());
}

template <class T>
void EW_triplet_slha_io::fill_slhaea(
   SLHAea::Coll& slhaea, const EW_triplet_slha<T>& model,
   const softsusy::QedQcd& qedqcd, const EW_triplet_scales& scales,
   const EW_triplet_observables& observables)
{
   EW_triplet_slha_io slha_io;
   const EW_triplet_input_parameters& input = model.get_input();
   const Problems<EW_triplet_info::NUMBER_OF_PARTICLES>& problems
      = model.get_problems();
   const bool error = problems.have_problem();

   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(qedqcd);
   slha_io.set_minpar(input);
   slha_io.set_extpar(input);
   if (!error) {
      slha_io.set_spectrum(model);
      slha_io.set_extra(model, scales, observables);
   }

   slhaea = slha_io.get_slha_io().get_data();
}

template <class T>
SLHAea::Coll EW_triplet_slha_io::fill_slhaea(
   const EW_triplet_slha<T>& model, const softsusy::QedQcd& qedqcd,
   const EW_triplet_scales& scales, const EW_triplet_observables& observables)
{
   SLHAea::Coll slhaea;
   EW_triplet_slha_io::fill_slhaea(slhaea, model, qedqcd, scales, observables);

   return slhaea;
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class T>
void EW_triplet_slha_io::set_model_parameters(const EW_triplet_slha<T>& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "gY")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("Yu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block("Yd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block("Ye", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   {
      std::ostringstream block;
      block << "Block SM Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(mu2)), "mu2")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(3, (MODELPARAMETER(v)), "v")
      ;
      slha_io.set_block(block);
   }


}

/**
 * Writes extra SLHA blocks
 *
 * @param model model class
 * @param scales struct of boundary condition scales
 * @param observables struct of observables
 */
template <class T>
void EW_triplet_slha_io::set_extra(
   const EW_triplet_slha<T>& model, const EW_triplet_scales& scales,
   const EW_triplet_observables& observables)
{
   const EW_triplet_physical physical(model.get_physical_slha());


}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void EW_triplet_slha_io::set_spectrum(const EW_triplet<T>& model)
{
   const EW_triplet_slha<T> model_slha(model);
   set_spectrum(model_slha);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class T>
void EW_triplet_slha_io::set_spectrum(const EW_triplet_slha<T>& model)
{
   const EW_triplet_physical physical(model.get_physical_slha());
   const bool write_sm_masses = model.do_calculate_sm_pole_masses();

   set_model_parameters(model);
   set_mass(physical, write_sm_masses);
   set_mixing_matrices(physical, write_sm_masses);

   if (slha_io.get_modsel().quark_flavour_violated)
      set_ckm(model.get_ckm_matrix(), model.get_scale());

   if (slha_io.get_modsel().lepton_flavour_violated)
      set_pmns(model.get_pmns_matrix(), model.get_scale());
}

} // namespace flexiblesusy

#undef Pole
#undef PHYSICAL
#undef PHYSICAL_SLHA
#undef LOCALPHYSICAL
#undef MODEL
#undef MODELPARAMETER
#undef OBSERVABLES
#undef LowEnergyConstant
#undef SCALES

#endif
