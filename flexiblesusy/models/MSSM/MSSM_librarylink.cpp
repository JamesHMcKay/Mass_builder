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

// File generated at Sat 2 Sep 2017 19:01:06

#include "MSSM_info.hpp"
#include "MSSM_input_parameters.hpp"
#include "MSSM_observables.hpp"
#include "MSSM_physical.hpp"
#include "MSSM_spectrum_generator.hpp"
#include "MSSM_two_scale_model.hpp"
#include "MSSM_two_scale_model_slha.hpp"

#include "error.hpp"
#include "physical_input.hpp"
#include "spectrum_generator_settings.hpp"
#include "lowe.h"

#include "mathlink.h"
#include "mathlink_utils.hpp"
#include "WolframLibrary.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>

#define INPUTPARAMETER(p) data.input.p
#define MODELPARAMETER(p) model.get_##p()
#define PHYSICALPARAMETER(p) model.get_physical().p
#define OBSERVABLE(o) observables.o

namespace MSSM_librarylink {

using namespace flexiblesusy;

typedef Two_scale algorithm_type;
typedef mint Handle_id;

class Redirect_output {
public:
   explicit Redirect_output(MLINK link_)
      : link(link_)
      , buffer()
      , old_cout(std::cout.rdbuf(buffer.rdbuf()))
      , old_cerr(std::cerr.rdbuf(buffer.rdbuf()))
      {}

   ~Redirect_output() {
      std::cout.rdbuf(old_cout);
      std::cerr.rdbuf(old_cerr);
      flush();
   }

private:
   MLINK link;               ///< redirect to this link
   std::stringstream buffer; ///< buffer caching stdout
   std::streambuf* old_cout; ///< original stdout buffer
   std::streambuf* old_cerr; ///< original stderr buffer

   void flush() {
      std::string line;
      while (std::getline(buffer, line)) {
         MLPutFunction(link, "CompoundExpression", 2);
         MLPutFunction(link, "FSMSSMMessage", 1);
         MLPutString(link, line.c_str());
      }
   }
};

class EUnknownHandle : public Error {
public:
   explicit EUnknownHandle(Handle_id hid_) : hid(hid_) {}
   virtual ~EUnknownHandle() {}
   virtual std::string what() const {
      return "Unknown handle: " + ToString(hid);
   }
   Handle_id hid;
};

class ENotEnoughFreeHandles : public Error {
public:
   explicit ENotEnoughFreeHandles(std::size_t max_handles_)
      : max_handles(max_handles_) {}
   virtual ~ENotEnoughFreeHandles() {}
   virtual std::string what() const {
      return "Maximum number of open handles reached: "
         + ToString(max_handles) + ".  Please close some handles!";
   }
   std::size_t max_handles;
};

class EWrongNumberOfParameters : public Error {
public:
   EWrongNumberOfParameters(unsigned pars_, unsigned expected_)
      : pars(pars_), expected(expected_) {}
   virtual ~EWrongNumberOfParameters() {}
   virtual std::string what() const {
      return "Wrong number of arguments: " + ToString(pars)
         + ".  Expected: " + ToString(expected);
   }
   unsigned pars, expected;
};

class EInvalidSpectrum : public Error {
public:
   EInvalidSpectrum() {}
   virtual ~EInvalidSpectrum() {}
   virtual std::string what() const { return "Invalid spectrum"; }
};

struct Model_data {
   Model_data()
      : input()
      , physical_input()
      , qedqcd()
      , settings()
      , parameter_output_scale(0.)
      , model()
   {}

   MSSM_input_parameters input;     ///< model input parameters
   Physical_input physical_input;          ///< extra non-SLHA physical input
   softsusy::QedQcd qedqcd;                ///< SLHA physical input
   Spectrum_generator_settings settings;   ///< spectrum generator settings
   double parameter_output_scale;          ///< output scale for running parameters
   MSSM_slha<algorithm_type> model; ///< running parameters and pole masses
};

/// current handles
typedef std::map<Handle_id, Model_data> Handle_map;
Handle_map handles;

/******************************************************************/

Handle_id get_new_handle()
{
   static const std::size_t max_handles =
      static_cast<std::size_t>(std::exp2(8*sizeof(Handle_id)) - 1);

   if (handles.size() >= max_handles)
      throw ENotEnoughFreeHandles(handles.size());

   Handle_id hid = 0;

   while (handles.find(hid) != handles.end())
      hid++;

   return hid;
}

/******************************************************************/

Model_data find_data(Handle_id hid)
{
   const Handle_map::iterator handle = handles.find(hid);

   if (handle == handles.end())
      throw EUnknownHandle(hid);

   return handle->second;
}

/******************************************************************/

static Handle_id get_handle_from(MLINK link)
{
   Handle_id hid;
   MLGet(link, &hid);

   return hid;
}

/******************************************************************/

static long number_of_args(MLINK link, const std::string& head)
{
   long argc;

   if (!MLCheckFunction(link, head.c_str(), &argc))
      std::cerr << "Error: argument is not a " << head << std::endl;

   return argc;
}

/******************************************************************/

bool check_number_of_args(MLINK link, unsigned number_of_arguments,
                          const std::string& function_name)
{
   const long n_given = number_of_args(link, "List");
   const bool ok = n_given == number_of_arguments;

   if (!ok) {
      std::cerr << "Error: " << function_name << " expects "
                << ToString(number_of_arguments) << " argument ("
                << n_given << " given)." << std::endl;
   }

   return ok;
}

/******************************************************************/

static void put_error_output(MLINK link)
{
   MLPutSymbol(link, "$Failed");
}

/******************************************************************/

void put_message(MLINK link,
                 const std::string& function_name,
                 const std::string& message_tag,
                 const std::string& message_str)
{
   MLPutFunction(link, "CompoundExpression", 2);
   MLPutFunction(link, "Message", 2);
   MLPutFunction(link, "MessageName", 2);
   MLPutSymbol(link, function_name.c_str());
   MLPutString(link, message_tag.c_str());
   MLPutString(link, message_str.c_str());
}

/******************************************************************/

void put_settings(const Model_data& data, MLINK link)
{
   MLPutFunction(link, "List", Spectrum_generator_settings::NUMBER_OF_OPTIONS - 1);

   MLPutRuleTo(link, data.settings.get(Spectrum_generator_settings::precision), "precisionGoal");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::max_iterations), "maxIterations");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::calculate_sm_masses), "calculateStandardModelMasses");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::pole_mass_loop_order), "poleMassLoopOrder");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::ewsb_loop_order), "ewsbLoopOrder");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::beta_loop_order), "betaFunctionLoopOrder");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::threshold_corrections_loop_order), "thresholdCorrectionsLoopOrder");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::higgs_2loop_correction_at_as), "higgs2loopCorrectionAtAs");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::higgs_2loop_correction_ab_as), "higgs2loopCorrectionAbAs");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::higgs_2loop_correction_at_at), "higgs2loopCorrectionAtAt");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::higgs_2loop_correction_atau_atau), "higgs2loopCorrectionAtauAtau");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::force_output), "forceOutput");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::top_pole_qcd_corrections), "topPoleQCDCorrections");
   MLPutRuleTo(link, data.settings.get(Spectrum_generator_settings::beta_zero_threshold), "betaZeroThreshold");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::force_positive_masses), "forcePositiveMasses");
   MLPutRuleTo(link, data.settings.get(Spectrum_generator_settings::pole_mass_scale), "poleMassScale");
   MLPutRuleTo(link, data.settings.get(Spectrum_generator_settings::eft_pole_mass_scale), "eftPoleMassScale");
   MLPutRuleTo(link, data.settings.get(Spectrum_generator_settings::eft_matching_scale), "eftMatchingScale");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::eft_matching_loop_order_up), "eftMatchingLoopOrderUp");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::eft_matching_loop_order_down), "eftMatchingLoopOrderDown");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::eft_higgs_index), "eftHiggsIndex");
   MLPutRuleTo(link, (int)data.settings.get(Spectrum_generator_settings::calculate_bsm_masses), "calculateBSMMasses");
   MLPutRuleTo(link, data.parameter_output_scale, "parameterOutputScale");

   MLEndPacket(link);
}

/******************************************************************/

void put_sm_input_parameters(const Model_data& data, MLINK link)
{
   MLPutFunction(link, "List",softsusy::NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS
                              + Physical_input::NUMBER_OF_INPUT_PARAMETERS);

   MLPutRuleTo(link, data.qedqcd.displayAlphaEmInput(), "alphaEmMZ");
   MLPutRuleTo(link, data.qedqcd.displayFermiConstant(), "GF");
   MLPutRuleTo(link, data.qedqcd.displayAlphaSInput(), "alphaSMZ");
   MLPutRuleTo(link, data.qedqcd.displayPoleMZ(), "MZ");
   MLPutRuleTo(link, data.qedqcd.displayMbMb(), "mbmb");
   MLPutRuleTo(link, data.qedqcd.displayPoleMt(), "Mt");
   MLPutRuleTo(link, data.qedqcd.displayPoleMtau(), "Mtau");
   MLPutRuleTo(link, data.qedqcd.displayNeutrinoPoleMass(3), "Mv3");
   MLPutRuleTo(link, data.qedqcd.displayPoleMW(), "MW");
   MLPutRuleTo(link, data.qedqcd.displayPoleMel(), "Me");
   MLPutRuleTo(link, data.qedqcd.displayNeutrinoPoleMass(1), "Mv1");
   MLPutRuleTo(link, data.qedqcd.displayPoleMmuon(), "Mm");
   MLPutRuleTo(link, data.qedqcd.displayNeutrinoPoleMass(2), "Mv2");
   MLPutRuleTo(link, data.qedqcd.displayMd2GeV(), "md2GeV");
   MLPutRuleTo(link, data.qedqcd.displayMu2GeV(), "mu2GeV");
   MLPutRuleTo(link, data.qedqcd.displayMs2GeV(), "ms2GeV");
   MLPutRuleTo(link, data.qedqcd.displayMcMc(), "mcmc");

   const flexiblesusy::CKM_parameters ckm(data.qedqcd.displayCKM());
   MLPutRuleTo(link, ckm.theta_12, "CKMTheta12");
   MLPutRuleTo(link, ckm.theta_13, "CKMTheta13");
   MLPutRuleTo(link, ckm.theta_23, "CKMTheta23");
   MLPutRuleTo(link, ckm.delta   , "CKMDelta");

   const flexiblesusy::PMNS_parameters pmns(data.qedqcd.displayPMNS());
   MLPutRuleTo(link, pmns.theta_12, "PMNSTheta12");
   MLPutRuleTo(link, pmns.theta_13, "PMNSTheta13");
   MLPutRuleTo(link, pmns.theta_23, "PMNSTheta23");
   MLPutRuleTo(link, pmns.delta   , "PMNSDelta");
   MLPutRuleTo(link, pmns.alpha_1 , "PMNSAlpha1");
   MLPutRuleTo(link, pmns.alpha_2 , "PMNSAlpha2");

   MLPutRuleTo(link, data.physical_input.get(Physical_input::alpha_em_0), "alphaEm0");
   MLPutRuleTo(link, data.physical_input.get(Physical_input::mh_pole), "Mh");

   MLEndPacket(link);
}

/******************************************************************/

void put_input_parameters(const Model_data& data, MLINK link)
{
   MLPutFunction(link, "List", 17);

   MLPutRuleTo(link, INPUTPARAMETER(TanBeta), "TanBeta");
   MLPutRuleTo(link, INPUTPARAMETER(SignMu), "SignMu");
   MLPutRuleTo(link, INPUTPARAMETER(MassWBIN), "MassWBIN");
   MLPutRuleTo(link, INPUTPARAMETER(Qin), "Qin");
   MLPutRuleTo(link, INPUTPARAMETER(QEWSB), "QEWSB");
   MLPutRuleTo(link, INPUTPARAMETER(mHd2IN), "mHd2IN");
   MLPutRuleTo(link, INPUTPARAMETER(mHu2IN), "mHu2IN");
   MLPutRuleTo(link, INPUTPARAMETER(Aeij), "Aeij");
   MLPutRuleTo(link, INPUTPARAMETER(Adij), "Adij");
   MLPutRuleTo(link, INPUTPARAMETER(Auij), "Auij");
   MLPutRuleTo(link, INPUTPARAMETER(mq2Input), "mq2Input");
   MLPutRuleTo(link, INPUTPARAMETER(ml2Input), "ml2Input");
   MLPutRuleTo(link, INPUTPARAMETER(md2Input), "md2Input");
   MLPutRuleTo(link, INPUTPARAMETER(mu2Input), "mu2Input");
   MLPutRuleTo(link, INPUTPARAMETER(me2Input), "me2Input");
   MLPutRuleTo(link, INPUTPARAMETER(MassBInput), "MassBInput");
   MLPutRuleTo(link, INPUTPARAMETER(MassGInput), "MassGInput");


   MLEndPacket(link);
}

/******************************************************************/

void put_spectrum(const MSSM_slha<algorithm_type>& model, MLINK link)
{
   MLPutFunction(link, "List", 94);

   MLPutRuleTo(link, MODELPARAMETER(MVG), "VG", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MGlu), "Glu", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MFv), "Fv", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MSd), "Sd", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MSv), "Sv", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MSu), "Su", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MSe), "Se", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(Mhh), "hh", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MAh), "Ah", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MHpm), "Hpm", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MChi), "Chi", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MCha), "Cha", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MFe), "Fe", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MFd), "Fd", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MFu), "Fu", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MVWm), "VWm", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MVP), "VP", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MVZ), "VZ", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(ZD), "ZD");
   MLPutRuleTo(link, MODELPARAMETER(ZV), "ZV");
   MLPutRuleTo(link, MODELPARAMETER(ZU), "ZU");
   MLPutRuleTo(link, MODELPARAMETER(ZE), "ZE");
   MLPutRuleTo(link, MODELPARAMETER(ZH), "ZH");
   MLPutRuleTo(link, MODELPARAMETER(ZA), "ZA");
   MLPutRuleTo(link, MODELPARAMETER(ZP), "ZP");
   MLPutRuleTo(link, MODELPARAMETER(ZN), "ZN");
   MLPutRuleTo(link, MODELPARAMETER(UM), "UM");
   MLPutRuleTo(link, MODELPARAMETER(UP), "UP");
   MLPutRuleTo(link, MODELPARAMETER(ZEL), "ZEL");
   MLPutRuleTo(link, MODELPARAMETER(ZER), "ZER");
   MLPutRuleTo(link, MODELPARAMETER(ZDL), "ZDL");
   MLPutRuleTo(link, MODELPARAMETER(ZDR), "ZDR");
   MLPutRuleTo(link, MODELPARAMETER(ZUL), "ZUL");
   MLPutRuleTo(link, MODELPARAMETER(ZUR), "ZUR");
   MLPutRuleTo(link, MODELPARAMETER(ZZ), "ZZ");
   MLPutRuleTo(link, PHYSICALPARAMETER(MVG), "VG", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MGlu), "Glu", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MFv), "Fv", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MSd), "Sd", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MSv), "Sv", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MSu), "Su", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MSe), "Se", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Mhh), "hh", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MAh), "Ah", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MHpm), "Hpm", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MChi), "Chi", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MCha), "Cha", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MFe), "Fe", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MFd), "Fd", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MFu), "Fu", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MVWm), "VWm", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MVP), "VP", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MVZ), "VZ", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZD), "ZD", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZV), "ZV", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZU), "ZU", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZE), "ZE", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZH), "ZH", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZA), "ZA", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZP), "ZP", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZN), "ZN", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(UM), "UM", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(UP), "UP", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZEL), "ZEL", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZER), "ZER", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZDL), "ZDL", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZDR), "ZDR", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZUL), "ZUL", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZUR), "ZUR", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZZ), "ZZ", {"Pole"});
   MLPutRuleTo(link, MODELPARAMETER(Yd), "Yd");
   MLPutRuleTo(link, MODELPARAMETER(Ye), "Ye");
   MLPutRuleTo(link, MODELPARAMETER(Yu), "Yu");
   MLPutRuleTo(link, MODELPARAMETER(Mu), "\u03bc");
   MLPutRuleTo(link, MODELPARAMETER(g1), "g1");
   MLPutRuleTo(link, MODELPARAMETER(g2), "g2");
   MLPutRuleTo(link, MODELPARAMETER(g3), "g3");
   MLPutRuleTo(link, MODELPARAMETER(vd), "vd");
   MLPutRuleTo(link, MODELPARAMETER(vu), "vu");
   MLPutRuleTo(link, MODELPARAMETER(TYd), "Yd", {"T"});
   MLPutRuleTo(link, MODELPARAMETER(TYe), "Ye", {"T"});
   MLPutRuleTo(link, MODELPARAMETER(TYu), "Yu", {"T"});
   MLPutRuleTo(link, MODELPARAMETER(BMu), "\u03bc", {"B"});
   MLPutRuleTo(link, MODELPARAMETER(mq2), "mq2");
   MLPutRuleTo(link, MODELPARAMETER(ml2), "ml2");
   MLPutRuleTo(link, MODELPARAMETER(mHd2), "mHd2");
   MLPutRuleTo(link, MODELPARAMETER(mHu2), "mHu2");
   MLPutRuleTo(link, MODELPARAMETER(md2), "md2");
   MLPutRuleTo(link, MODELPARAMETER(mu2), "mu2");
   MLPutRuleTo(link, MODELPARAMETER(me2), "me2");
   MLPutRuleTo(link, MODELPARAMETER(MassB), "MassB");
   MLPutRuleTo(link, MODELPARAMETER(MassWB), "MassWB");
   MLPutRuleTo(link, MODELPARAMETER(MassG), "MassG");
   MLPutRuleTo(link, MODELPARAMETER(scale), "SCALE");


   MLEndPacket(link);
}

/******************************************************************/

void put_observables(const MSSM_observables& observables, MLINK link)
{
   MLPutFunction(link, "List", 0);



   MLEndPacket(link);
}

/******************************************************************/

void check_spectrum(const Model_data& data, MLINK link)
{
   const Problems<MSSM_info::NUMBER_OF_PARTICLES>& problems
      = data.model.get_problems();

   if (problems.have_problem()) {
      std::ostringstream msg;
      problems.print_problems(msg);
      put_message(link, "FSMSSMCalculateSpectrum", "error", msg.str());
   }

   if (problems.have_warning()) {
      std::ostringstream msg;
      problems.print_warnings(msg);
      put_message(link, "FSMSSMCalculateSpectrum", "warning", msg.str());
   }

   if (problems.have_problem() &&
       !data.settings.get(Spectrum_generator_settings::force_output))
      throw EInvalidSpectrum();
}

/******************************************************************/

void calculate_spectrum(Model_data& data, MLINK link)
{
   softsusy::QedQcd qedqcd(data.qedqcd);

   try {
      qedqcd.to(qedqcd.displayPoleMZ());
   } catch (const flexiblesusy::Error& e) {
      put_message(link, "FSMSSMCalculateSpectrum", "error", e.what());
      throw EInvalidSpectrum();
   }

   MSSM_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_settings(data.settings);
   spectrum_generator.set_parameter_output_scale(data.parameter_output_scale);
   spectrum_generator.run(qedqcd, data.input);

   data.model = MSSM_slha<algorithm_type>(
      spectrum_generator.get_model(),
      data.settings.get(Spectrum_generator_settings::force_positive_masses) == 0.);
}

/******************************************************************/

Model_data make_data(double* pars, mint npars)
{
   Model_data data;

   const mint n_settings = Spectrum_generator_settings::NUMBER_OF_OPTIONS - 1,
      n_sm_parameters = softsusy::NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS
                        + Physical_input::NUMBER_OF_INPUT_PARAMETERS,
      n_input_pars = 81;
   const mint n_total = n_settings + n_sm_parameters + n_input_pars;

   if (npars != n_total)
      throw EWrongNumberOfParameters(npars, n_total);

   mint c = 0; // counter

   data.settings.set(Spectrum_generator_settings::precision, pars[c++]);
   data.settings.set(Spectrum_generator_settings::max_iterations, pars[c++]);
   data.settings.set(Spectrum_generator_settings::calculate_sm_masses, pars[c++]);
   data.settings.set(Spectrum_generator_settings::pole_mass_loop_order, pars[c++]);
   data.settings.set(Spectrum_generator_settings::ewsb_loop_order, pars[c++]);
   data.settings.set(Spectrum_generator_settings::beta_loop_order, pars[c++]);
   data.settings.set(Spectrum_generator_settings::threshold_corrections_loop_order, pars[c++]);
   data.settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, pars[c++]);
   data.settings.set(Spectrum_generator_settings::higgs_2loop_correction_ab_as, pars[c++]);
   data.settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, pars[c++]);
   data.settings.set(Spectrum_generator_settings::higgs_2loop_correction_atau_atau, pars[c++]);
   data.settings.set(Spectrum_generator_settings::force_output, pars[c++]);
   data.settings.set(Spectrum_generator_settings::top_pole_qcd_corrections, pars[c++]);
   data.settings.set(Spectrum_generator_settings::beta_zero_threshold, pars[c++]);
   data.settings.set(Spectrum_generator_settings::force_positive_masses, pars[c++]);
   data.settings.set(Spectrum_generator_settings::pole_mass_scale, pars[c++]);
   data.settings.set(Spectrum_generator_settings::eft_pole_mass_scale, pars[c++]);
   data.settings.set(Spectrum_generator_settings::eft_matching_scale, pars[c++]);
   data.settings.set(Spectrum_generator_settings::eft_matching_loop_order_up, pars[c++]);
   data.settings.set(Spectrum_generator_settings::eft_matching_loop_order_down, pars[c++]);
   data.settings.set(Spectrum_generator_settings::eft_higgs_index, pars[c++]);
   data.settings.set(Spectrum_generator_settings::calculate_bsm_masses, pars[c++]);
   data.parameter_output_scale = pars[c++];

   data.qedqcd.setAlpha(softsusy::ALPHA, pars[c]);
   data.qedqcd.setAlphaEmInput(pars[c++]);
   data.qedqcd.setFermiConstant(pars[c++]);
   data.qedqcd.setAlpha(softsusy::ALPHAS, pars[c]);
   data.qedqcd.setAlphaSInput(pars[c++]);
   data.qedqcd.setPoleMZ(pars[c]);
   data.qedqcd.setMu(pars[c++]);
   data.qedqcd.setMass(softsusy::mBottom, pars[c]);
   data.qedqcd.setMbMb(pars[c++]);
   data.qedqcd.setPoleMt(pars[c++]);
   data.qedqcd.setMass(softsusy::mTau, pars[c]);
   data.qedqcd.setPoleMtau(pars[c++]);
   data.qedqcd.setNeutrinoPoleMass(3, pars[c++]);
   data.qedqcd.setPoleMW(pars[c++]);
   data.qedqcd.setMass(softsusy::mElectron, pars[c]);
   data.qedqcd.setPoleMel(pars[c++]);
   data.qedqcd.setNeutrinoPoleMass(1, pars[c++]);
   data.qedqcd.setMass(softsusy::mMuon, pars[c]);
   data.qedqcd.setPoleMmuon(pars[c++]);
   data.qedqcd.setNeutrinoPoleMass(2, pars[c++]);
   data.qedqcd.setMass(softsusy::mDown, pars[c]);
   data.qedqcd.setMd2GeV(pars[c++]);
   data.qedqcd.setMass(softsusy::mUp, pars[c]);
   data.qedqcd.setMu2GeV(pars[c++]);
   data.qedqcd.setMass(softsusy::mStrange, pars[c]);
   data.qedqcd.setMs2GeV(pars[c++]);
   data.qedqcd.setMass(softsusy::mCharm, pars[c]);
   data.qedqcd.setMcMc(pars[c++]);

   {
      flexiblesusy::CKM_parameters ckm;
      ckm.theta_12 = pars[c++];
      ckm.theta_13 = pars[c++];
      ckm.theta_23 = pars[c++];
      ckm.delta    = pars[c++];
      data.qedqcd.setCKM(ckm);
   }

   {
      flexiblesusy::PMNS_parameters pmns;
      pmns.theta_12 = pars[c++];
      pmns.theta_13 = pars[c++];
      pmns.theta_23 = pars[c++];
      pmns.delta    = pars[c++];
      pmns.alpha_1  = pars[c++];
      pmns.alpha_2  = pars[c++];
      data.qedqcd.setPMNS(pmns);
   }

   data.physical_input.set(Physical_input::alpha_em_0, pars[c++]);
   data.physical_input.set(Physical_input::mh_pole, pars[c++]);

   INPUTPARAMETER(TanBeta) = pars[c++];
   INPUTPARAMETER(SignMu) = pars[c++];
   INPUTPARAMETER(MassWBIN) = pars[c++];
   INPUTPARAMETER(Qin) = pars[c++];
   INPUTPARAMETER(QEWSB) = pars[c++];
   INPUTPARAMETER(mHd2IN) = pars[c++];
   INPUTPARAMETER(mHu2IN) = pars[c++];
   INPUTPARAMETER(Aeij(0,0)) = pars[c++];
   INPUTPARAMETER(Aeij(0,1)) = pars[c++];
   INPUTPARAMETER(Aeij(0,2)) = pars[c++];
   INPUTPARAMETER(Aeij(1,0)) = pars[c++];
   INPUTPARAMETER(Aeij(1,1)) = pars[c++];
   INPUTPARAMETER(Aeij(1,2)) = pars[c++];
   INPUTPARAMETER(Aeij(2,0)) = pars[c++];
   INPUTPARAMETER(Aeij(2,1)) = pars[c++];
   INPUTPARAMETER(Aeij(2,2)) = pars[c++];
   INPUTPARAMETER(Adij(0,0)) = pars[c++];
   INPUTPARAMETER(Adij(0,1)) = pars[c++];
   INPUTPARAMETER(Adij(0,2)) = pars[c++];
   INPUTPARAMETER(Adij(1,0)) = pars[c++];
   INPUTPARAMETER(Adij(1,1)) = pars[c++];
   INPUTPARAMETER(Adij(1,2)) = pars[c++];
   INPUTPARAMETER(Adij(2,0)) = pars[c++];
   INPUTPARAMETER(Adij(2,1)) = pars[c++];
   INPUTPARAMETER(Adij(2,2)) = pars[c++];
   INPUTPARAMETER(Auij(0,0)) = pars[c++];
   INPUTPARAMETER(Auij(0,1)) = pars[c++];
   INPUTPARAMETER(Auij(0,2)) = pars[c++];
   INPUTPARAMETER(Auij(1,0)) = pars[c++];
   INPUTPARAMETER(Auij(1,1)) = pars[c++];
   INPUTPARAMETER(Auij(1,2)) = pars[c++];
   INPUTPARAMETER(Auij(2,0)) = pars[c++];
   INPUTPARAMETER(Auij(2,1)) = pars[c++];
   INPUTPARAMETER(Auij(2,2)) = pars[c++];
   INPUTPARAMETER(mq2Input(0,0)) = pars[c++];
   INPUTPARAMETER(mq2Input(0,1)) = pars[c++];
   INPUTPARAMETER(mq2Input(0,2)) = pars[c++];
   INPUTPARAMETER(mq2Input(1,0)) = pars[c++];
   INPUTPARAMETER(mq2Input(1,1)) = pars[c++];
   INPUTPARAMETER(mq2Input(1,2)) = pars[c++];
   INPUTPARAMETER(mq2Input(2,0)) = pars[c++];
   INPUTPARAMETER(mq2Input(2,1)) = pars[c++];
   INPUTPARAMETER(mq2Input(2,2)) = pars[c++];
   INPUTPARAMETER(ml2Input(0,0)) = pars[c++];
   INPUTPARAMETER(ml2Input(0,1)) = pars[c++];
   INPUTPARAMETER(ml2Input(0,2)) = pars[c++];
   INPUTPARAMETER(ml2Input(1,0)) = pars[c++];
   INPUTPARAMETER(ml2Input(1,1)) = pars[c++];
   INPUTPARAMETER(ml2Input(1,2)) = pars[c++];
   INPUTPARAMETER(ml2Input(2,0)) = pars[c++];
   INPUTPARAMETER(ml2Input(2,1)) = pars[c++];
   INPUTPARAMETER(ml2Input(2,2)) = pars[c++];
   INPUTPARAMETER(md2Input(0,0)) = pars[c++];
   INPUTPARAMETER(md2Input(0,1)) = pars[c++];
   INPUTPARAMETER(md2Input(0,2)) = pars[c++];
   INPUTPARAMETER(md2Input(1,0)) = pars[c++];
   INPUTPARAMETER(md2Input(1,1)) = pars[c++];
   INPUTPARAMETER(md2Input(1,2)) = pars[c++];
   INPUTPARAMETER(md2Input(2,0)) = pars[c++];
   INPUTPARAMETER(md2Input(2,1)) = pars[c++];
   INPUTPARAMETER(md2Input(2,2)) = pars[c++];
   INPUTPARAMETER(mu2Input(0,0)) = pars[c++];
   INPUTPARAMETER(mu2Input(0,1)) = pars[c++];
   INPUTPARAMETER(mu2Input(0,2)) = pars[c++];
   INPUTPARAMETER(mu2Input(1,0)) = pars[c++];
   INPUTPARAMETER(mu2Input(1,1)) = pars[c++];
   INPUTPARAMETER(mu2Input(1,2)) = pars[c++];
   INPUTPARAMETER(mu2Input(2,0)) = pars[c++];
   INPUTPARAMETER(mu2Input(2,1)) = pars[c++];
   INPUTPARAMETER(mu2Input(2,2)) = pars[c++];
   INPUTPARAMETER(me2Input(0,0)) = pars[c++];
   INPUTPARAMETER(me2Input(0,1)) = pars[c++];
   INPUTPARAMETER(me2Input(0,2)) = pars[c++];
   INPUTPARAMETER(me2Input(1,0)) = pars[c++];
   INPUTPARAMETER(me2Input(1,1)) = pars[c++];
   INPUTPARAMETER(me2Input(1,2)) = pars[c++];
   INPUTPARAMETER(me2Input(2,0)) = pars[c++];
   INPUTPARAMETER(me2Input(2,1)) = pars[c++];
   INPUTPARAMETER(me2Input(2,2)) = pars[c++];
   INPUTPARAMETER(MassBInput) = pars[c++];
   INPUTPARAMETER(MassGInput) = pars[c++];


   if (npars != c)
      throw EWrongNumberOfParameters(npars, c);

   return data;
}

} // namespace MSSM_librarylink

extern "C" {

/******************************************************************/

DLLEXPORT mint WolframLibrary_getVersion()
{
   return WolframLibraryVersion;
}

/******************************************************************/

DLLEXPORT int WolframLibrary_initialize(WolframLibraryData /* libData */)
{
   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSMSSMGetSettings(WolframLibraryData /* libData */, MLINK link)
{
   using namespace MSSM_librarylink;

   if (!check_number_of_args(link, 1, "FSMSSMGetSettings"))
      return LIBRARY_TYPE_ERROR;

   const Handle_id hid = get_handle_from(link);

   try {
      const Model_data data = find_data(hid);
      put_settings(data, link);
   } catch (const flexiblesusy::Error& e) {
      std::cerr << e.what() << std::endl;
      put_error_output(link);
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSMSSMGetSMInputParameters(WolframLibraryData /* libData */, MLINK link)
{
   using namespace MSSM_librarylink;

   if (!check_number_of_args(link, 1, "FSMSSMGetSMInputParameters"))
      return LIBRARY_TYPE_ERROR;

   const Handle_id hid = get_handle_from(link);

   try {
      const Model_data data = find_data(hid);
      put_sm_input_parameters(data, link);
   } catch (const flexiblesusy::Error& e) {
      std::cerr << e.what() << std::endl;
      put_error_output(link);
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSMSSMGetInputParameters(WolframLibraryData /* libData */, MLINK link)
{
   using namespace MSSM_librarylink;

   if (!check_number_of_args(link, 1, "FSMSSMGetInputParameters"))
      return LIBRARY_TYPE_ERROR;

   const Handle_id hid = get_handle_from(link);

   try {
      const Model_data data = find_data(hid);
      put_input_parameters(data, link);
   } catch (const flexiblesusy::Error& e) {
      std::cerr << e.what() << std::endl;
      put_error_output(link);
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSMSSMOpenHandle(
   WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
   using namespace MSSM_librarylink;

   if (Argc != 1)
      return LIBRARY_TYPE_ERROR;

   MTensor pars = MArgument_getMTensor(Args[0]);

   if (libData->MTensor_getType(pars) != MType_Real ||
       libData->MTensor_getRank(pars) != 1)
      return LIBRARY_TYPE_ERROR;

   try {
      Model_data data = make_data(
         libData->MTensor_getRealData(pars),
         libData->MTensor_getDimensions(pars)[0]);

      const Handle_id hid = get_new_handle();

      handles.insert(std::make_pair(hid, std::move(data)));

      MArgument_setInteger(Res, hid);
   } catch (const flexiblesusy::Error& e) {
      std::cerr << e.what() << std::endl;
      return LIBRARY_FUNCTION_ERROR;
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSMSSMCloseHandle(
   WolframLibraryData /* libData */, mint Argc, MArgument* Args, MArgument /* Res */)
{
   using namespace MSSM_librarylink;

   if (Argc != 1)
      return LIBRARY_TYPE_ERROR;

   const Handle_id hid = MArgument_getInteger(Args[0]);

   const Handle_map::iterator handle = handles.find(hid);

   if (handle != handles.end())
      handles.erase(handle);

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSMSSMSet(
   WolframLibraryData libData, mint Argc, MArgument* Args, MArgument /* Res */)
{
   using namespace MSSM_librarylink;

   if (Argc != 2)
      return LIBRARY_TYPE_ERROR;

   const Handle_id hid = MArgument_getInteger(Args[0]);
   MTensor pars = MArgument_getMTensor(Args[1]);

   if (libData->MTensor_getType(pars) != MType_Real ||
       libData->MTensor_getRank(pars) != 1)
      return LIBRARY_TYPE_ERROR;

   const Handle_map::iterator handle = handles.find(hid);

   if (handle == handles.end()) {
      std::cerr << "Error: FSMSSMSet: Unknown handle: "
                << hid << std::endl;
      return LIBRARY_FUNCTION_ERROR;
   }

   try {
      handle->second = make_data(
         libData->MTensor_getRealData(pars),
         libData->MTensor_getDimensions(pars)[0]);
   } catch (const flexiblesusy::Error& e) {
      std::cerr << e.what() << std::endl;
      return LIBRARY_FUNCTION_ERROR;
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSMSSMCalculateSpectrum(
   WolframLibraryData /* libData */, MLINK link)
{
   using namespace MSSM_librarylink;

   if (!check_number_of_args(link, 1, "FSMSSMCalculateSpectrum"))
      return LIBRARY_TYPE_ERROR;

   const Handle_id hid = get_handle_from(link);

   try {
      Model_data data = find_data(hid);

      {
         Redirect_output crd(link);
         calculate_spectrum(data, link);
      }

      check_spectrum(data, link);
      put_spectrum(data.model, link);

      handles[hid] = std::move(data);
   } catch (const flexiblesusy::Error&) {
      put_error_output(link);
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSMSSMCalculateObservables(
   WolframLibraryData /* libData */, MLINK link)
{
   using namespace MSSM_librarylink;

   if (!check_number_of_args(link, 1, "FSMSSMCalculateObservables"))
      return LIBRARY_TYPE_ERROR;

   const Handle_id hid = get_handle_from(link);

   try {
      Model_data data = find_data(hid);

      if (data.model.get_scale() == 0.) {
         put_message(link,
            "FSMSSMCalculateObservables", "warning",
            "Renormalization scale is 0.  Did you run "
            "FSMSSMCalculateSpectrum[]?");
      }

      MSSM_observables observables;

      {
         Redirect_output crd(link);
         observables =
            calculate_observables(data.model, data.qedqcd, data.physical_input);
      }

      put_observables(observables, link);
   } catch (const flexiblesusy::Error& e) {
      put_message(link, "FSMSSMCalculateObservables", "error", e.what());
      put_error_output(link);
   }

   return LIBRARY_NO_ERROR;
}

} // extern "C"
