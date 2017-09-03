DIR          := src
MODNAME      := src
WITH_$(MODNAME) := yes

LIBFLEXI_MK  := \
		$(DIR)/module.mk

LIBFLEXI_SRC := \
		$(DIR)/betafunction.cpp \
		$(DIR)/build_info.cpp \
		$(DIR)/ckm.cpp \
		$(DIR)/command_line_options.cpp \
		$(DIR)/database.cpp \
		$(DIR)/def.cpp \
		$(DIR)/dilog.cpp \
		$(DIR)/dilogc.f \
		$(DIR)/effective_couplings.cpp \
		$(DIR)/error.cpp \
		$(DIR)/gm2calc_interface.cpp \
		$(DIR)/gsl_utils.cpp \
		$(DIR)/linalg.cpp \
		$(DIR)/lowe.cpp \
		$(DIR)/sfermions.cpp \
		$(DIR)/mssm_twoloophiggs.cpp \
		$(DIR)/mssm_twoloophiggs_impl.f \
		$(DIR)/nmssm_twoloophiggs.cpp \
		$(DIR)/nmssm2loop.f \
		$(DIR)/numerics.cpp \
		$(DIR)/numerics2.cpp \
		$(DIR)/spectrum_generator_settings.cpp \
		$(DIR)/physical_input.cpp \
		$(DIR)/pmns.cpp \
		$(DIR)/pv.cpp \
		$(DIR)/rge.cpp \
		$(DIR)/rk.cpp \
		$(DIR)/scan.cpp \
		$(DIR)/slha_io.cpp \
		$(DIR)/sm_twoloophiggs.cpp \
		$(DIR)/split_threeloophiggs.cpp \
		$(DIR)/splitmssm_thresholds.cpp \
		$(DIR)/standard_model.cpp \
		$(DIR)/standard_model_effective_couplings.cpp \
		$(DIR)/standard_model_physical.cpp \
		$(DIR)/standard_model_two_scale_convergence_tester.cpp \
		$(DIR)/standard_model_two_scale_low_scale_constraint.cpp \
		$(DIR)/standard_model_two_scale_model.cpp \
		$(DIR)/threshold_loop_functions.cpp \
		$(DIR)/utils.cpp \
		$(DIR)/weinberg_angle.cpp \
		$(DIR)/wrappers.cpp

LIBFLEXI_HDR := \
		$(DIR)/betafunction.hpp \
		$(DIR)/build_info.hpp \
		$(DIR)/cextensions.hpp \
		$(DIR)/ckm.hpp \
		$(DIR)/command_line_options.hpp \
		$(DIR)/compare.hpp \
		$(DIR)/composite_convergence_tester.hpp \
		$(DIR)/compound_constraint.hpp \
		$(DIR)/constraint.hpp \
		$(DIR)/convergence_tester.hpp \
		$(DIR)/convergence_tester_drbar.hpp \
		$(DIR)/coupling_monitor.hpp \
		$(DIR)/database.hpp \
		$(DIR)/def.h \
		$(DIR)/derivative.hpp \
		$(DIR)/dilog.hpp \
		$(DIR)/effective_couplings.hpp \
		$(DIR)/eigen_utils.hpp \
		$(DIR)/eigen_tensor.hpp \
		$(DIR)/error.hpp \
		$(DIR)/ew_input.hpp \
		$(DIR)/ewsb_solver.hpp \
		$(DIR)/fixed_point_iterator.hpp \
		$(DIR)/functors.hpp \
		$(DIR)/gm2calc_interface.hpp \
		$(DIR)/gsl_utils.hpp \
		$(DIR)/gut_scale_calculator.hpp \
		$(DIR)/two_loop_corrections.hpp \
		$(DIR)/initial_guesser.hpp \
		$(DIR)/linalg2.hpp \
		$(DIR)/linalg.h \
		$(DIR)/logger.hpp \
		$(DIR)/lowe.h \
		$(DIR)/matching.hpp \
		$(DIR)/mathlink_utils.hpp \
		$(DIR)/minimizer.hpp \
		$(DIR)/mssm_twoloophiggs.h \
		$(DIR)/mssm_twoloophiggs.hpp \
		$(DIR)/mycomplex.h \
		$(DIR)/nmssm_twoloophiggs.hpp \
		$(DIR)/nmssm2loop.h \
		$(DIR)/numerics.h \
		$(DIR)/numerics2.hpp \
		$(DIR)/physical_input.hpp \
		$(DIR)/parallel.hpp \
		$(DIR)/pmns.hpp \
		$(DIR)/problems.hpp \
		$(DIR)/pv.hpp \
		$(DIR)/rge.h \
		$(DIR)/rg_flow.hpp \
		$(DIR)/rk.hpp \
		$(DIR)/root_finder.hpp \
		$(DIR)/scan.hpp \
		$(DIR)/sfermions.hpp \
		$(DIR)/slha_io.hpp \
		$(DIR)/sm_twoloophiggs.hpp \
		$(DIR)/split_threeloophiggs.hpp \
		$(DIR)/splitmssm_thresholds.hpp \
		$(DIR)/spectrum_generator_settings.hpp \
		$(DIR)/standard_model.hpp \
		$(DIR)/standard_model_convergence_tester.hpp \
		$(DIR)/standard_model_effective_couplings.hpp \
		$(DIR)/standard_model_low_scale_constraint.hpp \
		$(DIR)/standard_model_physical.hpp \
		$(DIR)/standard_model_two_scale_convergence_tester.hpp \
		$(DIR)/standard_model_two_scale_low_scale_constraint.hpp \
		$(DIR)/standard_model_two_scale_model.hpp \
		$(DIR)/sum.hpp \
		$(DIR)/threshold_loop_functions.hpp \
		$(DIR)/utils.h \
		$(DIR)/weinberg_angle.hpp \
		$(DIR)/wrappers.hpp \
		$(DIR)/xpr-base.h \
		$(DIR)/xpr-matrix.h \
		$(DIR)/xpr-vector.h

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBFLEXI_SRC += \
		$(DIR)/two_scale_composite_convergence_tester.cpp \
		$(DIR)/two_scale_convergence_tester.cpp \
		$(DIR)/two_scale_running_precision.cpp \
		$(DIR)/two_scale_solver.cpp

LIBFLEXI_HDR += \
		$(DIR)/two_scale_composite_convergence_tester.hpp \
		$(DIR)/two_scale_constraint.hpp \
		$(DIR)/two_scale_convergence_tester.hpp \
		$(DIR)/two_scale_convergence_tester_drbar.hpp \
		$(DIR)/two_scale_initial_guesser.hpp \
		$(DIR)/two_scale_matching.hpp \
		$(DIR)/two_scale_model.hpp \
		$(DIR)/two_scale_running_precision.hpp \
		$(DIR)/two_scale_solver.hpp
endif

LIBFLEXI_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBFLEXI_SRC))) \
		$(patsubst %.c, %.o, $(filter %.c, $(LIBFLEXI_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBFLEXI_SRC)))

LIBFLEXI_DEP := \
		$(LIBFLEXI_OBJ:.o=.d)

LIBFLEXI     := $(DIR)/libflexisusy$(MODULE_LIBEXT)

LIBFLEXI_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIBFLEXI)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(LIBFLEXI_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBFLEXI_SRC) $(LIBFLEXI_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBFLEXI_HDR) $(LIBFLEXI_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBFLEXI_MK) $(LIBFLEXI_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBFLEXI_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBFLEXI)

clean-$(MODNAME)-obj:
		-rm -f $(LIBFLEXI_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBFLEXI_DEP) $(LIBFLEXI_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(SQLITEFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBFLEXI_DEP) $(LIBFLEXI_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_SHARED_LIBS),yes)
$(LIBFLEXI): $(LIBFLEXI_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^ $(BOOSTTHREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS)
else
$(LIBFLEXI): $(LIBFLEXI_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^
endif

ALLDEP += $(LIBFLEXI_DEP)
ALLLIB += $(LIBFLEXI)
