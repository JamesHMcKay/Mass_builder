DIR          := models/MSSM
MODNAME      := MSSM
SARAH_MODEL  := MSSM
WITH_$(MODNAME) := yes

MSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

MSSM_MK     := \
		$(DIR)/module.mk

MSSM_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

MSSM_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

MSSM_TWO_SCALE_MK := \
		$(MSSM_TWO_SCALE_SUSY_MK) \
		$(MSSM_TWO_SCALE_SOFT_MK)

MSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MSSM_generated \
		$(DIR)/LesHouches.in.MSSM

MSSM_GNUPLOT := \
		$(DIR)/MSSM_plot_rgflow.gnuplot \
		$(DIR)/MSSM_plot_spectrum.gnuplot

MSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBMSSM_SRC :=
EXEMSSM_SRC :=
LLMSSM_LIB  :=
LLMSSM_OBJ  :=
LLMSSM_SRC  :=
LLMSSM_MMA  :=

LIBMSSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBMSSM_SRC += \
		$(DIR)/MSSM_effective_couplings.cpp \
		$(DIR)/MSSM_mass_eigenstates.cpp \
		$(DIR)/MSSM_info.cpp \
		$(DIR)/MSSM_input_parameters.cpp \
		$(DIR)/MSSM_observables.cpp \
		$(DIR)/MSSM_slha_io.cpp \
		$(DIR)/MSSM_physical.cpp \
		$(DIR)/MSSM_utilities.cpp \
		$(DIR)/MSSM_standard_model_matching.cpp \
		$(DIR)/MSSM_standard_model_two_scale_matching.cpp \
		$(DIR)/MSSM_two_scale_convergence_tester.cpp \
		$(DIR)/MSSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/MSSM_two_scale_initial_guesser.cpp \
		$(DIR)/MSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/MSSM_two_scale_model.cpp \
		$(DIR)/MSSM_two_scale_model_slha.cpp \
		$(DIR)/MSSM_two_scale_susy_parameters.cpp \
		$(DIR)/MSSM_two_scale_soft_parameters.cpp \
		$(DIR)/MSSM_two_scale_susy_scale_constraint.cpp
EXEMSSM_SRC += \
		$(DIR)/run_MSSM.cpp \
		$(DIR)/run_cmd_line_MSSM.cpp \
		$(DIR)/scan_MSSM.cpp
LIBMSSM_HDR += \
		$(DIR)/MSSM_convergence_tester.hpp \
		$(DIR)/MSSM_effective_couplings.hpp \
		$(DIR)/MSSM_high_scale_constraint.hpp \
		$(DIR)/MSSM_mass_eigenstates.hpp \
		$(DIR)/MSSM_info.hpp \
		$(DIR)/MSSM_initial_guesser.hpp \
		$(DIR)/MSSM_input_parameters.hpp \
		$(DIR)/MSSM_low_scale_constraint.hpp \
		$(DIR)/MSSM_model.hpp \
		$(DIR)/MSSM_model_slha.hpp \
		$(DIR)/MSSM_observables.hpp \
		$(DIR)/MSSM_physical.hpp \
		$(DIR)/MSSM_slha_io.hpp \
		$(DIR)/MSSM_spectrum_generator_interface.hpp \
		$(DIR)/MSSM_spectrum_generator.hpp \
		$(DIR)/MSSM_standard_model_matching.hpp \
		$(DIR)/MSSM_standard_model_two_scale_matching.hpp \
		$(DIR)/MSSM_susy_scale_constraint.hpp \
		$(DIR)/MSSM_utilities.hpp \
		$(DIR)/MSSM_two_scale_convergence_tester.hpp \
		$(DIR)/MSSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/MSSM_two_scale_initial_guesser.hpp \
		$(DIR)/MSSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/MSSM_two_scale_model.hpp \
		$(DIR)/MSSM_two_scale_model_slha.hpp \
		$(DIR)/MSSM_two_scale_soft_parameters.hpp \
		$(DIR)/MSSM_two_scale_susy_parameters.hpp \
		$(DIR)/MSSM_two_scale_susy_scale_constraint.hpp
LLMSSM_SRC  += \
		$(DIR)/MSSM_librarylink.cpp

LLMSSM_MMA  += \
		$(DIR)/MSSM_librarylink.m \
		$(DIR)/run_MSSM.m

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(MSSM_TWO_SCALE_SUSY_MK)
-include $(MSSM_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MSSM_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSM_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

endif

# remove duplicates in case all algorithms are used
LIBMSSM_SRC := $(sort $(LIBMSSM_SRC))
EXEMSSM_SRC := $(sort $(EXEMSSM_SRC))

LIBMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSM_SRC)))

EXEMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMSSM_SRC)))

EXEMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMSSM_SRC)))

LIBMSSM_DEP := \
		$(LIBMSSM_OBJ:.o=.d)

EXEMSSM_DEP := \
		$(EXEMSSM_OBJ:.o=.d)

LLMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMSSM_SRC)))

LLMSSM_OBJ  := $(LLMSSM_SRC:.cpp=.o)
LLMSSM_LIB  := $(LLMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMSSM) $(EXEMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(MSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSM_SRC) $(MSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSM_HDR) $(MSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEMSSM_SRC) $(MSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMSSM_SRC) $(MSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMSSM_MMA) $(MSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(MSSM_MK) $(MSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(MSSM_TWO_SCALE_MK) $(MSSM_INSTALL_DIR)
ifneq ($(MSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(MSSM_SLHA_INPUT) $(MSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(MSSM_GNUPLOT) $(MSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBMSSM_DEP)
		-rm -f $(EXEMSSM_DEP)
		-rm -f $(LLMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBMSSM)
		-rm -f $(LLMSSM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBMSSM_OBJ)
		-rm -f $(EXEMSSM_OBJ)
		-rm -f $(LLMSSM_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(MSSM_TARBALL) \
		$(LIBMSSM_SRC) $(LIBMSSM_HDR) \
		$(EXEMSSM_SRC) \
		$(LLMSSM_SRC) $(LLMSSM_MMA) \
		$(MSSM_MK) $(MSSM_TWO_SCALE_MK) \
		$(MSSM_SLHA_INPUT) $(MSSM_GNUPLOT)

$(LIBMSSM_SRC) $(LIBMSSM_HDR) $(EXEMSSM_SRC) $(LLMSSM_SRC) $(LLMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_MSSM)"
		@echo "Note: to regenerate MSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MSSM):
		@true
endif

$(LIBMSSM_DEP) $(EXEMSSM_DEP) $(LLMSSM_DEP) $(LIBMSSM_OBJ) $(EXEMSSM_OBJ) $(LLMSSM_OBJ) $(LLMSSM_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMSSM_DEP) $(EXEMSSM_DEP) $(LLMSSM_DEP) $(LIBMSSM_OBJ) $(EXEMSSM_OBJ) $(LLMSSM_OBJ) $(LLMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMSSM_OBJ) $(LLMSSM_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBMSSM): $(LIBMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLMSSM_LIB): $(LLMSSM_OBJ) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBMSSM_DEP) $(EXEMSSM_DEP)
ALLSRC += $(LIBMSSM_SRC) $(EXEMSSM_SRC)
ALLLIB += $(LIBMSSM)
ALLEXE += $(EXEMSSM_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMSSM_DEP)
ALLSRC += $(LLMSSM_SRC)
ALLLL  += $(LLMSSM_LIB)
endif
