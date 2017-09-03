DIR          := models/EW_triplet
MODNAME      := EW_triplet
SARAH_MODEL  := EW_triplet
WITH_$(MODNAME) := yes

EW_triplet_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

EW_triplet_MK     := \
		$(DIR)/module.mk

EW_triplet_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

EW_triplet_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

EW_triplet_TWO_SCALE_MK := \
		$(EW_triplet_TWO_SCALE_SUSY_MK) \
		$(EW_triplet_TWO_SCALE_SOFT_MK)

EW_triplet_SLHA_INPUT := \
		$(DIR)/LesHouches.in.EW_triplet_generated \
		$(DIR)/LesHouches.in.EW_triplet

EW_triplet_GNUPLOT := \
		$(DIR)/EW_triplet_plot_rgflow.gnuplot \
		$(DIR)/EW_triplet_plot_spectrum.gnuplot

EW_triplet_TARBALL := \
		$(MODNAME).tar.gz

LIBEW_triplet_SRC :=
EXEEW_triplet_SRC :=
LLEW_triplet_LIB  :=
LLEW_triplet_OBJ  :=
LLEW_triplet_SRC  :=
LLEW_triplet_MMA  :=

LIBEW_triplet_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBEW_triplet_SRC += \
		$(DIR)/EW_triplet_effective_couplings.cpp \
		$(DIR)/EW_triplet_mass_eigenstates.cpp \
		$(DIR)/EW_triplet_info.cpp \
		$(DIR)/EW_triplet_input_parameters.cpp \
		$(DIR)/EW_triplet_observables.cpp \
		$(DIR)/EW_triplet_slha_io.cpp \
		$(DIR)/EW_triplet_physical.cpp \
		$(DIR)/EW_triplet_utilities.cpp \
		$(DIR)/EW_triplet_standard_model_matching.cpp \
		$(DIR)/EW_triplet_standard_model_two_scale_matching.cpp \
		$(DIR)/EW_triplet_two_scale_convergence_tester.cpp \
		$(DIR)/EW_triplet_two_scale_high_scale_constraint.cpp \
		$(DIR)/EW_triplet_two_scale_initial_guesser.cpp \
		$(DIR)/EW_triplet_two_scale_low_scale_constraint.cpp \
		$(DIR)/EW_triplet_two_scale_model.cpp \
		$(DIR)/EW_triplet_two_scale_model_slha.cpp \
		$(DIR)/EW_triplet_two_scale_susy_parameters.cpp \
		$(DIR)/EW_triplet_two_scale_soft_parameters.cpp \
		$(DIR)/EW_triplet_two_scale_susy_scale_constraint.cpp
EXEEW_triplet_SRC += \
		$(DIR)/run_EW_triplet.cpp \
		$(DIR)/run_cmd_line_EW_triplet.cpp \
		$(DIR)/scan_EW_triplet.cpp
LIBEW_triplet_HDR += \
		$(DIR)/EW_triplet_convergence_tester.hpp \
		$(DIR)/EW_triplet_effective_couplings.hpp \
		$(DIR)/EW_triplet_high_scale_constraint.hpp \
		$(DIR)/EW_triplet_mass_eigenstates.hpp \
		$(DIR)/EW_triplet_info.hpp \
		$(DIR)/EW_triplet_initial_guesser.hpp \
		$(DIR)/EW_triplet_input_parameters.hpp \
		$(DIR)/EW_triplet_low_scale_constraint.hpp \
		$(DIR)/EW_triplet_model.hpp \
		$(DIR)/EW_triplet_model_slha.hpp \
		$(DIR)/EW_triplet_observables.hpp \
		$(DIR)/EW_triplet_physical.hpp \
		$(DIR)/EW_triplet_slha_io.hpp \
		$(DIR)/EW_triplet_spectrum_generator_interface.hpp \
		$(DIR)/EW_triplet_spectrum_generator.hpp \
		$(DIR)/EW_triplet_standard_model_matching.hpp \
		$(DIR)/EW_triplet_standard_model_two_scale_matching.hpp \
		$(DIR)/EW_triplet_susy_scale_constraint.hpp \
		$(DIR)/EW_triplet_utilities.hpp \
		$(DIR)/EW_triplet_two_scale_convergence_tester.hpp \
		$(DIR)/EW_triplet_two_scale_high_scale_constraint.hpp \
		$(DIR)/EW_triplet_two_scale_initial_guesser.hpp \
		$(DIR)/EW_triplet_two_scale_low_scale_constraint.hpp \
		$(DIR)/EW_triplet_two_scale_model.hpp \
		$(DIR)/EW_triplet_two_scale_model_slha.hpp \
		$(DIR)/EW_triplet_two_scale_soft_parameters.hpp \
		$(DIR)/EW_triplet_two_scale_susy_parameters.hpp \
		$(DIR)/EW_triplet_two_scale_susy_scale_constraint.hpp
LLEW_triplet_SRC  += \
		$(DIR)/EW_triplet_librarylink.cpp

LLEW_triplet_MMA  += \
		$(DIR)/EW_triplet_librarylink.m \
		$(DIR)/run_EW_triplet.m

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(EW_triplet_TWO_SCALE_SUSY_MK)
-include $(EW_triplet_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(EW_triplet_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(EW_triplet_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBEW_triplet_SRC := $(sort $(LIBEW_triplet_SRC))
EXEEW_triplet_SRC := $(sort $(EXEEW_triplet_SRC))

LIBEW_triplet_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBEW_triplet_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBEW_triplet_SRC)))

EXEEW_triplet_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEEW_triplet_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEEW_triplet_SRC)))

EXEEW_triplet_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEEW_triplet_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEEW_triplet_SRC)))

LIBEW_triplet_DEP := \
		$(LIBEW_triplet_OBJ:.o=.d)

EXEEW_triplet_DEP := \
		$(EXEEW_triplet_OBJ:.o=.d)

LLEW_triplet_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLEW_triplet_SRC)))

LLEW_triplet_OBJ  := $(LLEW_triplet_SRC:.cpp=.o)
LLEW_triplet_LIB  := $(LLEW_triplet_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBEW_triplet     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_EW_triplet := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_EW_triplet := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBEW_triplet) $(EXEEW_triplet_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(EW_triplet_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBEW_triplet_SRC) $(EW_triplet_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBEW_triplet_HDR) $(EW_triplet_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEEW_triplet_SRC) $(EW_triplet_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLEW_triplet_SRC) $(EW_triplet_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLEW_triplet_MMA) $(EW_triplet_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(EW_triplet_MK) $(EW_triplet_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(EW_triplet_TWO_SCALE_MK) $(EW_triplet_INSTALL_DIR)
ifneq ($(EW_triplet_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(EW_triplet_SLHA_INPUT) $(EW_triplet_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(EW_triplet_GNUPLOT) $(EW_triplet_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBEW_triplet_DEP)
		-rm -f $(EXEEW_triplet_DEP)
		-rm -f $(LLEW_triplet_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBEW_triplet)
		-rm -f $(LLEW_triplet_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBEW_triplet_OBJ)
		-rm -f $(EXEEW_triplet_OBJ)
		-rm -f $(LLEW_triplet_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEEW_triplet_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(EW_triplet_TARBALL) \
		$(LIBEW_triplet_SRC) $(LIBEW_triplet_HDR) \
		$(EXEEW_triplet_SRC) \
		$(LLEW_triplet_SRC) $(LLEW_triplet_MMA) \
		$(EW_triplet_MK) $(EW_triplet_TWO_SCALE_MK) \
		$(EW_triplet_SLHA_INPUT) $(EW_triplet_GNUPLOT)

$(LIBEW_triplet_SRC) $(LIBEW_triplet_HDR) $(EXEEW_triplet_SRC) $(LLEW_triplet_SRC) $(LLEW_triplet_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_EW_triplet)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_EW_triplet): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_EW_triplet)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_EW_triplet)"
		@echo "Note: to regenerate EW_triplet source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_EW_triplet)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_EW_triplet):
		@true
endif

$(LIBEW_triplet_DEP) $(EXEEW_triplet_DEP) $(LLEW_triplet_DEP) $(LIBEW_triplet_OBJ) $(EXEEW_triplet_OBJ) $(LLEW_triplet_OBJ) $(LLEW_triplet_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBEW_triplet_DEP) $(EXEEW_triplet_DEP) $(LLEW_triplet_DEP) $(LIBEW_triplet_OBJ) $(EXEEW_triplet_OBJ) $(LLEW_triplet_OBJ) $(LLEW_triplet_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLEW_triplet_OBJ) $(LLEW_triplet_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBEW_triplet): $(LIBEW_triplet_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBEW_triplet) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLEW_triplet_LIB): $(LLEW_triplet_OBJ) $(LIBEW_triplet) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBEW_triplet_DEP) $(EXEEW_triplet_DEP)
ALLSRC += $(LIBEW_triplet_SRC) $(EXEEW_triplet_SRC)
ALLLIB += $(LIBEW_triplet)
ALLEXE += $(EXEEW_triplet_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLEW_triplet_DEP)
ALLSRC += $(LLEW_triplet_SRC)
ALLLL  += $(LLEW_triplet_LIB)
endif
