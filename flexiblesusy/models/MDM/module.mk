DIR          := models/MDM
MODNAME      := MDM
SARAH_MODEL  := MDM
WITH_$(MODNAME) := yes

MDM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

MDM_MK     := \
		$(DIR)/module.mk

MDM_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

MDM_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

MDM_TWO_SCALE_MK := \
		$(MDM_TWO_SCALE_SUSY_MK) \
		$(MDM_TWO_SCALE_SOFT_MK)

MDM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MDM_generated \
		$(DIR)/LesHouches.in.MDM

MDM_GNUPLOT := \
		$(DIR)/MDM_plot_rgflow.gnuplot \
		$(DIR)/MDM_plot_spectrum.gnuplot

MDM_TARBALL := \
		$(MODNAME).tar.gz

LIBMDM_SRC :=
EXEMDM_SRC :=
LLMDM_LIB  :=
LLMDM_OBJ  :=
LLMDM_SRC  :=
LLMDM_MMA  :=

LIBMDM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBMDM_SRC += \
		$(DIR)/MDM_effective_couplings.cpp \
		$(DIR)/MDM_mass_eigenstates.cpp \
		$(DIR)/MDM_info.cpp \
		$(DIR)/MDM_input_parameters.cpp \
		$(DIR)/MDM_observables.cpp \
		$(DIR)/MDM_slha_io.cpp \
		$(DIR)/MDM_physical.cpp \
		$(DIR)/MDM_utilities.cpp \
		$(DIR)/MDM_standard_model_matching.cpp \
		$(DIR)/MDM_standard_model_two_scale_matching.cpp \
		$(DIR)/MDM_two_scale_convergence_tester.cpp \
		$(DIR)/MDM_two_scale_high_scale_constraint.cpp \
		$(DIR)/MDM_two_scale_initial_guesser.cpp \
		$(DIR)/MDM_two_scale_low_scale_constraint.cpp \
		$(DIR)/MDM_two_scale_model.cpp \
		$(DIR)/MDM_two_scale_model_slha.cpp \
		$(DIR)/MDM_two_scale_susy_parameters.cpp \
		$(DIR)/MDM_two_scale_soft_parameters.cpp \
		$(DIR)/MDM_two_scale_susy_scale_constraint.cpp
EXEMDM_SRC += \
		$(DIR)/run_MDM.cpp \
		$(DIR)/run_cmd_line_MDM.cpp \
		$(DIR)/scan_MDM.cpp
LIBMDM_HDR += \
		$(DIR)/MDM_convergence_tester.hpp \
		$(DIR)/MDM_effective_couplings.hpp \
		$(DIR)/MDM_high_scale_constraint.hpp \
		$(DIR)/MDM_mass_eigenstates.hpp \
		$(DIR)/MDM_info.hpp \
		$(DIR)/MDM_initial_guesser.hpp \
		$(DIR)/MDM_input_parameters.hpp \
		$(DIR)/MDM_low_scale_constraint.hpp \
		$(DIR)/MDM_model.hpp \
		$(DIR)/MDM_model_slha.hpp \
		$(DIR)/MDM_observables.hpp \
		$(DIR)/MDM_physical.hpp \
		$(DIR)/MDM_slha_io.hpp \
		$(DIR)/MDM_spectrum_generator_interface.hpp \
		$(DIR)/MDM_spectrum_generator.hpp \
		$(DIR)/MDM_standard_model_matching.hpp \
		$(DIR)/MDM_standard_model_two_scale_matching.hpp \
		$(DIR)/MDM_susy_scale_constraint.hpp \
		$(DIR)/MDM_utilities.hpp \
		$(DIR)/MDM_two_scale_convergence_tester.hpp \
		$(DIR)/MDM_two_scale_high_scale_constraint.hpp \
		$(DIR)/MDM_two_scale_initial_guesser.hpp \
		$(DIR)/MDM_two_scale_low_scale_constraint.hpp \
		$(DIR)/MDM_two_scale_model.hpp \
		$(DIR)/MDM_two_scale_model_slha.hpp \
		$(DIR)/MDM_two_scale_soft_parameters.hpp \
		$(DIR)/MDM_two_scale_susy_parameters.hpp \
		$(DIR)/MDM_two_scale_susy_scale_constraint.hpp
LLMDM_SRC  += \
		$(DIR)/MDM_librarylink.cpp

LLMDM_MMA  += \
		$(DIR)/MDM_librarylink.m \
		$(DIR)/run_MDM.m

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(MDM_TWO_SCALE_SUSY_MK)
-include $(MDM_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MDM_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MDM_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBMDM_SRC := $(sort $(LIBMDM_SRC))
EXEMDM_SRC := $(sort $(EXEMDM_SRC))

LIBMDM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMDM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMDM_SRC)))

EXEMDM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMDM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMDM_SRC)))

EXEMDM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMDM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMDM_SRC)))

LIBMDM_DEP := \
		$(LIBMDM_OBJ:.o=.d)

EXEMDM_DEP := \
		$(EXEMDM_OBJ:.o=.d)

LLMDM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMDM_SRC)))

LLMDM_OBJ  := $(LLMDM_SRC:.cpp=.o)
LLMDM_LIB  := $(LLMDM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMDM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MDM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MDM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMDM) $(EXEMDM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(MDM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMDM_SRC) $(MDM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMDM_HDR) $(MDM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEMDM_SRC) $(MDM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMDM_SRC) $(MDM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMDM_MMA) $(MDM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(MDM_MK) $(MDM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(MDM_TWO_SCALE_MK) $(MDM_INSTALL_DIR)
ifneq ($(MDM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(MDM_SLHA_INPUT) $(MDM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(MDM_GNUPLOT) $(MDM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBMDM_DEP)
		-rm -f $(EXEMDM_DEP)
		-rm -f $(LLMDM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBMDM)
		-rm -f $(LLMDM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBMDM_OBJ)
		-rm -f $(EXEMDM_OBJ)
		-rm -f $(LLMDM_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEMDM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(MDM_TARBALL) \
		$(LIBMDM_SRC) $(LIBMDM_HDR) \
		$(EXEMDM_SRC) \
		$(LLMDM_SRC) $(LLMDM_MMA) \
		$(MDM_MK) $(MDM_TWO_SCALE_MK) \
		$(MDM_SLHA_INPUT) $(MDM_GNUPLOT)

$(LIBMDM_SRC) $(LIBMDM_HDR) $(EXEMDM_SRC) $(LLMDM_SRC) $(LLMDM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MDM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MDM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MDM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_MDM)"
		@echo "Note: to regenerate MDM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MDM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MDM):
		@true
endif

$(LIBMDM_DEP) $(EXEMDM_DEP) $(LLMDM_DEP) $(LIBMDM_OBJ) $(EXEMDM_OBJ) $(LLMDM_OBJ) $(LLMDM_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMDM_DEP) $(EXEMDM_DEP) $(LLMDM_DEP) $(LIBMDM_OBJ) $(EXEMDM_OBJ) $(LLMDM_OBJ) $(LLMDM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMDM_OBJ) $(LLMDM_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBMDM): $(LIBMDM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMDM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLMDM_LIB): $(LLMDM_OBJ) $(LIBMDM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBMDM_DEP) $(EXEMDM_DEP)
ALLSRC += $(LIBMDM_SRC) $(EXEMDM_SRC)
ALLLIB += $(LIBMDM)
ALLEXE += $(EXEMDM_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMDM_DEP)
ALLSRC += $(LLMDM_SRC)
ALLLL  += $(LLMDM_LIB)
endif
