# work_standard_makefile

# change this as necessary.  see utils/makefile_header for definitions.

include $(MESA_DIR)/utils/makefile_header

LOAD_MESA = $(LOAD_EXTRAS1) $(LOAD_MESA_STAR) $(LOAD_EXTRAS)

#################################################################

ifndef STAR
STAR = star
endif

# STAR_WORK_OBJS defined in makefile_header
OBJS = $(STAR_WORK_OBJS) run.o agn_stars.o mass_changes.o rotation.o atmosphere.o

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src
STAR_JOB_DIR = $(MESA_DIR)/star/job

$(STAR) : $(OBJS)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(STAR) $(OBJS) $(LOAD_MESA)

#################################################################

# WORK_COMPILE defined in makefile_header



ifeq ($(SKIP_NAN_TRAPS),YES)
run.o: $(WORK_SRC_DIR)/run.f90
	$(WORK_COMPILE) $<
else
run.o: $(WORK_SRC_DIR)/run.f90
	$(WORK_COMPILE) $(FCtrapNANs) $<
endif

agn_stars.o: $(WORK_SRC_DIR)/agn_stars.f90
	$(WORK_COMPILE) $(WORK_SRC_DIR)/agn_stars.f90

atmosphere.o: $(WORK_SRC_DIR)/atmosphere.f90 agn_stars.o mass_changes.o
	$(WORK_COMPILE) $(WORK_SRC_DIR)/atmosphere.f90

mass_changes.o: $(WORK_SRC_DIR)/mass_changes.f90 agn_stars.o rotation.o
	$(WORK_COMPILE) $(WORK_SRC_DIR)/mass_changes.f90

rotation.o: $(WORK_SRC_DIR)/rotation.f90 agn_stars.o
	$(WORK_COMPILE) $(WORK_SRC_DIR)/rotation.f90

run_star_extras.o: atmosphere.o rotation.o mass_changes.o agn_stars.o $(WORK_SRC_DIR)/run_star_extras.f90
	$(WORK_COMPILE) $(WORK_SRC_DIR)/run_star_extras.f90

%.o: $(STAR_JOB_DIR)/%.f90
	$(WORK_COMPILE) $<

clean:
	-@rm -f *.o *.mod $(WORK_DIR)/$(STAR)

remk:
	-@rm -f run.o $(WORK_DIR)/$(STAR)









