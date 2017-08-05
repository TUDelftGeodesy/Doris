#########################################
### Makefile to compile envisat tools ###
### based on BEST toolbox of ESA      ###
### see source code for info          ###
### i did a quick hack, but am not    ###
### familar with c (but did not want  ###
### to mix it with c++ for such a     ###
### simple program).                  ###
### Bert Kampes, 15-JUN-2003          ###
###                                   ###
### Support for ERS in Envisat Format ###
### Mahmut Arikan, 19-Oct-2010        ###
###                                   ###
#########################################


### Check if you agree with the options below (install_dir OK?)
### then simply type: "make install"
SHELL = 	/bin/sh
CC = 		gcc
#CFLAGS =	-O3
#CFLAGS =        -m32         # for 64-bit systems, it requires compatibility lib32, no need for epr_api v2.2
LFLAGS =	-lm
INSTALL_DIR =	/home/gertmulder/bin/doris/doris_v5_wu_branch



#############################################################
### THERE SHOULD BE NOTHING YOU HAVE TO CHANGE BELOW HERE ###
#############################################################
EPR_API_SRCDIR = epr_api-2.2/src
LIB_SRC = $(EPR_API_SRCDIR)/epr_api.c \
		$(EPR_API_SRCDIR)/epr_band.c \
		$(EPR_API_SRCDIR)/epr_bitmask.c \
		$(EPR_API_SRCDIR)/epr_core.c \
		$(EPR_API_SRCDIR)/epr_dataset.c \
		$(EPR_API_SRCDIR)/epr_dddb.c \
		$(EPR_API_SRCDIR)/epr_dsd.c \
		$(EPR_API_SRCDIR)/epr_dump.c \
		$(EPR_API_SRCDIR)/epr_field.c \
		$(EPR_API_SRCDIR)/epr_msph.c \
		$(EPR_API_SRCDIR)/epr_param.c \
		$(EPR_API_SRCDIR)/epr_product.c \
		$(EPR_API_SRCDIR)/epr_ptrarray.c \
		$(EPR_API_SRCDIR)/epr_record.c \
		$(EPR_API_SRCDIR)/epr_string.c \
		$(EPR_API_SRCDIR)/epr_swap.c \
		$(EPR_API_SRCDIR)/epr_typconv.c
LIB_OBJ = $(LIB_SRC:.c=.o)

### targets ###
default: 	envisat_dump_header envisat_dump_data envisat_dump_HH envisat_dump_VV  
### compile envisat_dump_header program (dumps envisat header info to stdout) ###
envisat_dump_header: libepr_api.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) \
	-I$(EPR_API_SRCDIR) -L. \
	-o envisat_dump_header \
	envisat_dump_header.c \
	-lepr_api $(LFLAGS)

### compile envisat_dump_data program (dumps envisat slc short int data to file) ###
envisat_dump_data: libepr_api.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) \
	-I$(EPR_API_SRCDIR) -L. \
	-o envisat_dump_data \
	envisat_dump_data.c \
	-lepr_api $(LFLAGS)

### compile envisat_dump_data program (dumps envisat slc short int data to file) ###
envisat_dump_HH: libepr_api.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) \
	-I$(EPR_API_SRCDIR) -L. \
	-o envisat_dump_HH \
	envisat_dump_HH.c \
	-lepr_api $(LFLAGS)

### compile envisat_dump_data program (dumps envisat slc short int data to file) ###
envisat_dump_VV: libepr_api.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) \
	-I$(EPR_API_SRCDIR) -L. \
	-o envisat_dump_VV \
	envisat_dump_VV.c \
	-lepr_api $(LFLAGS)



### create object code BEST library functions ###
libepr_api.a: $(LIB_OBJ)
	ar rcs libepr_api.a $(LIB_OBJ)


### installer ###
PROGS =		envisat_dump_header envisat_dump_data envisat_dump_HH envisat_dump_VV
SCRIPTS =	envisat_dump_header2doris.csh
install:	$(PROGS) $(SCRIPTS)
		@echo "installing envisat tools..."
		mkdir -p $(INSTALL_DIR)
		cp -f $(PROGS) $(SCRIPTS) $(INSTALL_DIR)/.
		#if test -d /cygdrive/c/cygwin; then \
                #  echo "* cygwin detected"; \
		#  cp -f envisat_dump_header $(INSTALL_DIR)/.; \
		#  cp -f envisat_dump_data $(INSTALL_DIR)/.; \
		#  cp -f envisat_dump_header.exe $(INSTALL_DIR)/.; \
		#  cp -f envisat_dump_data.exe $(INSTALL_DIR)/.; \
		#fi;
		@echo " "
		@echo "------------------------------------"
		@echo "FINISHED INSTALLING ENVISAT TOOLS..."
		@echo " executables in: $(INSTALL_DIR)"
		@echo "------------------------------------"
### cleaners ###
uninstall:	cleaner
		@echo "un-installing envisat tools in $(INSTALL_DIR)..."
		(cd $(INSTALL_DIR); rm -f $(PROGS) $(SCRIPTS))
		@echo "now you can rm this directory yourself"
		@echo "for a full uninstall"
		@echo " "
clean:		
		@echo "cleaning envisat tools..."
		rm -f *.o *.a $(LIB_OBJ)
cleaner:	clean
		rm -f $(PROGS)

### EOF.
