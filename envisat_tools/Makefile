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
CC = 		g++
#CFLAGS =	-O3
#CFLAGS =        -m32         # for 64-bit systems, it requires compatibility lib32, no need for epr_api v2.2
LFLAGS =	-lm
# INSTALL_DIR =	/home/gertmulder/bin/doris/doris_v5_wu_branch
INSTALL_DIR = /usr/local/bin


#############################################################
### THERE SHOULD BE NOTHING YOU HAVE TO CHANGE BELOW HERE ###
#############################################################
LIB_SRC = 	epr_api-2.2/src/epr_api.c \
		epr_api-2.2/src/epr_band.c \
		epr_api-2.2/src/epr_bitmask.c \
		epr_api-2.2/src/epr_core.c \
		epr_api-2.2/src/epr_dataset.c \
		epr_api-2.2/src/epr_dddb.c \
		epr_api-2.2/src/epr_dsd.c \
		epr_api-2.2/src/epr_dump.c \
		epr_api-2.2/src/epr_field.c \
		epr_api-2.2/src/epr_msph.c \
		epr_api-2.2/src/epr_param.c \
		epr_api-2.2/src/epr_product.c \
		epr_api-2.2/src/epr_ptrarray.c \
		epr_api-2.2/src/epr_record.c \
		epr_api-2.2/src/epr_string.c \
		epr_api-2.2/src/epr_swap.c \
		epr_api-2.2/src/epr_typconv.c

### targets ###
default: 	envisat_dump_header envisat_dump_data envisat_dump_HH envisat_dump_VV  
### compile envisat_dump_header program (dumps envisat header info to stdout) ###
envisat_dump_header:	BEST
		$(CC) $(CFLAGS) $(LFLAGS) \
		-o envisat_dump_header \
		epr_api-2.2/src/*.o \
		envisat_dump_header.c

### compile envisat_dump_data program (dumps envisat slc short int data to file) ###
envisat_dump_data:	BEST
		$(CC) $(CFLAGS) $(LFLAGS) \
		-o envisat_dump_data \
		epr_api-2.2/src/*.o \
		envisat_dump_data.c

### compile envisat_dump_data program (dumps envisat slc short int data to file) ###
envisat_dump_HH:  BEST
		$(CC) $(CFLAGS) $(LFLAGS) \
    -o envisat_dump_HH \
    epr_api-2.2/src/*.o \
    envisat_dump_HH.c

### compile envisat_dump_data program (dumps envisat slc short int data to file) ###
envisat_dump_VV:  BEST
		$(CC) $(CFLAGS) $(LFLAGS) \
    -o envisat_dump_VV \
    epr_api-2.2/src/*.o \
    envisat_dump_VV.c



### create object code BEST library functions ###
BEST:           $(LIB_SRC)
		(cd epr_api-2.2/src/; make)


### installer ###
PROGS =		envisat_dump_header envisat_dump_data envisat_dump_HH envisat_dump_VV
SCRIPTS =	envisat_dump_header2doris.csh
install:	$(PROGS) $(SCRIPTS)
		@echo "installing envisat tools..."
		cp -f $(PROGS) $(SCRIPTS) $(INSTALL_DIR)/.
		#if test -d /cygdrive/c/cygwin; then \
                #  echo "* cygwin detected"; \
		#  cp -f envisat_dump_header $(INSTALL_DIR)/.; \
		#  cp -f envisat_dump_data $(INSTALL_DIR)/.; \
		#  cp -f envisat_dump_header.exe $(INSTALL_DIR)/.; \
		#  cp -f envisat_dump_data.exe $(INSTALL_DIR)/.; \
		#fi;
		make cleaner;
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
		(cd epr_api-2.2/src/; make clean)
		rm -f *.o
cleaner:	clean
		rm -f $(PROGS)

### EOF.
