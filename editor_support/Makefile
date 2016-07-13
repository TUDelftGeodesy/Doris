#########################################################################
# Makefile for the Doris syntax hilighting using vim.			#
# 22-Mar-2010                                                           #
#                                                                       #
# Delft University of Technology                                        #
# Delft Institute of Earth Observation and Space Systems                #
# http://enterprise.lr.tudelft.nl/doris/				#
# See also the user manual, annex installation.				#
#									#
# Version 1.0 - Mahmut Arikan   					#
#									#
#########################################################################
###################################################################
# Please change if required: 					###
###################################################################
### The shell used by this makefile ###
SHELL   = /bin/sh

### Specify compiler/installation directory ###
SCRIPTSDIR = vim/syntax
INSTALLDIR = ~/.vim/syntax

### Define statements controlling compilation ###
SYNTAX	   = 	doris.vim


#####################################################
### No need to change anything below here... ####
#####################################################
default:	vim

vim:		install

### Install executable in installdir ###
install:	$(SCRIPTSDIR)/$(SYNTAX)
		@echo "* Installing $(SYNTAX) in: $(INSTALLDIR)"
		@mkdir -p $(INSTALLDIR)
		@( cd $(SCRIPTSDIR) 1>/dev/null ;  cp -f $(SYNTAX) $(INSTALLDIR)/. )
		@echo " "
		@echo "*******************************"
		@echo "* ...Installation finished... *"
		@echo "*******************************"
		@echo " "
		@echo "* Check that you have enabled syntax hilighting in .vimrc by adding: "
		@echo "*  syntax on"
		@echo "* "
		@echo "* Type :set syntax=doris to check whether it is working or not. "
		@echo "* "
		@echo " "


#####################################################
### Cleaners ###
uninstall:	
		@rm -f $(INSTALLDIR)/$(SYNTAX)
		@echo "* Removed syntax definition $(SYNTAX) from dir: $(INSTALLDIR)."

### TUDelft 22-03-2010 ###
