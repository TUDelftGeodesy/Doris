                                 README

                     ENVISAT Product Reader API for C
                              Version 2.2

                              31. July 2010


Contents

  o Introduction
  o The EPR-API for C
  o Zip Archive Contents
  o Documentation
  o Release Notes
  o Bug Reports and Feedback


Introduction

    Thank you for downloading this release of the EPR API package.
    The ENVISAT Product Reader API is a set of C-source code files
    supporting developers who want to use MERIS, AATSR, and ASAR data
    products of the ESA ENVISAT satellite in their software.

    The ENVISAT Product Reader API is developed and distributed under the
    terms of open source software development.

    The main use case for the C API is the ingestion of ENVISAT data
    into
      o new scientific algorithms developed for the MERIS, AATSR or ASAR
        sensors or even all of them,
      o existing scientific software packages written in C or C++,
        or COTS software systems which allow for extension using a
        C or C++ interface.


The EPR-API for C (epr-c-api)

    The epr-c-api is a rich library of functions for reading data from ENVISAT
    MERIS, AATSR Level 1b and Level 2 and also ASAR data products.
    The API provides access to the data either on a geophysical (decoded,
    ready-to-use pixel samples) or on a raw data layer. The raw data access
    makes it possible to read any data field contained in a product file.

    The epr-c-api is written in pure ANSI C and should compile with every
    ANSI conformant C-compiler. Just copy all '*.h' and '*.c' files contained in
	the 'src' folder of the distribution to your own source folder of your project
    and include them in your development environment.

    Examples:
	
       You'll find some starting points for your own applications in the examples
       directory:
	   
       epr_api-XXX/
         |- src/examples/         - contains source code examples for developers
           |- write_bands.c
           |- write_bitmask.c
           |- write_ndvi.c

       Compile the examples with console:
          
		  To compile the examples, use the scripts in:
		  
          epr_api-XXX/
            |- prj/unix

          Windows users can use msys (http://www.mingw.org/wiki/msys) or
          cygwin (http://www.cygwin.com/) to compile the examples.

       Compile the examples with an IDE:
          
		  Also all users can use Code::Blocks (http://www.codeblocks.org/) a cross
          platform  IDE (integrated developmenht engine) to compile an work with the
          examples and the api-source code.

          Load the project settings from:
          
		  epr_api-XXX/
            |- prj/codeblocks

		Build a library of the epr-c-api:
		
			A library of the epr-c-api can be built by using the 'makefile' which is
			provided with the distribution.  Note that the 'makefile' is provided as
			is. The manner how to build a library depends on the operating system as
			well as on the type of library being built.
			
			Settings are provided for building shared (or dynamic) libraries for Linux,
			Solaris and Mac OS X. For other operating systems please consider your
			system guidelines in order to find out the correct settings. 
			

Release Notes

    Although this is rather advanced release of the EPR API Software,
    it is still under development.
    Naturally, this software will improve and grow through
    extensive testing and user feedback.


Documentation

    The documentation is located at:

    epr_api-XXX/
      |- docs/
        |- html/    - API documentation for C

    Aditional documentation listing the supported ENVISAT product tables (required as
	a reference  for dataset-, field-, band- and flag-names) can be found at:

    http://www.brockmann-consult.de/beam/doc/help/general/envisat-products/index.html


Zip Archive Contents

    After a successfull unzipping of the epr_api-XXX.zip archive you should get
    a folder named epr_api-XXX with following structure:

    epr_api-XXX/
     |- docs/             - contains the api documentation as html
     |- prj/              - provides build scripts for the C examples and the
                            Code::Blocks IDE workspace
     |- src/                  - epr-c-api C source code
     |- src/examples/         - contains source code examples for C
     |- README.txt            - this file
     |- VERSION.txt           - the version string
     |- CHANGELOG.txt         - the changes from version to version
     |- LICENSE.txt           - the GNU public license
     |- makefile          - the makefile to buld the epr-c-api static library
                            ready to use in your projects


Bug Reports and Feedback

    Please submit your bug reports via e-mail to
    beam-issues@brockmann-consult.de
    and tell us which bug fixes matter most to you.

    You can also send comments directly to the software engineering team
    norman.fomferra@brockmann-consult.de
    sabine.embacher@brockmann-consult.de

--------------------------------------------------------------------------------
Thank you for using the EPR API Software.
2009-04-28 Geesthacht, Germany
                                                      Brockmann Consult
                                                  www.brockmann-consult.com

--------------------------------------------------------------------------------
The EPR API Software is licensed under the terms of the GNU General Public License
aggreement (see the accompanying LICENSE.txt file).
--------------------------------------------------------------------------------
Brockmann Consult 2010
