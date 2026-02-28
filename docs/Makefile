# Makefile for MOM6 documentation (doxygen)

# You can set these variables from the command line.
BUILDDIR         ?= _build
DOXYGEN_BIN      ?= doxygen
DOXYGEN_CONF     ?= Doxyfile_rtd
DOXYGEN_RELEASE  ?= Release_1_9_8

# Set DOXYGENBIN variable if the executable is found. If not, set it to
# a makefile rule that will clone and build doxygen from source (for
# developers without root access to install a system package).
ifeq ($(shell which $(DOXYGEN_BIN) >/dev/null 2>&1; echo $$?), 0)
DOXYGENBIN = $(shell which $(DOXYGEN_BIN))
else
DOXYGENBIN = doxygen/bin/doxygen
endif

.PHONY: help clean html nortd nortd_latex

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html         to make HTML files using doxygen (RTD configuration)"
	@echo "  nortd        to make HTML files using doxygen (standalone configuration)"
	@echo "  nortd_latex  to make PDF via doxygen LaTeX output"
	@echo "  clean        to remove generated files"

clean:
	rm -rf $(BUILDDIR) html xml MOM6.tags

$(BUILDDIR):
	mkdir -p $@

html: $(DOXYGENBIN)
	$(DOXYGENBIN) $(DOXYGEN_CONF)
	@echo
	@echo "Build finished. The HTML pages are in html/."

doxygen:
	git clone -b $(DOXYGEN_RELEASE) --depth 1 https://github.com/doxygen/doxygen

doxygen/bin/doxygen: doxygen
	(cd $(<); cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -G "Unix Makefiles" .)
	(cd $(<); $(MAKE))

nortd: $(DOXYGENBIN) $(BUILDDIR)
ifdef DOXYGEN_CONF
	$(DOXYGENBIN) $(DOXYGEN_CONF)
else
	$(DOXYGENBIN) Doxyfile_nortd
endif
ifeq ($(UPDATEHTMLEQS), Y)
	@echo "Post processing equations."
	@python3 ./postProcessEquations.py -d $(BUILDDIR) -p APIs -b doxygen -s index.html $(UPDATEHTMLEQSVERBOSE)
endif
	@echo "Build finished. The html files are in $(BUILDDIR)/APIs."

nortd_latex: $(DOXYGENBIN) $(BUILDDIR)
ifdef DOXYGEN_CONF
	$(DOXYGENBIN) $(DOXYGEN_CONF)
else
	$(DOXYGENBIN) Doxyfile_nortd_latex
endif
	# Make sure we are in nonstopmode for generating the pdf by patching the Makefile
	# Older doxygen versions do not have LATEX_CMD
	-grep -q "pdflatex refman" $(BUILDDIR)/doxygen_latex/Makefile && sed -i.bak 's/pdflatex refman/pdflatex -interaction=nonstopmode refman/' $(BUILDDIR)/doxygen_latex/Makefile
	make -C $(BUILDDIR)/doxygen_latex LATEX_CMD="pdflatex -interaction=nonstopmode"

Clean: clean

CLEAN: Clean
	rm -rf doxygen
