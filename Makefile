PD = ../packages
PN = FLa4a
SD = $(PD)/$(PN)/inst/admb
BD = $(PD)/$(PN)/inst/bin/linux
DD = test
LD = /home/millaco/R/x86_64-pc-linux-gnu-library/3.0/FLa4a/R
sourcefiles := $(wildcard $(PD)/$(PN)/R/*.R)

.PHONY = roxygen install compile run data clean showlog

all: install

clean:
	rm -rf test

compile: $(BD)/a4a
$(BD)/a4a: $(SD)/a4a.tpl $(SD)/nLogNormal.h 
	rm -rf _tmp
	mkdir _tmp
	cp $(SD)/a4a.tpl $(SD)/nLogNormal.h _tmp
	cd _tmp; admb -s a4a 
	cp _tmp/a4a $(BD)
	rm -rf _tmp

roxygen: $(PD)/$(PN)/man/FLa4a-package.Rd
$(PD)/$(PN)/man/FLa4a-package.Rd: $(sourcefiles)
	rm -f $(PD)/$(PN)/man/*
	cd $(PD); echo 'library(roxygen2); roxygenize("FLa4a", roclets = c("namespace","rd"))' | R --vanilla --slave

install: $(LD)/FLa4a 
$(LD)/FLa4a: $(PD)/$(PN)/man/FLa4a-package.Rd $(BD)/a4a
	cd $(PD); R CMD INSTALL FLa4a

data: $(DD)/a4a.dat
$(DD)/a4a.dat: $(LD)/FLa4a makeTest.R
	rm -rf test
	echo 'source("makeTest.R")' | R --vanilla --slave

run: $(DD)/program.log
$(DD)/program.log: $(DD)/a4a.dat 
	cp $(BD)/a4a $(DD)
	cd $(DD); ./a4a

showlog: $(DD)/program.log
	less $(DD)/program.log
