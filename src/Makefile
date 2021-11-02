ADRM_DIR = adapterremoval
HISAT_DIR = hisat
SHELL=cmd.exe
CXX = g++

.PHONY: all clean

all: 
	mkdir  -p ../inst
	(cd	$(ADRM_DIR) && ($(MAKE) CXX="$(CXX) -std=gnu++11" -f Makefile))
	(cd	$(ADRM_DIR) && ($(MAKE) clean_o -f Makefile))
	(cd $(HISAT_DIR) && ($(MAKE) CXX="$(CXX)" -f Makefile))
	(cd	$(HISAT_DIR) && ($(MAKE) move -f Makefile))
	(cd	$(HISAT_DIR) && ($(MAKE) clean_dSYM -f Makefile))
	$(CXX) -shared -fPIC -o Rhisat.so  version_info.cpp


clean:
	(cd	$(ADRM_DIR) && ($(MAKE) clean -f Makefile))
	(cd	$(HISAT_DIR) && ($(MAKE) clean -f Makefile))
	rm -f *.so
	rm -f *.o




