MODULES = handAna.so 
EXECUTABLES = 

all: $(MODULES) $(EXECUTABLES)
mods : $(MODULES)
exes : $(EXECUTABLES)

LIBS = -L$(BELLE_RUN_DIR)/lib/so -ltuple -lip -lkid -leid -lparticle \
       -lkfitter -lbenergy -lbelleCLHEP -lcrypt -lm \
       -lgcc $(shell root-config --libs)
SOFLAGS = -shared -Wl,-export-dynamic $(shell root-config --ldflags) $(shell ../fastjet/fastjet-install/bin/fastjet-config --libs --plugins)

CPPFLAGS = -DHAVE_LIBCURSES=1 -DHAVE_LIBREADLINE=1 -DHAVE_POSTGRES=1 \
	   -DHAVE_LIBCURSES=1 -DHAVE_LIBTERMCAP=1 -DHAVE_LIBHISTORY=1 \
           -DHAVE_LIBREADLINE=1 -DHAVE_HISTORY=1 -DHAVE_LIBBSD=1 \
           -DHAVE_LIBM=1 -DHAVE_LIBDL=1 -DHAVE_LIBNSL=1 -DHAVE_LIBCRYPT=1 \
	   -DHAVE_LIBNSL=1 -DHAVE_LIBDL=1 -DFORTRAN_PPU=1 -DHAVE_LIBCRYPT=1 \
           -DCERNLIB_TYPE -DHEP_SHORT_NAMES -DBELLE_SHORT_NAMES \
           -DDSTXX_NOINLINE -DBELLE_TARGET_H=\"belle-default.h\" \
	   -fPIC -ggdb $(shell root-config --cflags) $(shell ../fastjet/fastjet-install/bin/fastjet-config --cxxflags --plugins)

INCLUDE = -I$(BELLE_TOP_DIR)/include -I. 

%.o : %.cc
	g++ -c $^ $(CPPFLAGS) $(INCLUDE)

%.so : *.o
	g++ $(SOFLAGS) -o $@ $< $(LIBS) 



install: $(MODULES) $(EXECUTABLES)
	/usr/bin/install $(MODULES) $(MY_TOP_DIR)/basf_modules
#	/usr/bin/install $(EXECUTABLES) $(MY_TOP_DIR)/bin	

.PHONY: clean

clean: 
	rm $(MODULES) $(patsubst %.so,%.o,$(MODULES)) *~


