PREFIX     = loader
THESYS     = `/bin/uname -s`
EXE        = ../bin/$(PREFIX).$(THESYS)

CXXOBJS  = loader.o

default: $(CXXOBJS)
	$(CXX) $(CXXFLAGS)  $(DYNLINKFLAGS) -o $(EXE) $(CXXOBJS) $(MPLIBS) $(DLLIBS)

clean:
	rm -f $(CXXOBJS) $(EXE)
