CXX     ?= g++

all:     OPTIMIZE_FLAGS build
debug:   DEBUG_FLAGS build
profile: PROFILE_FLAGS build
build:   clean-exe sam2m5 clean

LDFLAGS  = -lz #-static
CXXFLAGS = 

objects  = main.o

sam2m5: $(objects)
	$(CXX) $(objects) -o $@ ${LDFLAGS}

clean:
	@rm -f $(objects)

clean-exe:
	@rm -f sam2m5

DEBUG_FLAGS:
	$(eval CXXFLAGS = $(CXXFLAGS) -ggdb)
	$(eval LIBS = $(LIBS) -ggdb)

OPTIMIZE_FLAGS:
	$(eval CXXFLAGS = $(CXXFLAGS) -O2)

PROFILE_FLAGS:
		$(eval CXXFLAGS = $(CXXFLAGS) -pg -g)
	$(eval LIBS = $(LIBS) -pg -g)
