all: OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS build
profile: PROFILE_FLAGS build
build: clean-exe extract clean

LDFLAGS=#-static
LIBS=-lz# -lm -pthread -lpthread
CFLAGS=#-fno-pic 

objects=extract.o

# compile: $(objects)

extract: $(objects)
	$(CXX) $(objects) -o $@ ${LDFLAGS} ${LIBS}

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $<

clean:
	@rm -f $(objects)
	
clean-exe: clean
	@rm -f extract

DEBUG_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -ggdb)
	$(eval LIBS = $(LIBS) -ggdb)

OPTIMIZE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -O2)

PROFILE_FLAGS:
		$(eval CFLAGS = $(CFLAGS) -pg -g)
	$(eval LIBS = $(LIBS) -pg -g)
