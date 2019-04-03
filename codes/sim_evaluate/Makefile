CC=g++

all: OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS build
profile: PROFILE_FLAGS build
build: clean-exe evaluate clean

LDFLAGS=#-static
LIBS=-lz -lm
CFLAGS=#-fno-pic 

objects=evaluate.o

# compile: $(objects)

evaluate: $(objects)
	$(CC) $(objects) -o $@ ${LDFLAGS} ${LIBS}

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	@rm -f $(objects)
	
clean-exe:
	@rm -f evaluate

DEBUG_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -ggdb)
	$(eval LIBS = $(LIBS) -ggdb)

OPTIMIZE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -O2)

PROFILE_FLAGS:
		$(eval CFLAGS = $(CFLAGS) -pg -g)
	$(eval LIBS = $(LIBS) -pg -g)
