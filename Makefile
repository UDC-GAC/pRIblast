BUILD=release

OLEVEL=-O3

cflags.debug=$(OLEVEL) -g3 -march=native -D_DBG_
cflags.release=$(OLEVEL) -s -march=native

cxxflags.debug=$(OLEVEL) -g3 -march=native -D_DBG_
cxxflags.release=$(OLEVEL) -s -march=native

CC=mpicc
CFLAGS=$(cflags.$(BUILD)) -fopenmp

CXX=mpic++
CXXFLAGS=$(cxxflags.$(BUILD)) -fopenmp

LDXX=mpic++
LDXXFLAGS=-fopenmp

SRCSDIR=src
OBJSDIR=obj
DESTDIR=target

SRCS=$(wildcard $(SRCSDIR)/*.c)
SRCS+=$(wildcard $(SRCSDIR)/*.cpp)
TMP=$(patsubst $(SRCSDIR)/%.c,$(OBJSDIR)/$(BUILD)/%.o,$(SRCS))
OBJS=$(patsubst $(SRCSDIR)/%.cpp,$(OBJSDIR)/$(BUILD)/%.o,$(TMP))

pRIblast.$(BUILD): $(OBJS)
	@mkdir -p $(DESTDIR)
	$(LDXX) $(LDXXFLAGS) $^ -o $(DESTDIR)/$@ $(LDXXLIBS)

$(OBJSDIR)/$(BUILD)/%.o: $(SRCSDIR)/%.c
	@mkdir -p $(OBJSDIR)/$(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJSDIR)/$(BUILD)/%.o: $(SRCSDIR)/%.cpp
	@mkdir -p $(OBJSDIR)/$(BUILD)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJSDIR)/$(BUILD) $(DESTDIR)/pRIblast.$(BUILD)

cleanall:
	rm -rf $(OBJSDIR) $(DESTDIR)
