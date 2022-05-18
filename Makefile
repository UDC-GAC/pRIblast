BUILD=release

OLEVEL=-O3

cxxflags.debug=-D_DBG_
cxxflags.release=

CXX=mpic++
CXXFLAGS=$(OLEVEL) -march=native -s -fopenmp $(cxxflags.$(BUILD))

LD=mpic++
LDFLAGS=-fopenmp
LDLIBS=

SRCSDIR=src
OBJSDIR=obj
DESTDIR=target

SRCS=$(wildcard $(SRCSDIR)/*.c)
SRCS+=$(wildcard $(SRCSDIR)/*.cpp)
TMP=$(patsubst $(SRCSDIR)/%.c,$(OBJSDIR)/$(BUILD)/%.o,$(SRCS))
OBJS=$(patsubst $(SRCSDIR)/%.cpp,$(OBJSDIR)/$(BUILD)/%.o,$(TMP))

pRIblast.$(BUILD): $(OBJS)
	@mkdir -p $(DESTDIR)
	$(LD) $(LDFLAGS) -o $(DESTDIR)/$@ $^ $(LDLIBS)

$(OBJSDIR)/$(BUILD)/%.o: $(SRCSDIR)/%.c
	@mkdir -p $(OBJSDIR)/$(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(OBJSDIR)/$(BUILD)/%.o: $(SRCSDIR)/%.cpp
	@mkdir -p $(OBJSDIR)/$(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf $(OBJSDIR)/$(BUILD) $(DESTDIR)/pRIblast.$(BUILD)

cleanall:
	rm -rf $(OBJSDIR) $(DESTDIR)
