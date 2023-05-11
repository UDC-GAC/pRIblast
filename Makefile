BUILD=release

CXX=mpic++
CXXFLAGS+=-std=c++17 -Wall -Wextra -Wpedantic -Wno-unused-variable -Wno-cast-function-type

cxxflags.release=-fopenmp -flto -O3 -march=native -s

CXXFLAGS+=$(cxxflags.$(BUILD))

LD=mpic++

ldflags.release=-fopenmp -flto=auto

LDFLAGS+=$(ldflags.$(BUILD))

SRCSDIR=src
OBJSDIR=obj
DESTDIR=target

SRCS=$(wildcard $(SRCSDIR)/*.cpp)
OBJS=$(patsubst $(SRCSDIR)/%.cpp,$(OBJSDIR)/$(BUILD)/%.o,$(SRCS))
DEPS=$(patsubst $(SRCSDIR)/%.cpp,$(OBJSDIR)/$(BUILD)/%.d,$(SRCS))

pRIblast.$(BUILD): $(OBJS)
	@mkdir -p $(DESTDIR)
	$(LD) $(LDFLAGS) -o $(DESTDIR)/$@ $^ $(LDLIBS)

-include $(DEPS)

$(OBJSDIR)/$(BUILD)/%.o: $(SRCSDIR)/%.cpp
	@mkdir -p $(OBJSDIR)/$(BUILD)
	$(CXX) $(CXXFLAGS) -MMD -MP $(INCLUDES) -c $< -o $@

c: clean
clean:
	rm -rf $(OBJSDIR) $(DESTDIR)
