# Compiler
CXX := g++

# Sources / headers / target
SRCS := ./src/RKR_gui_top.cpp ./src/calculation_core.cpp ./src/calc_wrapper.hpp
HDRS := ./src/calculation_core.hpp  ./src/calc_wrapper.hpp ./src/SplineNaK.h ./src/SplineNaK2.h
OUT  := RKR_v1.exe
# wxWidgets (static)
WX_CXXFLAGS := $(shell wx-config --static --cxxflags)
WX_LIBS     := $(shell wx-config --static --libs)

# ---- Your local GSL ----
GSL_INC   := ./lib/gsl/include
GSL_LIB   := ./lib/gsl/lib

# ---- Header-only libs ----
#CSPL_INC  := ./lib/nak-spline/include
EIGEN_DIR := ./lib/eigen-master

# ============================================================
STATIC_WINPTHREAD ?= 1
FULLY_STATIC_EXE ?= 0
# ============================================================

# Include paths: header-only libs just need -I
INCLUDES := -I$(EIGEN_DIR) -I$(GSL_INC)  #-I$(CSPL_INC)

# Compile flags
CXXFLAGS := $(WX_CXXFLAGS) -O2 -Wall -pthread -fopenmp $(INCLUDES)

# Base linker/runtime static flags
BASE_STATIC_FLAGS := -static-libstdc++ -static-libgcc
ifeq ($(FULLY_STATIC_EXE),1)
  BASE_STATIC_FLAGS := -static -static-libstdc++ -static-libgcc
endif

# Libraries (WITH STATIC GSL AND OpenMP)
LDLIBS := $(WX_LIBS) \
          -Wl,--start-group -Wl,-Bstatic -L$(GSL_LIB) -lgsl -lgslcblas -Wl,-Bdynamic -lm -Wl,--end-group \
          -Wl,-Bstatic -lgomp -Wl,-Bdynamic

ifeq ($(STATIC_WINPTHREAD),1)
  LDLIBS += -Wl,-Bstatic -lstdc++ -lpthread -Wl,-Bdynamic
endif

.PHONY: all clean strip info

all: $(OUT)

$(OUT): $(SRCS) $(HDRS)
	$(CXX) $(SRCS) -o $(OUT) $(CXXFLAGS) $(BASE_STATIC_FLAGS) $(LDLIBS)

strip: $(OUT)
	strip $(OUT)

info:
	@echo "CXXFLAGS:" && echo $(CXXFLAGS)
	@echo "LDLIBS:" && echo $(LDLIBS)
	@echo "SRCS:" && echo $(SRCS)

clean:
	$(RM) $(OUT)
