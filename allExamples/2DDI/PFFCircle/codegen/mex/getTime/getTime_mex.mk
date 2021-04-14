MATLAB_ROOT = E:\MATLAB\R2018b
MAKEFILE = getTime_mex.mk

include getTime_mex.mki


SRC_FILES =  \
	matlab_c_get_current_time_in_sec.c \
	getTime_data.c \
	getTime_initialize.c \
	getTime_terminate.c \
	getTime.c \
	_coder_getTime_info.c \
	_coder_getTime_api.c \
	_coder_getTime_mex.c \
	c_mexapi_version.c

MEX_FILE_NAME_WO_EXT = getTime_mex
MEX_FILE_NAME = $(MEX_FILE_NAME_WO_EXT).mexw64
TARGET = $(MEX_FILE_NAME)

SYS_LIBS = 


#
#====================================================================
# gmake makefile fragment for building MEX functions using MinGW
# Copyright 2015-2017 The MathWorks, Inc.
#====================================================================
#

SHELL = cmd
LD = $(LINKER)
OBJEXT = o
.SUFFIXES: .$(OBJEXT)

OBJLISTC = $(SRC_FILES:.c=.$(OBJEXT))
OBJLISTCPP  = $(OBJLISTC:.cpp=.$(OBJEXT))
OBJLIST  = $(OBJLISTCPP:.cu=.$(OBJEXT))

target: $(TARGET)

ML_INCLUDES = -I "$(MATLAB_ROOT)/simulink/include"
ML_INCLUDES+= -I "$(MATLAB_ROOT)/toolbox/shared/simtargets"
SYS_INCLUDE = $(ML_INCLUDES)

# Additional includes

SYS_INCLUDE += -I "D:\Documents\GitHub\kinoRRTstar\allExamples\2DDI\PFFCircle\codegen\mex\getTime"
SYS_INCLUDE += -I "D:\Documents\GitHub\kinoRRTstar\allExamples\2DDI\PFFCircle"
SYS_INCLUDE += -I ".\interface"
SYS_INCLUDE += -I "$(MATLAB_ROOT)\extern\include"
SYS_INCLUDE += -I "."

EML_LIBS = -llibemlrt -llibcovrt -llibut -llibmwmathutil 
SYS_LIBS += $(CLIBS) $(EML_LIBS)

EXPORTFILE = $(MEX_FILE_NAME_WO_EXT)_mex.map
EXPORTOPT = -Wl,--version-script,$(EXPORTFILE)
LINK_FLAGS = $(filter-out /export:mexFunction, $(LINKFLAGS))
COMP_FLAGS = $(CFLAGS) $(OMPFLAGS) -D__USE_MINGW_ANSI_STDIO=1
CXX_FLAGS = $(CXXFLAGS) $(OMPFLAGS) -D__USE_MINGW_ANSI_STDIO=1
LINK_FLAGS = $(LINKFLAGS) 
LINK_FLAGS += $(OMPLINKFLAGS)
ifeq ($(EMC_CONFIG),optim)
  COMP_FLAGS += $(OPTIMFLAGS)
  CXX_FLAGS += $(OPTIMFLAGS)
  LINK_FLAGS += $(LINKOPTIMFLAGS)
else
  COMP_FLAGS += $(DEBUGFLAGS)
  CXX_FLAGS += $(DEBUGFLAGS)
  LINK_FLAGS += $(LINKDEBUGFLAGS)
endif
LINK_FLAGS += -o $(TARGET)
LINK_FLAGS += 

CCFLAGS = $(COMP_FLAGS) -std=c99   $(USER_INCLUDE) $(SYS_INCLUDE)
CPPFLAGS = $(CXX_FLAGS) -std=c++11   $(USER_INCLUDE) $(SYS_INCLUDE)

%.$(OBJEXT) : %.c
	$(CC) $(CCFLAGS) "$<"

%.$(OBJEXT) : %.cpp
	$(CXX) $(CPPFLAGS) "$<"

# Additional sources

%.$(OBJEXT) : /%.c
	$(CC) $(CCFLAGS) "$<"

%.$(OBJEXT) : D:\Documents\GitHub\kinoRRTstar\allExamples\2DDI\PFFCircle/%.c
	$(CC) $(CCFLAGS) "$<"

%.$(OBJEXT) : D:\Documents\GitHub\kinoRRTstar\allExamples\2DDI\PFFCircle\codegen\mex\getTime/%.c
	$(CC) $(CCFLAGS) "$<"

%.$(OBJEXT) : interface/%.c
	$(CC) $(CCFLAGS) "$<"



%.$(OBJEXT) : /%.cpp
	$(CXX) $(CPPFLAGS) "$<"

%.$(OBJEXT) : D:\Documents\GitHub\kinoRRTstar\allExamples\2DDI\PFFCircle/%.cpp
	$(CXX) $(CPPFLAGS) "$<"

%.$(OBJEXT) : D:\Documents\GitHub\kinoRRTstar\allExamples\2DDI\PFFCircle\codegen\mex\getTime/%.cpp
	$(CXX) $(CPPFLAGS) "$<"

%.$(OBJEXT) : interface/%.cpp
	$(CXX) $(CPPFLAGS) "$<"



%.$(OBJEXT) : /%.cu
	$(CC) $(CCFLAGS) "$<"

%.$(OBJEXT) : D:\Documents\GitHub\kinoRRTstar\allExamples\2DDI\PFFCircle/%.cu
	$(CC) $(CCFLAGS) "$<"

%.$(OBJEXT) : D:\Documents\GitHub\kinoRRTstar\allExamples\2DDI\PFFCircle\codegen\mex\getTime/%.cu
	$(CC) $(CCFLAGS) "$<"

%.$(OBJEXT) : interface/%.cu
	$(CC) $(CCFLAGS) "$<"




$(TARGET): $(OBJLIST) $(MAKEFILE)
	$(LD) $(EXPORTOPT) $(OBJLIST) $(LINK_FLAGS) $(SYS_LIBS)
	@cmd /C "echo Build completed using compiler $(EMC_COMPILER)"

#====================================================================
