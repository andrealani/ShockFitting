##############################################################################
# include cmake macros
##############################################################################
INCLUDE(CheckIncludeFile)
INCLUDE(CheckIncludeFileCXX)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckSymbolExists)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckLibraryExists)
INCLUDE(CheckTypeSize)
INCLUDE(CheckCSourceCompiles)
INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckCCompilerFlag)
INCLUDE(CheckCXXCompilerFlag)


##############################################################################
# include ShockFitting macros
##############################################################################

INCLUDE(macros/SFVariables)
INCLUDE(macros/SFListOperations)
INCLUDE(macros/SFSearchPaths)
INCLUDE(macros/SFOptionAddSubdirectory)
INCLUDE(macros/SFSetUnion)
INCLUDE(macros/SFLogToFile)
INCLUDE(macros/SFDumpVariables)
INCLUDE(macros/SFBoolTo01)
INCLUDE(macros/SFSeparateSources)
INCLUDE(macros/SFAddLibrary)
INCLUDE(macros/SFAddKernelLibrary)
INCLUDE(macros/SFAddPluginLibrary)
INCLUDE(macros/SFAddPluginApp)
INCLUDE(macros/SFWarnOrphanFiles)
INCLUDE(macros/SFCheckFileLength)
INCLUDE(macros/SFAddCompilationFlags)
