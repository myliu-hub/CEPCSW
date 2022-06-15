# belle2_utilities

gaudi_add_module(utilities
        SOURCES
            src/ColorPalette.cc
            src/Conversion.cc
            #src/EnvironmentVariables.cc
            src/FileSystem.cc
            src/FormulaParser.cc
            src/GeneralCut.cc
            src/HTML.cc
            src/IOIntercept.cc
            src/KeyValuePrinter.cc
            src/MakeROOTCompatible.cc
            #src/MultivariateNormalGenerator.cc
            src/NumberSequence.cc
            src/RegisterPythonModule.cc
            src/RootFileCreationManager.cc
            src/Spline.cc
            src/Stream.cc
            src/TestHelpers.cc
            src/Utils.cc
)

target_include_directories(utilities PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS utilities
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
