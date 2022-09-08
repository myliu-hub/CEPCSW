# belle2_utilities

gaudi_add_library(utilities
        SOURCES
            src/CompositeProcessingSignalListener.cc
            #src/ModuleParamInstances.cc
            #src/ParameterVariant.cc
            src/ProcessingSignalListener.cc
        LINK
            Python::Python
            Boost::boost
)

target_include_directories(utilities PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS utilities
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
