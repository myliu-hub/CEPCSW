# belle2_cdcentities

if (GenFit_FOUND)
gaudi_add_library(cdcentities
        SOURCES
            src/CDCCKFPath.cc
            src/CDCCKFState.cc
        LINK GenFit::genfit2
)

target_include_directories(cdcentities PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_definitions(cdcentities PRIVATE GENFIT_MY_DEBUG)

install(TARGETS cdcentities
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
endif()
