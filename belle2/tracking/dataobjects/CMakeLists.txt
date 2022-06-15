# belle2_trackingdataobjects

if (GenFit_FOUND)
gaudi_add_library(trackingdataobjects
        SOURCES
            src/RecoTrack.cc
        LINK GenFit::genfit2
            Python::Python
            Boost::boost
)

target_include_directories(trackingdataobjects PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_definitions(trackingdataobjects PRIVATE GENFIT_MY_DEBUG)

install(TARGETS trackingdataobjects
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
endif()
