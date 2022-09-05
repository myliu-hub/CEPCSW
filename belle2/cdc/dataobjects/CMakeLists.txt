# belle2_cdcdataobjects

if (GenFit_FOUND)
gaudi_add_library(cdcdataobjects
        SOURCES
            src/CDCHit.cc
            #src/CDCRawHit.cc
            src/CDCRawHitWaveForm.cc
            src/CDCRecoHit.cc
            src/WireID.cc
        LINK GenFit::genfit2
        logging
)

target_include_directories(cdcdataobjects PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_definitions(cdcdataobjects PRIVATE GENFIT_MY_DEBUG)

install(TARGETS cdcdataobjects
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
endif()
