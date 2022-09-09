# belle2_hits

gaudi_add_library(hits
        SOURCES
            src/CDCFacet.cc
            src/CDCRecoHit2D.cc
            src/CDCRecoHit3D.cc
            src/CDCRLWireHit.cc
            src/CDCRLWireHitPair.cc
            src/CDCRLWireHitTriple.cc
            src/CDCTangent.cc
            src/CDCWireHit.cc
        LINK ${ROOT_LIBRARIES}
            topology
            geometry
            trajectories
)

target_include_directories(hits PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS hits
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
