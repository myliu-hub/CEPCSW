# belle2_findlets

gaudi_add_library(findlets
        SOURCES
            #src/CDCCKFEclSeedCreator.cc
            src/CDCCKFResultStorer.cc
            src/CKFToCDCFindlet.cc
            #src/CKFToCDCFromEclFindlet.cc
        LINK utilities
             ${ROOT_LIBRARIES}
             GenFit::genfit2
             trajectories
             topology
             pathPairs
             geometry
             trackingdataobjects
             states
             paths
             cdcdataobjects
)

target_include_directories(findlets PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS findlets
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
