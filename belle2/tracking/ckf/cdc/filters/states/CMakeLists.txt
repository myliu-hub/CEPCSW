# belle2_states

gaudi_add_library(states
        SOURCES
            src/CDCfromEclStateTruthVarSet.cc
            src/CDCStateBasicVarSet.cc
            src/CDCStateFilterFactory.cc
            src/DistanceCDCStateFilter.cc
            src/ExtrapolateAndUpdateCDCStateFilter.cc
            #src/MCTruthCDCStateFilter.cc
            #src/MCTruthEclSeedFilter.cc
            #src/RoughCDCfromEclStateFilter.cc
            src/RoughCDCStateFilter.cc
        LINK
            GenFit::genfit2
            Eigen3::Eigen
)

target_include_directories(states PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../../../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS states
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
