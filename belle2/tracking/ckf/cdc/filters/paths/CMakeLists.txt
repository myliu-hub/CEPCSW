# belle2_paths

gaudi_add_library(paths
        SOURCES
            src/CDCPathBasicVarSet.cc
            src/CDCPathFilterFactory.cc
            src/CDCPathTruthVarSet.cc
)

target_include_directories(paths PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../../../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS paths
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
