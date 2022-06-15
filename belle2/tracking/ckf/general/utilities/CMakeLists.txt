# belle2_generalutilities

#if (GenFit_FOUND)
gaudi_add_library(generalutilities
        SOURCES
            src/Advancer.cc
        LINK
            GenFit::genfit2
            Eigen3::Eigen
)

target_include_directories(generalutilities PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

#target_compile_definitions(generalutilities PRIVATE GENFIT_MY_DEBUG)

install(TARGETS generalutilities
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)

#endif()
