# belle2_framework

gaudi_add_library(framework
        SOURCES
            src/geometry/BFieldManager.cc
            src/pcore/Mergeable.cc
            src/utilities/Conversion.cc
            src/utilities/FileSystem.cc
            src/utilities/Utils.cc
        LINK ${ROOT_LIBRARIES}
            Python::Python
            Boost::boost
            Boost::filesystem
)

target_include_directories(framework PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/core
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/gearbox
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/geometry
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/pcore
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/utilities
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS framework
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
