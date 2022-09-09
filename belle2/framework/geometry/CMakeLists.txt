# belle2_frameworkgeometry

gaudi_add_library(frameworkgeometry
        SOURCES
            src/BFieldManager.cc
        LINK ${ROOT_LIBRARIES}
)

target_include_directories(frameworkgeometry PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS frameworkgeometry
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
