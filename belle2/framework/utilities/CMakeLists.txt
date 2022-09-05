# belle2_Futilities

gaudi_add_library(Futilities
        SOURCES
            src/Conversion.cc
            src/FileSystem.cc
            src/Utils.cc
        LINK
            ${ROOT_LIBRARIES}
            Python::Python
            Boost::boost
            logging
            boost_filesystem

)

target_include_directories(Futilities PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS Futilities
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
