# belle2_logging

gaudi_add_library(logging
        SOURCES
            src/LogConfig.cc
            #src/LogConnectionConsole.cc
            src/LogConnectionFilter.cc
            src/LogConnectionJSON.cc
            src/LogConnectionTxtFile.cc
            src/LogConnectionUDP.cc
            src/LogMessage.cc
            src/LogMethod.cc
            src/LogSystem.cc
            src/LogVariableStream.cc
)

target_include_directories(logging PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS logging
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
