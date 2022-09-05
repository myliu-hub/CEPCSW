# belle2_core

gaudi_add_library(core
        SOURCES
            include/ModuleParam.templateDetails.h
            src/ModuleParamBase.cc
            #src/SubEventModule.cc
            #src/ModuleParamInfoPython.cc 
            src/ModuleParam.cc
            src/ModuleParamList.cc
        LINK
            Python::Python
            Boost::boost Boost::python
            gearbox
            pcore
            Futilities
            logging
)

target_include_directories(core PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS core
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
