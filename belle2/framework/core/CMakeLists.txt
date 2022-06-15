# belle2_core

gaudi_add_module(core
        SOURCES
            src/DataFlowVisualization.cc
            src/Environment.cc
            src/EventProcessor.cc
            src/FileCatalog.cc
            src/InputController.cc
            src/MetadataService.cc
            src/Module.cc
            src/ModuleCondition.cc
            src/ModuleManager.cc
            src/ModuleParamBase.cc
            src/ModuleParam.cc
            #src/ModuleParamInfoPython.cc
            src/ModuleParamList.cc
            src/Path.cc
            src/ProcessStatistics.cc
            src/RandomGenerator.cc
            #src/RandomNumbers.cc
            src/SubEventModule.cc
            src/SwitchDataStoreModule.cc
        LINK
            Python::Python
            Boost::boost
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
