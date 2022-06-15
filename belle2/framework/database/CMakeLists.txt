# belle2_database

gaudi_add_module(database
        SOURCES
            src/CentralMetadataProvider.cc
            #src/Configuration.cc
            #src/Database.cc
            src/DBImportBase.cc
            src/DBStore.cc
            src/DBStoreEntry.cc
            src/Downloader.cc
            src/EventDependency.cc
            src/IntervalOfValidity.cc
            src/LocalMetadataProvider.cc
            src/MetadataProvider.cc
            src/PayloadFile.cc
            src/PayloadProvider.cc
            src/TestingPayloadStorage.cc
)

target_include_directories(database PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS database
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
