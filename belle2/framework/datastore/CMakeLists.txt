# belle2_datastore

gaudi_add_module(datastore
        SOURCES
            src/DataStore.cc
            src/DependencyMap.cc
            src/RelationEntry.cc
            src/RelationIndexManager.cc
            src/RelationsObject.cc
            src/RelationVector.cc
            src/SelectSubset.cc
            src/StoreAccessorBase.cc
            src/StoreArray.cc
            src/StoreEntry.cc
            src/StoreObjPtr.cc
)

target_include_directories(datastore PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS datastore
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
