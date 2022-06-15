# belle2_dataobjects

gaudi_add_library(dataobjects
        SOURCES
            #src/BackgroundInfo.cc
            #src/BackgroundMetaData.cc
            #src/BinnedEventT0.cc
            #src/DisplayData.cc
            src/EventMetaData.cc
            #src/EventT0.cc
            #src/FileMetaData.cc
            #src/Helix.cc
            #src/MCInitialParticles.cc
            #src/ProfileInfo.cc
            #src/RelationContainer.cc
            #src/RelationElement.cc
            #src/UncertainHelix.cc
)

target_include_directories(dataobjects PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS dataobjects
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
