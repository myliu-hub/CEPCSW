# belle2_topology

gaudi_add_library(topology
        SOURCES
            #src/CDCGeometryLoader.cc
            src/CDCWire.cc
            src/CDCWireLayer.cc
            src/CDCWireSuperLayer.cc
            src/CDCWireTopology.cc
            #src/EWirePositionToCDC.cc
            src/ISuperLayer.cc
            src/WireLine.cc
            src/WireNeighborPair.cc
)

target_include_directories(topology PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS topology
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
