# belle2_geometry

gaudi_add_library(geometry
        SOURCES
            src/Circle2D.cc
            src/GeneralizedCircle.cc
            src/Helix.cc
            src/HelixParameters.cc
            src/Line2D.cc
            src/LineParameters.cc
            src/ParameterLine2D.cc
            src/PerigeeCircle.cc
            src/PerigeeParameters.cc
            src/SZLine.cc
            src/SZParameters.cc
            src/UncertainHelix.cc
            src/UncertainParameterLine2D.cc
            src/UncertainPerigeeCircle.cc
            src/UncertainSZLine.cc
            src/Vector2D.cc
            src/Vector3D.cc
        LINK
            Eigen3::Eigen
            numerics
            ${ROOT_LIBRARIES}
)

target_include_directories(geometry PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS geometry
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
