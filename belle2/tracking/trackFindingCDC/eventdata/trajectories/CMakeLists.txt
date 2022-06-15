# belle2_trajectories

gaudi_add_library(trajectories
        SOURCES
            src/CDCBFieldUtil.cc
            src/CDCTrajectory2D.cc
            src/CDCTrajectory3D.cc
            src/CDCTrajectorySZ.cc
        LINK
            GenFit::genfit2
            Eigen3::Eigen
)

target_include_directories(trajectories PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS trajectories
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
