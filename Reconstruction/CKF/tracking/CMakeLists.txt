# belle2_tracking

gaudi_add_library(tracking
        SOURCES
            src/ckf/Advancer.cc
            src/ckf/RecoTrack.cc
            src/ckf/CDCCKFPath.cc
            src/ckf/CDCCKFResultStorer.cc
            src/ckf/CDCCKFState.cc
            src/ckf/CDCPathBasicVarSet.cc
            src/ckf/CDCPathFilterFactory.cc
            src/ckf/CDCPathPairFilterFactory.cc
            src/ckf/CDCPathTruthVarSet.cc
            src/ckf/CDCStateBasicVarSet.cc
            src/ckf/CDCStateFilterFactory.cc
            src/ckf/CKFToCDCFindlet.cc
            src/ckf/DistanceCDCStateFilter.cc
            src/ckf/ExtrapolateAndUpdateCDCStateFilter.cc
            src/ckf/RoughCDCStateFilter.cc
            src/ca/AutomatonCell.cc
            src/eventdata/CDCFacet.cc
            src/eventdata/CDCRecoHit2D.cc
            src/eventdata/CDCRecoHit3D.cc
            src/eventdata/CDCRLWireHit.cc
            src/eventdata/CDCRLWireHitPair.cc
            src/eventdata/CDCRLWireHitTriple.cc
            src/eventdata/CDCTangent.cc
            src/eventdata/CDCWireHit.cc
            src/eventdata/CDCBFieldUtil.cc
            src/eventdata/CDCTrajectory2D.cc
            src/eventdata/CDCTrajectory3D.cc
            src/eventdata/CDCTrajectorySZ.cc
            src/geometry/Circle2D.cc
            src/geometry/GeneralizedCircle.cc
            src/geometry/Helix.cc
            src/geometry/HelixParameters.cc
            src/geometry/Line2D.cc
            src/geometry/LineParameters.cc
            src/geometry/ParameterLine2D.cc
            src/geometry/PerigeeCircle.cc
            src/geometry/PerigeeParameters.cc
            src/geometry/SZLine.cc
            src/geometry/SZParameters.cc
            src/geometry/UncertainHelix.cc
            src/geometry/UncertainParameterLine2D.cc
            src/geometry/UncertainPerigeeCircle.cc
            src/geometry/UncertainSZLine.cc
            src/geometry/Vector2D.cc
            src/geometry/Vector3D.cc
            src/filters/Filter.cc
            src/filters/FilterFactory.cc
            src/mva/PyEstimator.cc
            src/numerics/Median.cc
            src/numerics/SinEqLine.cc
            src/numerics/SpecialFunctions.cc
            src/numerics/Weight.cc
            src/topology/CDCWire.cc
            src/topology/CDCWireLayer.cc
            src/topology/CDCWireSuperLayer.cc
            src/topology/CDCWireTopology.cc
            src/topology/ISuperLayer.cc
            src/topology/WireLine.cc
            src/topology/WireNeighborPair.cc
            src/utilities/CompositeProcessingSignalListener.cc
            src/utilities/ProcessingSignalListener.cc
            src/varsets/NamedFloatTuple.cc
        LINK
            Eigen3::Eigen
            ${ROOT_LIBRARIES}
            Boost::boost
            GenFit::genfit2
            Python::Python
            cdcdataobjects
            boost_python
            framework
)

target_include_directories(tracking PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/ca
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/eventdata
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/filters
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/findlets
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/geometry
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/mva
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/numerics
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/topology
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/utilities
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/varsets
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include/ckf
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS tracking
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
