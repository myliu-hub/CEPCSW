# RecGenfitAlg

if (GenFit_FOUND)
gaudi_add_module(RecGenfitAlg
        SOURCES src/RecGenfitAlgSDT.cpp
                src/GenfitTrack.cpp
                src/GenfitHit.cpp
                src/GenfitField.cpp
                src/GenfitFitter.cpp
                src/GenfitMaterialInterface.cpp
                src/LSFitting.cpp
                src/WireMeasurementDC
                src/PlanarMeasurementSDT
                src/MBg.cpp
                src/DCTrackFinding.cpp
          LINK GearSvc
             Gaudi::GaudiAlgLib
             Gaudi::GaudiKernel
             ${GEAR_LIBRARIES}
             ${GSL_LIBRARIES}
             ${LCIO_LIBRARIES}
             DetSegmentation
             DetInterface
             DataHelperLib
             EDM4HEP::edm4hep
             EDM4HEP::edm4hepDict
             GenFit::genfit2
             cdcentities
             hits
             cdcdataobjects
             findlets
             pcore
             Boost::boost Boost::python
             #datastore
)

target_include_directories(RecGenfitAlg PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        
)

target_compile_definitions(RecGenfitAlg PRIVATE GENFIT_MY_DEBUG)

install(TARGETS RecGenfitAlg
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
endif()
