# belle2_pcore

gaudi_add_module(pcore
        SOURCES
            src/AsyncWrapper.cc
            src/DataStoreStreamer.cc
            src/EvtMessage.cc
            src/GlobalProcHandler.cc
            src/Mergeable.cc
            src/MsgHandler.cc
            #src/PathUtils.cc
            src/pEventProcessor.cc
            #src/ProcessMonitor.cc
            src/ProcHandler.cc
            src/RbTuple.cc
            src/RingBuffer.cc
            src/RxModule.cc
            src/SemaphoreLocker.cc
            src/SeqFile.cc
            src/TxModule.cc
            #src/ZMQEventProcessor.cc
)

target_include_directories(pcore PUBLIC
        #$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../../include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS pcore
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
