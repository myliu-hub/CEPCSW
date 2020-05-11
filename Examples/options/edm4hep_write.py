#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import CEPCDataSvc
dsvc = CEPCDataSvc("EventDataSvc")

from Configurables import Edm4hepWriteAlg
alg = Edm4hepWriteAlg("Edm4hepWriteAlg")
alg.HeaderCol.Path = "EventHeader"
alg.OutputCol.Path = "MCParticle"

from Configurables import PodioOutput
out = PodioOutput("out")
out.filename = "test.root"
out.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [alg, out],
                EvtSel = 'NONE',
                EvtMax = 10,
                ExtSvc=[dsvc],
                OutputLevel=DEBUG
)
