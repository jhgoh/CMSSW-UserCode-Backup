import FWCore.ParameterSet.Config as cms

from CondCore.DBCommon.CondDBSetup_cfi import *
from DQMServices.Core.DQM_cfg import *
from DQMServices.Components.DQMEnvironment_cfi import *
from DQMServices.Components.DQMStoreStats_cfi import *

DQMStore.verbose = 0
DQM.collectorHost = ''
dqmSaver.convention = 'Offline'
dqmSaver.workflow = '/DQM/Offline/FourVector'

from DQMOffline.Trigger.FourVectorHLTOffline_cfi import *
from DQMOffline.Trigger.FourVectorHLTOfflineClient_cfi import *

