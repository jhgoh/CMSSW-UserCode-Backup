from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
from PhysicsTools.PatAlgos.patEventContent_cff import patExtraAodEventContent

patDHEventContent = patEventContent+patExtraAodEventContent

# drop default cleanPatMuons and cleanPatElectrons
patDHEventContentWithLeptonFilter = patDHEventContent[:]
patDHEventContentWithLeptonFilter.remove('keep *_cleanPatElectrons_*_*')
patDHEventContentWithLeptonFilter.remove('keep *_cleanPatMuons_*_*')
patDHEventContentWithLeptonFilter.append('keep *_goodPatElectrons_*_*')
patDHEventContentWithLeptonFilter.append('keep *_goodPatMuons_*_*')
