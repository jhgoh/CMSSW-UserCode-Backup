from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
from PhysicsTools.PatAlgos.patEventContent_cff import patExtraAodEventContent

patDHEventContent = patEventContent+patExtraAodEventContent

# drop default cleanPatMuons and cleanPatElectrons
patDHEventContentWithLeptonFilter = patDHEventContent[:]
patDHEventContentWithLeptonFilter.remove('keep *_cleanPatElectrons*_*_*')
patDHEventContentWithLeptonFilter.remove('keep *_cleanPatMuons*_*_*')
patDHEventContentWithLeptonFilter.append('keep *_goodPatElectrons*_*_*')
patDHEventContentWithLeptonFilter.append('keep *_goodPatMuons*_*_*')
