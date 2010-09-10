#include "FWCore/Framework/interface/MakerMacros.h"

#include "HiggsAnalysis/DoublyChargedHiggs/plugins/DileptonProducer.h"
#include "HiggsAnalysis/DoublyChargedHiggs/plugins/DHGenEventFilter.h"
#include "HiggsAnalysis/DoublyChargedHiggs/plugins/MultipleCandCounterFilter.h"

DEFINE_FWK_MODULE(DileptonProducer);
DEFINE_FWK_MODULE(DHGenEventFilter);
DEFINE_FWK_MODULE(MultipleCandCounterFilter);

