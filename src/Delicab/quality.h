#ifndef QUALITY_H
#define QUALITY_H

#include <string>
void evaluateAntigen(std::string  antigenID, int nrAGs, int nrBCRperAG, bool startFromRandomAAseq, int NbMuts, std::string maskMut, std::string initialSeq = std::string(""));

void testQuality();
#endif // QUALITY_H
