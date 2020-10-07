#ifndef NUCLEOTIDES_H
#define NUCLEOTIDES_H

#include "../Ymir/proteins.h"

/// \file
/// \brief Manipulation of DNA sequences and transforming them into AA sequences.
/// \date 10th October 2019 \author Philippe A. Robert
/// \defgroup DNA Functions to using nucleotide sequences (nucleotides.h/cpp)

                    /// \brief Convert nucleotides into protein. Stop codonm becomes '!' \ingroup DNA
char getAA(string nucleotideSequence, int startPos = 0);
                    /// \brief Gives a random nucleotide \ingroup DNA
char randomNucleotide();
                    /// \brief Mutation. By default can be synonymous and induce stops, unless said opposite
string mutateDNA(string seq, bool onlyNonSynMuts = false, bool mutateAgainIfStop = false);
                    /// \brief Replaces stop codons in a nucleotide sequence by three other random nucleotides \ingroup DNA
string clearStops(string nucl);
                    /// \brief Creates a random DNA sequence, that might contain start/stop codons \ingroup DNA
string randomDNA(int size);
                    /// \brief Transforms a DNA sequence into protein from the beginning (no need of start codon). Input size should be multiple of 3. Stop codons become '!' \ingroup DNA
string convertToProtein(string nucl);
                    /// \brief Transforms a DNA sequence into protein from the beginning (no need of start codon).
                    /// \param nucle DNA sequence (string). Size should be multiple of 3.
                    /// \result AA sequence. Stop codons become '!'
                    /// \result flag containsStop is put to true if and only if a stop codon is found. \ingroup DNA
string convertToProtein(string nucl, bool &containsStop);

void testNucleotides();

#endif
