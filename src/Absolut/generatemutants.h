#ifndef GENERATEMUTANTS_H
#define GENERATEMUTANTS_H

#include <vector>
#include <string>
#include <map>
#include <algorithm>
using namespace std;

// This class takes two classes of sequences:
// - binders
// - nonbinders or anything else

// => aims:
//      - generate sequences with at least all mutations in each position separately + random mutations
//      - aims to generate mutants "In-between" the two pools of sequences (to go from one sequence to the other one) (cloud)
//              aims to find sequences with similar patterns in both pools,
class generateMutants
{
public:
    generateMutants();
};

/*
 Questions:

 Level 0: what is the vicinity of binding on the antibody side?

 generate mutants with at least all point-mutations tested, then random secondary or third mutations
 analysis:
      - how many mutations to exit the binding class we are in?
      - for each position, how far do we stay in a ncpr class?
      - among still binders, how far do we stay in the same structure
      - among still binders, how far do we stay in the same paratope/epitope (well, similar to structure)
      We could analyze vicinities of different affinity level, then different patterns

 Level 1: what is the vicinity of binding on the antigen side

 A priori, only the positions on existing hotspots are interesting, although we could want to make a new immunogenic antigen

 Mutate one position on the antigen (only on one hotspot), then two? which AAs do we test?
 note, maybe the structures become different.

 analysis is symmetrical.

 Level 2: linkage desequilibrium and co-evolution

 Question: if one mutation on one side, can we detect compensatory mutaations on the other sie
 (potentially to keep binding on the same place)

 idea: take mutated antigen, then "refine" the binders to become higher affinity or what?

 is a pattern predictive?
 Best would be:

 Level 3: Interpretability:

 3a: are motifs really the rule?

 take a type of motif, and generate sequences inside and outside the motif.
 problem: there are different ways to bind, so this is testing the vicinity of one motif?

 If we want to predict cross-reactivity, this is a bit different,



 1- cluster patterns of similarity (get all k-mers?)

 2 - from a motif, can we test if it is true?

 */
//
#endif // GENERATEMUTANTS_H
