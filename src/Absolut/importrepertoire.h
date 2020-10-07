#ifndef IMPORTREPERTOIRE_H
#define IMPORTREPERTOIRE_H

#include <vector>
#include <string>

//std::vector<std::string> slides(std::string largeSeq, int size);

// there is something wrong in this definition that crashes the compiling. don't know what.
//void analyzeRepertoire(std::string AntigenName, int receptorSize, int minInteract, int nParallel = 1,
//                       bool restartFromScratch = false, int IDprocess = 0);


class importRepertoire{
public:
    importRepertoire(){}
    importRepertoire(std::string fileName, bool thirdColumn = false);

    // Adding a sequence, like for merging repertoires together, it will check if the same ID or same sequence exist.
    // Adds a sequence. If a sequence is given whose ID already exists, tests the sequence. If sequence the same, do not add. If sequence is
    // different, then raises an error. If ID=-1, will find a new ID for it (caution, might collide with IDs of added sequences later that would have
    // this ID. If mergeIDsWhenSequenceIsEqual, a sequence already existing will not be added even if it has a non-used ID.
    void addSequence(std::string sequence, int ID = -1, bool treated = false, bool mergeIDsWhenSequenceIsEqual = false);
    std::vector<int> listIDs;
    std::vector<std::string> sequences;
    std::vector<bool> alreadyTreated;
    std::string nextSequence(bool setTreated = true, int lastLine = 1e10);
    std::pair<int, std::string> nextSequenceAndID(bool setTreated, int lastLine = 1e10);
    bool writeUpdatedRepertoire(std::string fileName = "");
    int nLines;

    int firstBlockedLine;
    int lastBlockedLine;
//private:
    std::string originalFileName;
    size_t currentIndex;
};

void testImportRepertoire();

#endif // IMPORTREPERTOIRE_H
