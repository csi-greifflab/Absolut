#ifndef HTML_H
#define HTML_H

#include <string>

std::string getPDB(std::string AGname);

std::string getChain(std::string AGname);

class html
{
    std::string AGname;
public:
    html(std::string _AGname);
    std::string analyze(std::string outputFoler = "", std::string folderWithBindings = "", std::string folderWithStructures = "");

    static std::string getHtmlHeader();
    static std::string getHtmlFooter();
};

#endif // HTML_H
