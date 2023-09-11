#ifndef TP1_WRITER_H
#define TP1_WRITER_H

#include <string>
#include <fstream>
#include <iostream>

class Writer
{
    private:
        std::ofstream file;

    public:
        Writer(std::string name)
        {
            file.open (name);
            if(!file.is_open()){
                std::cout << "Writer - Open File Error: " << name << std::endl;
            }
        }
        ~Writer(){}

        void write(std::string line)
        {
                file << line + "\n";
        }

        void close()
        {
            file.close();
        }


};

#endif
