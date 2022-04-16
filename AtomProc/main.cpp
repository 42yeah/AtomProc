// by 42yeah

#include <iostream>
#include <fstream>
#include <random>
#include "AtomProc.h"

#define BINDUMP_PATH "atomproc.bin"



int main(int argc, char *argv[]) {
    if (argc < 2 || argc >= 4) {
        std::cerr << "Usage: " << argv[0] << " <input .data file> [output .data file]" << std::endl;
        return 0;
    }

    std::string file = argv[1];
    std::string output = "result.data";
    if (argc == 3) {
        output = argv[2];
    }

    std::unique_ptr<AP::AtomProc> proc = AP::AtomProc::make_atom_proc(file, BINDUMP_PATH);
    if (proc == nullptr) {
        std::cerr << "AtomProc has failed to initialize. Please see the reason presented above." << std::endl;
        return -1;
    }
    std::optional<std::string> res = proc->run();
    if (res) {
        std::cout << "ERR! Error encountered while running: " << *res << std::endl;
        return -1;
    }
    res = proc->dump_result_to_data("header.txt", output);
    if (res) {
        std::cout << "ERR! Error encountered while dumping: " << *res << std::endl;
        return -2;
    }

    return 0;
}
