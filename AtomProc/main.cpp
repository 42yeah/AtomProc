// by 42yeah

#include <iostream>
#include <fstream>
#include <random>
#include "AtomProc.h"



int main() {
    std::string file = "packing-PBT-101_773repeated-units-PBT-Cu-medium-voids--longFilling-1800atm-packing-800-Cu-0K-2-20ps-2-No-Cu-H_O_1.data";
    std::string bindump_file = "atomproc.bin";

    std::unique_ptr<AP::AtomProc> proc = AP::AtomProc::make_atom_proc(file, bindump_file);
    if (proc == nullptr) {
        std::cerr << "AtomProc has failed to initialize. Please see the reason presented above." << std::endl;
        return -1;
    }
    proc->run();

    return 0;
}
