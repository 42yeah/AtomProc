#include "AtomProc.h"

#ifdef _DEBUG
#define THROW throw
#else
#define THROW return
#endif

using namespace AP;


std::unique_ptr<AtomProc> AtomProc::make_atom_proc(std::string data_path, std::string bin_path) {
    std::unique_ptr<AtomProc> proc(new AtomProc());
    
    if (bin_path != "" && std::filesystem::exists(bin_path)) {
        std::ifstream reader(bin_path, std::ios::binary);
        if (!reader.good()) {
            std::cerr << "Cannot find " << bin_path << ". Resorting to reading .data file..." << std::endl;
            proc->parse_mat_from_data(data_path);
        }
        proc->atom_data_molecule = mat_load(reader);
        proc->atom_bond = mat_load(reader);
        proc->atom_angle = mat_load(reader);
        proc->atom_impropers = mat_load(reader);
        proc->atom_dihedrals = mat_load(reader);
        reader.close();
    } else {
        proc->parse_mat_from_data(data_path);

        std::ofstream writer(bin_path, std::ios::binary);
        mat_dump(writer, proc->atom_data_molecule);
        mat_dump(writer, proc->atom_bond);
        mat_dump(writer, proc->atom_angle);
        mat_dump(writer, proc->atom_impropers);
        mat_dump(writer, proc->atom_dihedrals);
        writer.close();
    }
    
    return std::move(proc);
}

bool AP::AtomProc::parse_mat_from_data(std::string data_file) {
    std::ifstream dump(data_file);
    if (!dump.good()) {
        std::cerr << "Dumpfile not found!" << std::endl;
        return false;
    }

    char line[10240] = { 0 };
    while (!dump.eof()) {
        dump.getline(line, sizeof(line));
        std::string token(line);
        if (starts_with(token, "Atoms")) {
            // skip a line here
            dump.getline(line, sizeof(line));
            for (int j = 0; j < Natoms; j++) {
                // %d %d %d %f %f %f %f %d %d %d
                for (int i = 0; i < 10; i++) { 
                    dump >> atom_data_molecule(j, i);
                }
            }
        } else if (starts_with(token, "Bonds")) {
            dump.getline(line, sizeof(line));
            for (int j = 0; j < Bond_number; j++) {
                // %d %d %f %f
                for (int i = 0; i < 4; i++) {
                    dump >> atom_bond(j, i);
                }
            }
        } else if (starts_with(token, "Angles")) {
            dump.getline(line, sizeof(line));            
            for (int j = 0; j < Angle_number; j++) {
                // %d %d %f %f %f
                for (int i = 0; i < 5; i++) {
                    dump >> atom_angle(j, i);
                }
            }
        } else if (starts_with(token, "Impropers")) {
            dump.getline(line, sizeof(line));            
            for (int j = 0; j < Impropers_number; j++) {
                // %d %d %f %f %f %f
                for (int i = 0; i < 6; i++) {
                    dump >> atom_impropers(j, i);
                }
            }
        } else if (starts_with(token, "Dihedrals")) {
            dump.getline(line, sizeof(line));            
            for (int j = 0; j < Dihedrals_number; j++) {
                // %d %d %f %f %f %X
                for (int i = 0; i < 5; i++) {
                    dump >> atom_dihedrals(j, i);
                }
                double trash;
                dump >> trash;
            }
        }
    }

    dump.close();
    return true;
}

void AP::AtomProc::random_sample() {
    std::cout << "Atoms: " << std::endl;
    ::random_sample(atom_data_molecule, 5);
    
    std::cout << "Bonds: " << std::endl;
    ::random_sample(atom_bond, 5);
    
    std::cout << "Angles: " << std::endl;
    ::random_sample(atom_angle, 5);
    
    std::cout << "Impropers: " << std::endl;
    ::random_sample(atom_impropers, 5);
    
    std::cout << "Dihedrals: " << std::endl;
    ::random_sample(atom_dihedrals, 5);
}

std::optional<std::string> AP::AtomProc::run() {
    std::vector<int> position = locate(atom_data_molecule.col(2), (double) delete_type);
    if (position.size() == 0) {
        THROW "Cannot find delete_type in atom_data_molecule";
    }

    Eigen::MatrixXd ring_carbon = atom_data_molecule(position, Eigen::all);
    
    // 将在一个苯环上的芳香碳原子放在连续六行
    Eigen::MatrixXd ring_carbon_new = Eigen::MatrixXd::Zero(length(ring_carbon), 13);
    int mol_number = (int) ring_carbon.col(1).maxCoeff();

    // Parallel probability?
    for (int i = 0; i < mol_number; i++) {
        std::vector<int> mol = locate(ring_carbon.col(1), (double) (i + 1));
        Eigen::MatrixXd ring_carbon_mol = ring_carbon(mol, Eigen::all);
        for (int j = 0; j < length(ring_carbon_mol); j++) {
            Eigen::MatrixXd atom_id = Eigen::MatrixXd::Zero(3, 1);
            std::vector<int> atom_1 = locate(atom_bond.col(2), ring_carbon_mol(j, 0));
            std::vector<int> atom_2 = locate(atom_bond.col(3), ring_carbon_mol(j, 0));
            Eigen::VectorXd atom_1_id = atom_bond(atom_1, 3);
            Eigen::VectorXd atom_2_id = atom_bond(atom_2, 2);

            if (atom_1_id.size() + atom_2_id.size() == 2) {
                int tot_r = 0;
                for (int r = 0; r < atom_1_id.size(); r++) {
                    atom_id(tot_r++, 0) = atom_1_id(r);
                }
                for (int r = 0; r < atom_2_id.size(); r++) {
                    atom_id(tot_r++, 0) = atom_2_id(r);
                }
            } else {
                atom_id.resize(atom_1_id.size() + atom_2_id.size(), 1);
                int tot_r = 0;
                for (int r = 0; r < atom_1_id.size(); r++) {
                    atom_id(tot_r++, 0) = atom_1_id(r);
                }
                for (int r = 0; r < atom_2_id.size(); r++) {
                    atom_id(tot_r++, 0) = atom_2_id(r);
                }
            }
            ring_carbon_new.row(i * length(ring_carbon_mol) + j) << ring_carbon_mol.row(j), atom_id.transpose();
        }
    }


    return std::nullopt;
}

AtomProc::AtomProc() {
    atom_masses.col(0) = Eigen::VectorXd::LinSpaced(Natoms, 0, (double) Natoms - 1);
    atom_masses.col(1).array() = 12.0;
}


void AP::mat_dump(std::ofstream& out, Eigen::MatrixXd& mat) {
    int rows = (int) mat.rows(), cols = (int) mat.cols();
    out.put('a');
    out.write((const char *) &rows, sizeof(int));
    out.write((const char *) &cols, sizeof(int));
    out.write((const char *) &mat(0, 0), (long long) rows * cols * sizeof(double));
}

Eigen::MatrixXd AP::mat_load(std::ifstream& in) {
    Eigen::MatrixXd ret;
    int rows, cols;
    char align = in.get();
    if (align != 'a') {
        std::cerr << "ERR! Unalignment found." << std::endl;
        return ret;
    }
    in.read((char *) &rows, sizeof(int));
    in.read((char *) &cols, sizeof(int));
    ret.resize(rows, cols);
    in.read((char *) &ret(0, 0), (long long) rows * cols * sizeof(double));
    return ret;
}

void AP::random_sample(Eigen::MatrixXd mat, int n_samples) {
    std::uniform_int_distribution<int> row_distrib(0, (int) mat.rows() - 1);
    std::uniform_int_distribution<int> col_distrib(0, (int) mat.cols() - 1);
    std::random_device dev;
    for (int i = 0; i < n_samples; i++) {
        int row = row_distrib(dev), col = col_distrib(dev);
        std::cout << "Sample: (" << row << ", " << col << "): " << mat(row, col) << std::endl;
    }
    std::cout << std::endl;
}

int AP::length(Eigen::MatrixXd mat) {
    return (int) (mat.rows() > mat.cols() ? mat.rows() : mat.cols());
}

std::vector<int> AP::locate(Eigen::MatrixXd vec, double value) {
    std::vector<int> ret;
    for (int i = 0; i < vec.size(); i++) {
        if (vec(i) == value) {
            ret.push_back(i);
        }
    }
    return ret;
}

bool AP::starts_with(std::string what, std::string starts) {
    if (what.substr(0, starts.size()) == starts) {
        return true;
    }
    return false;
}
