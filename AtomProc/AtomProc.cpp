#include "AtomProc.h"

#include <filesystem>
#include <thread>
#include <Eigen/Dense>

#ifdef _DEBUG
#define THROW throw
#else
#define THROW return
#endif

#define RING_CARBON_NEW_DUMP_PATH ".atomproc_ring_carbon_new.bin"
#define RING_CARBON_FIND_DUMP_PATH ".atomproc_ring_carbon_find.bin"
#define RING_CARBON_ORDER_DUMP_PATH ".atomproc_ring_carbon_order.bin"
#define RING_CARBON_CONCAT_DUMP_PATH ".atomproc_ring_carbon_concat.bin"
#define DISTANCE_DUMP_PATH ".atomproc_distance.bin"
#define FLT_MIN std::numeric_limits<float>::min()
#define FLT_MAX std::numeric_limits<float>::max()

using namespace AP;


std::unique_ptr<AtomProc> AtomProc::make_atom_proc(std::string data_path, std::string bin_path, bool dump_temp_products) {
    std::unique_ptr<AtomProc> proc(new AtomProc());
    proc->process_results = std::nullopt;
    proc->dump_temp_products = dump_temp_products;
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

        if (dump_temp_products) {
            std::ofstream writer(bin_path, std::ios::binary);
            mat_dump(writer, proc->atom_data_molecule);
            mat_dump(writer, proc->atom_bond);
            mat_dump(writer, proc->atom_angle);
            mat_dump(writer, proc->atom_impropers);
            mat_dump(writer, proc->atom_dihedrals);
            writer.close();
        }
    }

    std::cout << "Report: " << std::endl;
    std::cout << "#atom_data_molecule: " << proc->atom_data_molecule.rows() << std::endl;
    std::cout << "#atom_bond: " << proc->atom_bond.rows() << std::endl;
    std::cout << "#atom_angle: " << proc->atom_angle.rows() << std::endl;
    std::cout << "#atom_impropers: " << proc->atom_impropers.rows() << std::endl;
    std::cout << "#atom_dihedrals: " << proc->atom_dihedrals.rows() << std::endl;
    
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
                for (int i = 0; i < 6; i++) {
                    dump >> atom_dihedrals(j, i);
                }
                // double trash;
                // dump >> trash;
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

std::optional<std::string> AP::AtomProc::dump_result_to_data(std::string header_file, std::string data_file) {
    std::ifstream header(header_file);
    if (!header.good()) {
        THROW std::string("Bad header: ") + header_file;
    }
    std::stringstream buf;
    buf << header.rdbuf();
    std::string header_str = buf.str();
    header.close();

    std::ofstream writer(data_file);
    if (!writer.good()) {
        THROW std::string("Bad path to write: ") + data_file;
    }
    writer << "LAMMPS data file via write_data, version 8 Apr 2021, timestep = 0" << std::endl << std::endl;
    writer << process_results->atom_data_molecule.rows() << " atoms" << std::endl;
    writer << "9 atom types" << std::endl;
    writer << process_results->atom_bond.rows() << " bonds" << std::endl;
    writer << "11 bond types" << std::endl;
    writer << process_results->atom_angle.rows() << " angles" << std::endl;
    writer << "17 angle types" << std::endl;
    writer << process_results->atom_dihedrals.rows() << " dihedrals" << std::endl;
    writer << "22 dihedral types" << std::endl;
    writer << process_results->atom_impropers.rows() << " impropers" << std::endl;
    writer << "11 improper types" << std::endl << std::endl;

    writer << header_str << std::endl;
    writer << "Atoms #full" << std::endl << std::endl;

    double trash;
    for (int r = 0; r < process_results->atom_data_molecule.rows(); r++) {
        for (int c = 0; c < process_results->atom_data_molecule.cols(); c++) {
            double val = process_results->atom_data_molecule(r, c);
            if (modf(val, &trash) < FLT_MIN && val > 100000.0) {
                writer << ((long) val) << " ";
            } else {
                writer << val << " ";
            }
        }
        writer << std::endl;
    }
    writer << std::endl;

    writer << "Bonds" << std::endl << std::endl;
    for (int r = 0; r < process_results->atom_bond.rows(); r++) {
        for (int c = 0; c < process_results->atom_bond.cols(); c++) {
            double val = process_results->atom_bond(r, c);
            if (modf(val, &trash) < FLT_MIN && val > 100000.0) {
                writer << ((long) val) << " ";
            } else {
                writer << val << " ";
            }
        }
        writer << std::endl;
    }
    writer << std::endl;

    writer << "Angles" << std::endl << std::endl;
    for (int r = 0; r < process_results->atom_angle.rows(); r++) {
        for (int c = 0; c < process_results->atom_angle.cols(); c++) {
            double val = process_results->atom_angle(r, c);
            if (modf(val, &trash) < FLT_MIN && val > 100000.0) {
                writer << ((long) val) << " ";
            } else {
                writer << val << " ";
            }
        }
        writer << std::endl;
    }
    writer << std::endl;

    writer << "Dihedrals" << std::endl << std::endl;
    for (int r = 0; r < process_results->atom_dihedrals.rows(); r++) {
        for (int c = 0; c < process_results->atom_dihedrals.cols(); c++) {
            double val = process_results->atom_dihedrals(r, c);
            if (modf(val, &trash) < FLT_MIN && val > 100000.0) {
                writer << ((long) val) << " ";
            } else {
                writer << val << " ";
            }
        }
        writer << std::endl;
    }
    writer << std::endl;

    writer << "Impropers" << std::endl << std::endl;
    for (int r = 0; r < process_results->atom_impropers.rows(); r++) {
        for (int c = 0; c < process_results->atom_impropers.cols(); c++) {
            double val = process_results->atom_impropers(r, c);
            if (modf(val, &trash) < FLT_MIN && val > 100000.0) {
                writer << ((long) val) << " ";
            } else {
                writer << val << " ";
            }
        }
        writer << std::endl;
    }
    writer << std::endl;

    writer.close();

    return std::nullopt;
}

std::optional<std::string> AP::AtomProc::run(bool cache_result) {
    std::vector<int> position = locate(atom_data_molecule.col(2), (double) delete_type);
    process_results = std::nullopt;
    if (position.size() == 0) {
        THROW "Cannot find delete_type in atom_data_molecule";
    }
    Eigen::MatrixXd ring_carbon = atom_data_molecule(position, Eigen::all);
    // 将在一个苯环上的芳香碳原子放在连续六行
    Eigen::MatrixXd ring_carbon_new = Eigen::MatrixXd::Zero(length(ring_carbon), 13);
    if (auto ring_carbon_new_opt = mat_load_if_exists(RING_CARBON_NEW_DUMP_PATH)) {
        ring_carbon_new = *ring_carbon_new_opt;
    } else {
        int mol_number = (int) ring_carbon.col(1).maxCoeff();
        int max_threads = min((int) std::thread::hardware_concurrency(), mol_number);
        // max_threads = 1;
        int p_per_thread = mol_number / max_threads;
        std::cout << "Parallel info: " << max_threads << " threads for " << mol_number << " molecules" << std::endl;
        std::cout << "Parallel info: " << p_per_thread << " molecules per thread" << std::endl;
        std::vector<std::thread> pool;
        for (int i = 0; i < max_threads + 1; i++) {
            int start = p_per_thread * i;
            int count = p_per_thread;
            pool.emplace_back(&AtomProc::process_ring_carbon, this, std::ref(ring_carbon), std::ref(ring_carbon_new), mol_number, start, count);
        }
        for (int i = 0; i < pool.size(); i++) {
            pool[i].join();
        }
        if (dump_temp_products) {
            mat_dump_to_path(RING_CARBON_NEW_DUMP_PATH, ring_carbon_new);
        }
    }

    Eigen::MatrixXd ring_carbon_find = Eigen::MatrixXd::Zero(length(ring_carbon_new), 19);
    if (auto ring_carbon_find_opt = mat_load_if_exists(RING_CARBON_FIND_DUMP_PATH)) {
        ring_carbon_find = *ring_carbon_find_opt;
    } else {
        int max_threads = min(std::thread::hardware_concurrency(), (unsigned int) ring_carbon_find.rows());
        int p_per_thread = (int) ring_carbon_find.rows() / max_threads;
        std::vector<std::thread> pool;
        for (int i = 0; i < max_threads + 1; i++) {
            int start = p_per_thread * i;
            int count = p_per_thread;
            pool.emplace_back(&AtomProc::process_ring_carbon_find, this, std::ref(ring_carbon_find), std::ref(ring_carbon_new), start, count);
        }
        for (int i = 0; i < pool.size(); i++) {
            pool[i].join();
        }
        std::vector<int> non_zeros;
        for (int i = 0; i < ring_carbon_find.rows(); i++) {
            if (ring_carbon_find(i, 0) != 0.0) {
                non_zeros.push_back(i);
            }
        }
        Eigen::MatrixXd ring_carbon_find_filtered = ring_carbon_find(non_zeros, Eigen::all);
        ring_carbon_find = ring_carbon_find_filtered;
        if (dump_temp_products) {
            mat_dump_to_path(RING_CARBON_FIND_DUMP_PATH, ring_carbon_find);
        }
    }

    // mat_csv_dump("ring_carbon_find.csv", ring_carbon_find);

    Eigen::MatrixXd ring_carbon_order = ring_carbon;
    if (auto ring_carbon_order_opt = mat_load_if_exists(RING_CARBON_ORDER_DUMP_PATH)) {
        ring_carbon_order = *ring_carbon_order_opt;
    } else {
        Eigen::MatrixXd order = ring_carbon_find(Eigen::all, Eigen::seq(13, 18));
        Eigen::MatrixXd ring = Eigen::MatrixXd::Zero(length(order), 6);
        parallel_for(0, length(order) - 1, [&](int i) {
            Eigen::VectorXd order_i = order.row(i);
            std::sort(order_i.begin(), order_i.end());
            ring.row((Eigen::Index) i) = order_i;
        });
        // mat_csv_dump("ring.csv", ring);
        std::vector<int> indx;
        for (int i = 0; i < ring.rows(); i++) {
            indx.push_back(i);
        }
        std::sort(indx.begin(), indx.end(), [&](auto a, auto b) {
            return ring(a, 0) < ring(b, 0);
        });
        Eigen::MatrixXd ring_sorted = ring(indx, Eigen::all);
        Eigen::MatrixXd ring_sort = ring_sorted(Eigen::seq(0, length(ring_sorted) - 1, 4), Eigen::all).transpose();
        // mat_csv_dump("ring_sort.csv", ring_sort);
        parallel_for(0, length(ring_carbon) - 1, [&](int i) {
            ring_carbon_order.row(i) = ring_carbon(locate_amount(ring_carbon.col(0), ring_sort(i), 1), Eigen::all);
        });
        if (dump_temp_products) {
            mat_dump_to_path(RING_CARBON_ORDER_DUMP_PATH, ring_carbon_order);
        }
    }
    // mat_csv_dump("ring_carbon_order.csv", ring_carbon_order);
    ring_carbon = ring_carbon_order;

    Eigen::MatrixXd ring_carbon_concat(ring_carbon.rows(), ring_carbon.cols() + 1);
    if (auto ring_carbon_concat_opt = mat_load_if_exists(RING_CARBON_CONCAT_DUMP_PATH)) {
        ring_carbon_concat = *ring_carbon_concat_opt;
    } else {
        Eigen::MatrixXd N_bond = Eigen::MatrixXd::Zero(length(ring_carbon), 1);
        Eigen::MatrixXi atom_bond_sub = atom_bond(Eigen::all, Eigen::seq(2, 3)).cast<int>();
        parallel_for(0, length(ring_carbon) - 1, [&](int i) {
            N_bond(i, 0) = atom_bond_sub.cwiseEqual((int) ring_carbon(i, 0)).cast<int>().sum();
        });

        ring_carbon_concat.block(0, 0, ring_carbon.rows(), ring_carbon.cols()) = ring_carbon;
        ring_carbon_concat.col(ring_carbon_concat.cols() - 1) = N_bond;
        if (dump_temp_products) {
            mat_dump_to_path(RING_CARBON_CONCAT_DUMP_PATH, ring_carbon_concat);
        }
    }
    // mat_csv_dump("ring_carbon_concat.csv", ring_carbon_concat);
    ring_carbon = ring_carbon_concat;
    std::vector<int> indx;
    for (int i = 0; i < ring_carbon.rows(); i++) {
        indx.push_back(i);
    }
    std::stable_sort(indx.begin(), indx.end(), [&](int a, int b) {
        // if (ring_carbon(a, 1) == ring_carbon(b, 1)) {
        //     return ring_carbon(a, 0) < ring_carbon(b, 0);
        // }
        return ring_carbon(a, 1) < ring_carbon(b, 1);
    });
    Eigen::MatrixXd ring_carbon_sort_del = ring_carbon(indx, Eigen::all);
    Eigen::VectorXi to_select = 1 - (ring_carbon_sort_del.col(10).array() == 3.0).cast<int>();
    std::vector<int> ring_carbon_del_sel;
    for (int i = 0; i < to_select.size(); i++) {
        if (to_select(i) == 1) {
            ring_carbon_del_sel.push_back(i);
        }
    }
    Eigen::MatrixXd ring_carbon_sort = ring_carbon_sort_del(ring_carbon_del_sel, Eigen::all);

    // mat_csv_dump("ring_carbon_sort.csv", ring_carbon_sort);

    Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(length(ring_carbon_sort), 1);
    if (auto distance_opt = mat_load_if_exists(DISTANCE_DUMP_PATH)) {
        distance = *distance_opt;
    } else {
        parallel_for(0, length(ring_carbon_sort) / 4 - 1, [&](int i) {
            Eigen::VectorXd atom_1 = ring_carbon_sort(4 * i, Eigen::seq(4, 6));
            for (int j = 1; j < 4; j++) {
                Eigen::VectorXd deta_d = (Eigen::VectorXd) (ring_carbon_sort(4 * i + j, Eigen::seq(4, 6))) - atom_1;
                if (abs(deta_d(0)) > 7.0) {
                    deta_d(0) = abs(deta_d(0)) - a;
                }
                if (abs(deta_d(1)) > 7.0) {
                    deta_d(1) = abs(deta_d(1)) - b;
                }
                if (abs(deta_d(2)) > 7.0) {
                    deta_d(2) = abs(deta_d(2)) - c;
                }
                distance(4 * i + j, 0) = sqrt(deta_d.dot(deta_d));
            }
        });
    }

    // mat_csv_dump("distance.csv", distance);
    
    double max_distance = distance.maxCoeff();
    Eigen::MatrixXd ring_carbon_sort_dist(ring_carbon_sort.rows(), ring_carbon_sort.cols() + 1);
    ring_carbon_sort_dist.block(0, 0, ring_carbon_sort.rows(), ring_carbon_sort.cols()) = ring_carbon_sort;
    ring_carbon_sort_dist.block(0, ring_carbon_sort.cols(), ring_carbon_sort.rows(), 1) = distance;
    ring_carbon_sort = ring_carbon_sort_dist;

    // mat_csv_dump("ring_carbon_sort.csv", ring_carbon_sort);

    Eigen::MatrixXd delete_atomid = Eigen::MatrixXd::Zero(length(ring_carbon_sort) / 2, 1);
    parallel_for(0, length(ring_carbon_sort) / 4 - 1, [&](int i) {
        delete_atomid(2 * i, 0) = ring_carbon_sort(4 * i, 0);
        Eigen::Index p, q;
        distance(Eigen::seq(4 * i + 1, 4 * i + 3), Eigen::all).minCoeff(&p, &q);
        delete_atomid(2 * i + 1, 0) = ring_carbon_sort(4 * i + p + 1, 0);
    });

    // mat_csv_dump("delete_atomid.csv", delete_atomid);

    // delete rows & create new respective final matrices
    // atom_molecule deletion
    Eigen::MatrixXd atom_data_molecule_new(atom_data_molecule.rows(), atom_data_molecule.cols());
    std::vector<int> keep;
    // if current row is not in delete_atomid, then this row is to be kept
    parallel_for(0, (int) atom_data_molecule.rows() - 1, [&](int i) {
        if (locate_amount(delete_atomid, atom_data_molecule(i, 0), 1).size() == 0) {
            std::unique_lock<std::mutex> lock(mutex);
            keep.push_back(i);
        }
    });
    atom_data_molecule_new = atom_data_molecule(keep, Eigen::all);

    // mat_csv_dump("atom_data_molecule_new.csv", atom_data_molecule_new);

    // atom_bond deletion
    keep.clear();
    Eigen::MatrixXd atom_bond_new(atom_bond.rows(), atom_bond.cols());
    parallel_for(0, (int) atom_bond.rows() - 1, [&](int i) {
        if (locate_amount_multiple(delete_atomid, { atom_bond(i, 2), atom_bond(i, 3) }, 1).size() == 0) {
            std::unique_lock<std::mutex> lock(mutex);
            keep.push_back(i);
        }
    });
    atom_bond_new = atom_bond(keep, Eigen::all);

    // mat_csv_dump("atom_bond_new.csv", atom_bond_new);

    // atom_angle deletion
    keep.clear();
    Eigen::MatrixXd atom_angle_new(atom_angle.rows(), atom_angle.cols());
    parallel_for(0, (int) atom_angle.rows() - 1, [&](int i) {
        if (locate_amount_multiple(delete_atomid, { atom_angle(i, 2), atom_angle(i, 3), atom_angle(i, 4) }, 1).size() == 0) {
            std::unique_lock<std::mutex> lock(mutex);
            keep.push_back(i);
        }
    });
    atom_angle_new = atom_angle(keep, Eigen::all);

    // mat_csv_dump("atom_angle_new.csv", atom_angle_new);

    // atom_dihedrals deletion
    keep.clear();
    Eigen::MatrixXd atom_dihedrals_new(atom_dihedrals.rows(), atom_dihedrals.cols());
    parallel_for(0, (int) atom_dihedrals.rows() - 1, [&](int i) {
        if (locate_amount_multiple(delete_atomid, { atom_dihedrals(i, 2), atom_dihedrals(i, 3), atom_dihedrals(i, 4), atom_dihedrals(i, 5) }, 1).size() == 0) {
            std::unique_lock<std::mutex> lock(mutex);
            keep.push_back(i);
        }
    });
    atom_dihedrals_new = atom_dihedrals(keep, Eigen::all);

    // mat_csv_dump("atom_dihedrals_new.csv", atom_dihedrals_new);

    // atom_impropers deletion
    keep.clear();
    Eigen::MatrixXd atom_impropers_new(atom_impropers.rows(), atom_impropers.cols());
    parallel_for(0, (int) atom_impropers.rows() - 1, [&](int i) {
        if (locate_amount_multiple(delete_atomid, { atom_impropers(i, 2), atom_impropers(i, 3), atom_impropers(i, 4), atom_impropers(i, 5) }, 1).size() == 0) {
            std::unique_lock<std::mutex> lock(mutex);
            keep.push_back(i);
        }
    });
    atom_impropers_new = atom_impropers(keep, Eigen::all);

    // mat_csv_dump("atom_impropers_new.csv", atom_impropers_new);

    process_results = {
        atom_data_molecule_new,
        atom_bond_new,
        atom_angle_new,
        atom_dihedrals_new,
        atom_impropers_new
    };

    return std::nullopt;
}

AtomProc::AtomProc() {
    atom_masses.col(0) = Eigen::VectorXd::LinSpaced(Natoms, 0, (double) Natoms - 1);
    atom_masses.col(1).array() = 12.0;
}

void AP::AtomProc::process_ring_carbon(Eigen::MatrixXd& ring_carbon, Eigen::MatrixXd& ring_carbon_new, int mol_number, int offset, int count) {
    for (int i = offset; i < min(offset + count, mol_number); i++) {
        std::vector<int> mol = locate(ring_carbon.col(1), (double) (i + 1));
        Eigen::MatrixXd ring_carbon_mol = ring_carbon(mol, Eigen::all);

        for (int j = 0; j < length(ring_carbon_mol); j++) {
            Eigen::MatrixXd atom_id = Eigen::MatrixXd::Zero(3, 1);
            std::vector<int> atom_1 = locate_amount(atom_bond.col(2), ring_carbon_mol(j, 0), 2);

            int atom_2_loc = atom_1.size() <= 1 ? 2 : 1;
            std::vector<int> atom_2 = locate_amount(atom_bond.col(3), ring_carbon_mol(j, 0), atom_2_loc);
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

            // std::cout << ring_carbon_mol.row(j) << std::endl;
            // std::cout << atom_id.transpose() << std::endl;
            ring_carbon_new.row(i * length(ring_carbon_mol) + j) << ring_carbon_mol.row(j), atom_id.transpose();
        }
    }
}

void AP::AtomProc::process_ring_carbon_find(Eigen::MatrixXd& ring_carbon_find, const Eigen::MatrixXd& ring_carbon_new, int offset, int count) {
    for (int k = offset; k < min(offset + count, length(ring_carbon_new)); k++) {
        if (ring_carbon_new(k, 12) != 0) {
            continue;
        }
        ring_carbon_find.block(k, 0, 1, 13) = ring_carbon_new.block(k, 0, 1, 13);
        ring_carbon_find(k, 13) = ring_carbon_new(k, 0);
        std::vector<int> atom_1_row = locate_amount(ring_carbon_new.col(0), ring_carbon_new(k, 10), 1);
        std::vector<int> atom_2_row = locate_amount(ring_carbon_new.col(0), ring_carbon_new(k, 11), 1);

        if (atom_1_row.size() > 0 && ring_carbon_new(atom_1_row[0], 12) > 0) {
            ring_carbon_find(k, 14) = ring_carbon_new(atom_1_row[0], 0);
        } 
        if (atom_2_row.size() > 0 && ring_carbon_new(atom_2_row[0], 12) > 0) {
            ring_carbon_find(k, 14) = ring_carbon_new(atom_2_row[0], 0);
        }

        // Eigen::MatrixXd sel_atom_3 = ring_carbon_new(atom_3_row, Eigen::seq(10, 12));
        int atom_3_row = locate_amount(ring_carbon_new.col(0), ring_carbon_find(k, 14), 1)[0];
        std::vector<double> sel_atom_3;
        for (int i = 10; i < 13; i++) {
            if (ring_carbon_new(atom_3_row, i) == ring_carbon_new(k, 0)) {
                continue;
            }
            sel_atom_3.push_back(ring_carbon_new(atom_3_row, i));
        }
        int sel_atom_indx = 0;
        if (locate_amount(ring_carbon_new.col(0), sel_atom_3[0], 1).size() == 0) {
            sel_atom_indx = 1;
        }
        ring_carbon_find(k, 15) = sel_atom_3[sel_atom_indx];

        int atom_4_row = locate_amount(ring_carbon_new.col(0), ring_carbon_find(k, 15), 1)[0];
        double sel_atom_4 = ring_carbon_new(atom_4_row, 10);
        if (sel_atom_4 == ring_carbon_find(k, 14)) {
            sel_atom_4 = ring_carbon_new(atom_4_row, 11);
        }
        ring_carbon_find(k, 16) = sel_atom_4;

        int atom_5_row = locate_amount(ring_carbon_new.col(0), ring_carbon_find(k, 16), 1)[0];
        double sel_atom_5 = ring_carbon_new(atom_5_row, 10);
        if (sel_atom_5 == ring_carbon_find(k, 15)) {
            sel_atom_5 = ring_carbon_new(atom_5_row, 11);
        }
        ring_carbon_find(k, 17) = sel_atom_5;

        int atom_6_row = locate_amount(ring_carbon_new.col(0), ring_carbon_find(k, 17), 1)[0];
        std::vector<double> sel_atom_6;
        for (int i = 10; i < 13; i++) {
            if (ring_carbon_new(atom_6_row, i) == ring_carbon_find(k, 16)) {
                continue;
            }
            sel_atom_6.push_back(ring_carbon_new(atom_6_row, i));
        }
        sel_atom_indx = 0;
        if (locate_amount(ring_carbon_new.col(0), sel_atom_6[0], 1).size() == 0) {
            sel_atom_indx = 1;
        }
        ring_carbon_find(k, 18) = sel_atom_6[sel_atom_indx];
    }
}


void AP::mat_dump(std::ofstream& out, const Eigen::MatrixXd& mat) {
    int rows = (int) mat.rows(), cols = (int) mat.cols();
    out.put('a');
    out.write((const char *) &rows, sizeof(int));
    out.write((const char *) &cols, sizeof(int));
    out.write((const char *) &mat(0, 0), (long long) rows * cols * sizeof(double));
}

void AP::mat_csv_dump(std::string path, const Eigen::MatrixXd& mat) {
    std::ofstream writer(path);
    if (!writer.good()) {
        std::cout << "ERR! Cannot dump .csv to " << path << std::endl;
        return;
    }
    // ;-separated CSV
    for (int r = 0; r < mat.rows(); r++) {
        for (int c = 0; c < mat.cols(); c++) {
            writer << std::fixed << mat(r, c) << ";";
        }
        writer << std::endl;
    }
    writer.close();
}

std::optional<Eigen::MatrixXd> AP::mat_load_if_exists(std::string path) {
    if (!std::filesystem::exists(path)) {
        return std::nullopt;
    }
    std::ifstream reader(path, std::ios::binary);
    Eigen::MatrixXd mat = mat_load(reader);
    reader.close();
    return mat;
}

void AP::mat_dump_to_path(std::string path, const Eigen::MatrixXd &mat) {
    std::ofstream writer(path, std::ios::binary);
    if (!writer.good()) {
        std::cerr << "ERR! Bad path to dump: " << path << std::endl;
        return;
    }
    mat_dump(writer, mat);
    writer.close();
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

void AP::parallel_for(int begin, int end, std::function<void(int)> func, std::optional<int> expected_threads) {
    std::vector<std::thread> pool;
    int num_loops = end - begin + 1;
    int max_threads = -1;
    if (expected_threads == std::nullopt) {
        max_threads = min((int) std::thread::hardware_concurrency(), num_loops);
    } else {
        max_threads = *expected_threads;   
    }
    int p_per_thread = num_loops / max_threads;

    for (int i = 0; i < max_threads + 1; i++) {
        int b = begin + (i * p_per_thread);
        int e = min(end + 1, b + p_per_thread);
        pool.emplace_back([](int b, int e, std::function<void(int)> func) {
            for (int i = b; i < e; i++) {
                func(i);
            }
        }, b, e, func);
    }
    for (int i = 0; i < pool.size(); i++) {
        pool[i].join();
    }
}

int AP::length(Eigen::MatrixXd mat) {
    return (int) (mat.rows() > mat.cols() ? mat.rows() : mat.cols());
}

std::vector<int> AP::locate(Eigen::MatrixXd vec, double value) {
    // TODO: find a faster method
    // auto begin = std::chrono::high_resolution_clock::now();
    std::vector<int> ret;
    for (int i = 0; i < vec.size(); i++) {
        if (vec(i) == value) {
            ret.push_back(i);
        }
    }
    // auto end = std::chrono::high_resolution_clock::now();
    // std::cout << "Profile: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;

    return ret;
}

std::vector<int> AP::locate_amount(Eigen::MatrixXd vec, double value, int n) {
    std::vector<int> ret;
    for (int i = 0; n != 0 && i < vec.size(); i++) {
        if (vec(i) == value) {
            ret.push_back(i);
            n--;
        }
        
    }
    return ret;
}

std::vector<int> AP::locate_amount_multiple(Eigen::MatrixXd vec, std::vector<double> value, int n) {
    std::vector<int> ret;
    for (int i = 0; n != 0 && i < vec.size(); i++) {
        for (double val : value) {
            if (vec(i) == val) {
                ret.push_back(i);
                n--;
                break;
            }
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
