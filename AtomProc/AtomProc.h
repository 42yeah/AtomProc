#pragma once

#include <iostream>
#include <fstream>
#include <random>
#include <optional>
#include <filesystem>
#include <map>
#include <mutex>
#include <chrono>
#include <functional>
#include <Eigen/Core>


namespace AP {

    class AtomProc {
    public:
        AtomProc(AtomProc &) = delete;

        static std::unique_ptr<AtomProc> make_atom_proc(std::string data_path, std::string bin_path, bool dump_temp_products = false);

        bool parse_mat_from_data(std::string data_file);

        void random_sample();

        std::optional<std::string> dump_result_to_data(std::string header_file, std::string data_file);

#ifdef _DEBUG
        std::optional<std::string> run(bool cache_result = true);
#else
        std::optional<std::string> run(bool cache_result = false);
#endif

        struct Result {
            Eigen::MatrixXd atom_data_molecule;
            Eigen::MatrixXd atom_bond;
            Eigen::MatrixXd atom_angle;
            Eigen::MatrixXd atom_dihedrals;
            Eigen::MatrixXd atom_impropers;
        };

        std::optional<Result> process_results;

    private:
        AtomProc();

        // process ring carbon in [offset, offset + count).
        void process_ring_carbon(Eigen::MatrixXd &ring_carbon, Eigen::MatrixXd &ring_carbon_new, int mol_number, int offset, int count);

        // process ring_carbon_find in [offset, offset + count).
        void process_ring_carbon_find(Eigen::MatrixXd &ring_carbon_find, const Eigen::MatrixXd &ring_carbon_new, int offset, int count);

        // default parameters
        long Natoms = 108220;
        long Bond_number = 115940;
        long Angle_number = 139120;
        long Dihedrals_number = 162300;
        long Impropers_number = 15460;
        int delete_type = 3;
        double a = 2590.601718131646 - 2330.316618086725;
        double b = 1822.6294254308316 - 1562.3443253858327;
        double c = 478.34097815406955 - 260.0027269264819;

        Eigen::MatrixXd atom_data_molecule = Eigen::MatrixXd::Zero(Natoms, 10);
        Eigen::MatrixXd atom_bond = Eigen::MatrixXd::Zero(Bond_number, 4);
        Eigen::MatrixXd atom_angle = Eigen::MatrixXd::Zero(Angle_number, 5);
        Eigen::MatrixXd atom_dihedrals = Eigen::MatrixXd::Zero(Dihedrals_number, 6);
        Eigen::MatrixXd atom_impropers = Eigen::MatrixXd::Zero(Impropers_number, 6);
        Eigen::MatrixXd atom_masses = Eigen::MatrixXd::Zero(Natoms, 2);

        bool dump_temp_products;

        std::mutex mutex;
    };

    // length returns max(size(mat)) for nonempty mat
    int length(Eigen::MatrixXd mat);

    std::vector<int> locate(Eigen::MatrixXd vec, double value);

    std::vector<int> locate_amount(Eigen::MatrixXd vec, double value, int n);

    std::vector<int> locate_amount_multiple(Eigen::MatrixXd vec, std::vector<double> value, int n);

    bool starts_with(std::string what, std::string starts);

    void mat_dump(std::ofstream& out, const Eigen::MatrixXd& mat);

    void mat_csv_dump(std::string path, const Eigen::MatrixXd &mat);

    std::optional<Eigen::MatrixXd> mat_load_if_exists(std::string path);

    void mat_dump_to_path(std::string path, const Eigen::MatrixXd &mat);

    Eigen::MatrixXd mat_load(std::ifstream& in);

    void random_sample(Eigen::MatrixXd mat, int n_samples);

    void parallel_for(int begin, int end, std::function<void(int)> func, std::optional<int> expected_threads = std::nullopt);

    template<typename T>
    T min(T a, T b) {
        if (a < b) {
            return a;
        }
        return b;
    }

    template<typename T>
    T max(T a, T b) {
        if (a > b) {
            return a;
        }
        return b;
    }

}

