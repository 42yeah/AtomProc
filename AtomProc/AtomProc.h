#pragma once

#include <iostream>
#include <fstream>
#include <random>
#include <optional>
#include <filesystem>
#include <Eigen/Core>


namespace AP {

    class AtomProc {
    public:
        AtomProc(AtomProc &) = delete;

        static std::unique_ptr<AtomProc> make_atom_proc(std::string data_path, std::string bin_path);

        bool parse_mat_from_data(std::string data_file);

        void random_sample();

        std::optional<std::string> run();

    private:
        AtomProc();

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
    };

    // length returns max(size(mat)) for nonempty mat
    int length(Eigen::MatrixXd mat);

    std::vector<int> locate(Eigen::MatrixXd vec, double value);

    bool starts_with(std::string what, std::string starts);

    void mat_dump(std::ofstream& out, Eigen::MatrixXd& mat);

    Eigen::MatrixXd mat_load(std::ifstream& in);

    void random_sample(Eigen::MatrixXd mat, int n_samples);

}

