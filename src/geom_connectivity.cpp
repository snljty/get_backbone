#include "get_backbone.hpp"

Backbone_extracter::Backbone_extracter(std::string_view ifilename) : mol_origin(ifilename) {}

void Backbone_extracter::get_connections() {
    Eigen::MatrixXd atomic_covalence_radius_sum = 
        mol_origin.atomic_covalence_radius.replicate(1, mol_origin.natoms) + 
        mol_origin.atomic_covalence_radius.transpose().replicate(mol_origin.natoms, 1);

    // column-major
    Eigen::MatrixXd coords_expanded_row = mol_origin.coordinates.replicate(mol_origin.natoms, 1);
    coords_expanded_row.resize(ncoords, mol_origin.natoms * mol_origin.natoms);
    Eigen::MatrixXd coords_expanded_col = mol_origin.coordinates.replicate(1, mol_origin.natoms);
    
    Eigen::MatrixXd diff = coords_expanded_col - coords_expanded_row;
    Eigen::MatrixXd distance = diff.colwise().norm().reshaped(mol_origin.natoms, mol_origin.natoms);

    constexpr double scaler_tol = 1.15;

    connections = (distance.array() <= (atomic_covalence_radius_sum * scaler_tol).array()).cast<int>();
    connections.diagonal().setZero();
}

void Backbone_extracter::set_connect(int iatom, int jatom) {
    if (iatom == jatom) throw std::invalid_argument("Error: cannot set an atom connect to itself.");
    connections(iatom - 1, jatom - 1) = connections(jatom - 1, iatom - 1) = 1;
}

void Backbone_extracter::set_disconnect(int iatom, int jatom) {
    if (iatom == jatom) return;
    connections(iatom - 1, jatom - 1) = connections(jatom - 1, iatom - 1) = 0;
}

