#include "get_backbone.hpp"

void Backbone_extracter::get_backbone() {
    get_connections();
    Eigen::VectorXi num_connections = connections.rowwise().sum();
    // get alkyl
    alkyl.reserve(mol_origin.natoms);
    std::vector<bool> in_alkyl(mol_origin.natoms, false);
    for (int iatom = 0; iatom < mol_origin.natoms; ++ iatom) {
        if (mol_origin.elements[iatom] == "C" && num_connections(iatom) == 4) {
            alkyl.push_back(iatom);
            in_alkyl[iatom] = true;
            for (int jatom = 0; jatom < mol_origin.natoms; ++ jatom) {
                if (num_connections[jatom] == 1 && connections(iatom, jatom)) {
                    alkyl.push_back(jatom);
                    in_alkyl[jatom] = true;
                }
            }
        }
    }
    std::sort(alkyl.begin(), alkyl.end());
    // get backbone
    backbone.clear();
    for (int iatom = 0; iatom < mol_origin.natoms; ++ iatom) {
        if (!in_alkyl[iatom]) backbone.push_back(iatom);
    }

    // get connection site of alkyl, and atoms in alkyl that connects to them.
    alkyl_connection_site.clear();
    std::vector<std::vector<int> > alkyl_connection_site_connected_alkyl;
    int natoms_truncated_methyl_total = 0;
    for (int iatom = 0; iatom < mol_origin.natoms; ++ iatom) {
        if (in_alkyl[iatom]) {
            bool being_alkyl_connection_site = false;
            for (int jatom = 0; jatom < mol_origin.natoms; ++ jatom) {
                if (!in_alkyl[jatom] && connections(iatom, jatom)) {
                    being_alkyl_connection_site = true;
                    break;
                }
            }
            if (being_alkyl_connection_site) {
                alkyl_connection_site.push_back(iatom);
                alkyl_connection_site_connected_alkyl.emplace_back();
                for (int jatom = 0; jatom < mol_origin.natoms; ++ jatom) {
                    if (in_alkyl[jatom] && connections(iatom, jatom)) {
                        alkyl_connection_site_connected_alkyl.back().push_back(jatom);
                    }
                }
                natoms_truncated_methyl_total += 1 + alkyl_connection_site_connected_alkyl.back().size();
            }
        }
    }

    constexpr double carbon_hydrogen_bond_length = 1.09105;
    int natoms_trimmed = backbone.size() + natoms_truncated_methyl_total;
    hydrogens_to_optimize.clear();
    mol_trimmed.set_natoms(natoms_trimmed);
    mol_trimmed.filename_prefix = mol_origin.filename_prefix + "_trimmed";
    for (size_t i = 0; i < backbone.size(); ++ i) {
        mol_trimmed.elements[i] = mol_origin.elements[backbone[i]];
        mol_trimmed.coordinates.col(i) = mol_origin.coordinates.col(backbone[i]);
    }
    int iatom_in_methyl = 0;
    for (size_t imethyl = 0; imethyl < alkyl_connection_site.size(); ++ imethyl) {
        int iatom = alkyl_connection_site[imethyl];
        mol_trimmed.elements[backbone.size() + iatom_in_methyl] = mol_origin.elements[iatom];
        mol_trimmed.coordinates.col(backbone.size() + iatom_in_methyl) = mol_origin.coordinates.col(iatom);
        ++ iatom_in_methyl;
        for (size_t i = 0; i < alkyl_connection_site_connected_alkyl[imethyl].size(); ++ i) {
            int jatom = alkyl_connection_site_connected_alkyl[imethyl][i];
            mol_trimmed.elements[backbone.size() + iatom_in_methyl + i] = "H";
            Eigen::Vector3d current_bond = mol_origin.coordinates.col(jatom) - mol_origin.coordinates.col(iatom);
            mol_trimmed.coordinates.col(backbone.size() + iatom_in_methyl + i) = 
                mol_origin.coordinates.col(iatom) + current_bond / current_bond.norm() * carbon_hydrogen_bond_length;
            hydrogens_to_optimize.push_back(backbone.size() + iatom_in_methyl + i);
        }
        iatom_in_methyl += alkyl_connection_site_connected_alkyl[imethyl].size();
    }

    backbone_no_hydrogen.clear();
    backbone_no_hydrogen.reserve(backbone.size());
    for (int iatom : backbone) {
        if (mol_origin.elements[iatom] != "H") backbone_no_hydrogen.push_back(iatom);
    }

    mol_backbone_no_hydrogen.set_natoms(backbone_no_hydrogen.size());
    mol_backbone_no_hydrogen.filename_prefix = mol_origin.filename_prefix + "_backbone_only";
    for (size_t i = 0; i < backbone_no_hydrogen.size(); ++ i) {
        mol_backbone_no_hydrogen.elements[i] = mol_origin.elements[backbone_no_hydrogen[i]];
        mol_backbone_no_hydrogen.coordinates.col(i) = mol_origin.coordinates.col(backbone_no_hydrogen[i]);
    }
}

void Backbone_extracter::write_results() const {
    mol_trimmed.write_gjf();
    mol_backbone_no_hydrogen.write_gjf();
    fmt::print("{:s}\n", list_to_indices_str_from_1(backbone));
    fmt::print("{:s}\n", list_to_indices_str_from_1(hydrogens_to_optimize));
}

