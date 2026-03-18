#pragma once
#ifndef __GET_BACKBONE_HPP__
#define __GET_BACKBONE_HPP__

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <string>
#include <unordered_map>
#include <regex>
#include <string_view>
#include <algorithm>
#include <stdexcept>

#include <Eigen/Core>

#include <fmt/format.h>
#include <fmt/os.h>

using coord_direction = enum {coord_x, coord_y, coord_z};

constexpr int ncoords = 3;

constexpr size_t max_element_num = 118;

constexpr std::array<std::string_view, max_element_num + 1> element_names = {
    "",
    "H" , "He", 
    "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne",
    "Na", "Mg", "Al", "Si", "P" , "S" , "Cl", "Ar",
    "K" , "Ca", "Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y" , "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I" , "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", 
    "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U" , "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", 
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nb", "Fl", "Mc", "Lv", "Ts", "Og"
};

constexpr std::array<float, max_element_num + 1> element_covalence_radius = {
    0.1, 
    0.31, 0.28, 1.28, 0.96, 0.84, 0.76,
    0.71, 0.66, 0.57, 0.58, 1.66, 1.41,
    1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
    2.03, 1.76, 1.70, 1.60, 1.53, 1.39,
    1.39, 1.32, 1.26, 1.24, 1.32, 1.22,
    1.22, 1.20, 1.19, 1.20, 1.20, 1.16,
    2.20, 1.95, 1.90, 1.75, 1.64, 1.54,
    1.47, 1.46, 1.42, 1.39, 1.45, 1.44,
    1.42, 1.39, 1.39, 1.38, 1.39, 1.40,
    2.44, 2.15, 2.07, 2.04, 2.03, 2.01,
    1.99, 1.98, 1.98, 1.96, 1.94, 1.92,
    1.92, 1.89, 1.90, 1.87, 1.87, 1.75,
    1.70, 1.62, 1.51, 1.44, 1.41, 1.36,
    1.36, 1.32, 1.45, 1.46, 1.48, 1.40,
    1.50, 1.50, 2.60, 2.21, 2.15, 2.06,
    2.00, 1.96, 1.90, 1.87, 1.80, 1.69,
    1.68, 1.68, 1.65, 1.67, 1.73, 1.76,
    1.61, 1.57, 1.49, 1.43, 1.41, 1.34,
    1.29, 1.28, 1.21, 1.22, 1.5 , 1.5 ,
    1.5 , 1.5 , 1.5 , 1.5
};

class Backbone_extracter;

class PeriodicTable {
public:
    static const PeriodicTable& get();
    int get_atomic_number(std::string_view symbol) const;

private:
    std::unordered_map<std::string_view, int> elements_map;
    PeriodicTable();
};

class Backbone_extracter;

class Geom_file {
private:
    int natoms;
    const PeriodicTable& periodic_table = PeriodicTable::get();
    std::vector<std::string> elements;
    Eigen::MatrixXd coordinates;
    Eigen::VectorXd atomic_covalence_radius;
    std::string filename_prefix, filename_suffix;
public:
    Geom_file() = default;
    Geom_file(std::string_view ifilename);
private:
    void set_natoms(int natoms);
    void read_gjf();
    void read_xyz();
    void write_gjf(std::string_view keywords="#P") const;
    void write_xyz() const;
    friend class Backbone_extracter;
public:
    void read(std::string_view ifilename);
    void write(std::string_view keywords="#P B3LYP/6-31G** EmpiricalDispersion=GD3BJ") const;
};

class Backbone_extracter {
private:
    Geom_file mol_origin;
    Geom_file mol_trimmed;
    Geom_file mol_backbone_no_hydrogen;
    Eigen::MatrixXi connections;
    std::vector<int> alkyl_connection_site, backbone, backbone_no_hydrogen, alkyl, hydrogens_to_optimize;

public:
    Backbone_extracter(std::string_view ifilename);
    void set_connect(int iatom, int jatom); // counts from 1
    void set_disconnect(int iatom, int jatom); // counts from 1
    void get_backbone();
    void write_results() const;
private:
    void get_connections();
};

std::vector<int> indices_str_to_list_from_0(const std::string& indices);
std::string list_to_indices_str_from_1(const std::vector<int>& nums);

#endif // __GET_BACKBONE_HPP__
