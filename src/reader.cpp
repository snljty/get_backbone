#include "get_backbone.hpp"

Geom_file::Geom_file(std::string_view ifilename) {
    read(ifilename);
}

void Geom_file::set_natoms(int natoms) {
    this->natoms = natoms;
    elements.resize(natoms);
    coordinates.resize(ncoords, natoms);
    atomic_covalence_radius.resize(natoms);
}

void Geom_file::read(std::string_view ifilename) {
    std::string::size_type pos = ifilename.rfind('.');
    if (pos == std::string::npos || (ifilename.substr(pos) != ".gjf" && ifilename.substr(pos) != ".xyz")) {
        throw std::invalid_argument("Error: the suffix of the given file is not recognizable.");
    }

    filename_prefix = ifilename.substr(0, pos);
    filename_suffix = ifilename.substr(pos);
    if (filename_suffix == ".gjf") {
        read_gjf();
    } else if (filename_suffix == ".xyz") {
        read_xyz();
    } else {
        throw std::invalid_argument("Error: the suffix of the given file is not recognizable.");
    }
}

void Geom_file::read_gjf() {
    if (filename_suffix != ".gjf") throw std::invalid_argument("Error: the suffix should be \".gjf\".");
    std::ifstream ifile(filename_prefix + filename_suffix, std::ios_base::binary);
    if (!ifile) throw std::ios_base::failure(fmt::format("Error: cannot open \"{:s}{:s}\" for reading.", filename_prefix, filename_prefix));

    std::regex blank_line_pattern("\\s*");
    std::string line;

    bool found_blank_line = false;
    while (std::getline(ifile, line)) {
        if (std::regex_match(line, blank_line_pattern)) {
            found_blank_line = true;
            break;
        }
    }
    if (!found_blank_line) throw std::ios_base::failure("Error: cannot read title.");

    found_blank_line = false;
    while (std::getline(ifile, line)) {
        if (std::regex_match(line, blank_line_pattern)) {
            found_blank_line = true;
            break;
        }
    }
    if (!found_blank_line) throw std::ios_base::failure("Error: cannot read charge and multiplicity.");

    if (!std::getline(ifile, line)) throw std::ios_base::failure("Error: cannot read charge and multiplicity.");
    std::istringstream iss0(line);
    int charge, multiplicity;

    if (!(iss0 >> charge >> multiplicity)) throw std::ios_base::failure("Error: cannot read charge and multiplicity.");

    std::streampos file_pos = ifile.tellg();

    int natoms_now = 0;
    while (std::getline(ifile, line)) {
        if (std::regex_match(line, blank_line_pattern)) break;
        ++ natoms_now;
    }
    set_natoms(natoms_now);
    ifile.seekg(file_pos);

    for (int i = 0; i < natoms; ++ i) {
        std::getline(ifile, line);
        std::istringstream iss(line);
        if (!(iss >> elements[i] >> coordinates(coord_x, i) >> coordinates(coord_y, i) >> coordinates(coord_z, i))) {
            throw std::ios_base::failure(fmt::format("Error: cannot read atom {:d}.", i + 1));
        }
        atomic_covalence_radius[i] = element_covalence_radius[periodic_table.get_atomic_number(elements[i])];
    }

    ifile.close();
}

void Geom_file::read_xyz() {
    if (filename_suffix != ".xyz") throw std::invalid_argument("Error: the suffix should be \".xyz\".");
    std::ifstream ifile(filename_prefix + filename_suffix);
    if (!ifile) throw std::ios_base::failure(fmt::format("Error: cannot open \"{:s}{:s}\" for reading.", filename_prefix, filename_prefix));

    std::string line;

    if (!std::getline(ifile, line)) throw std::ios_base::failure("Error: cannot read amount of atoms.");
    std::istringstream iss0(line);
    int natoms_now;
    if (!(iss0 >> natoms_now)) throw std::ios_base::failure("Error: cannot read amount of atoms.");
    set_natoms(natoms_now);

    if (!std::getline(ifile, line)) throw std::ios_base::failure("Error: cannot get comment line.");

    for (int i = 0; i < natoms; ++ i) {
        if (!std::getline(ifile, line)) {
            throw std::ios_base::failure(fmt::format("Error: cannot read atom {:d}.", i + 1));
        }
        std::istringstream iss(line);
        if (!(iss >> elements[i] >> coordinates(coord_x, i) >> coordinates(coord_y, i) >> coordinates(coord_z, i))) {
            throw std::ios_base::failure(fmt::format("Error: cannot read atom {:d}.", i + 1));
        }
        atomic_covalence_radius[i] = element_covalence_radius[periodic_table.get_atomic_number(elements[i])];
    }

    ifile.close();
}

