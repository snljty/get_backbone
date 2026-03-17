#include "get_backbone.hpp"

Gjf_file::Gjf_file(std::string_view ifilename) {
    read_gjf(ifilename);
}

void Gjf_file::set_natoms(int natoms) {
    this->natoms = natoms;
    elements.resize(natoms);
    coordinates.resize(ncoords, natoms);
    atomic_covalence_radius.resize(natoms);
}

void Gjf_file::read_gjf(std::string_view ifilename) {
    std::string::size_type pos = ifilename.rfind('.');
    if (pos == std::string::npos || ifilename.substr(pos) != ".gjf") {
        throw std::invalid_argument("Error: the suffix of the given gjf file is incorrect.");
    }
    filename_prefix = ifilename.substr(0, pos);

    std::ifstream ifile(std::string(ifilename), std::ios_base::binary);
    if (!ifile) throw std::ios_base::failure(fmt::format("Error: cannot open \"{:s}\" for reading.", ifilename));

    std::regex blank_line_pattern("\\s*");
    std::string line;

    while (std::getline(ifile, line) && !std::regex_match(line, blank_line_pattern)) {;}
    while (std::getline(ifile, line) && !std::regex_match(line, blank_line_pattern)) {;}

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

