#include "get_backbone.hpp"

void Geom_file::write(std::string_view keyword) const {
    if (filename_suffix == ".gjf") {
        write_gjf(keyword);
    } else if (filename_suffix == ".xyz") {
        write_xyz();
    } else {
        throw std::invalid_argument("Error: the suffix of the given file is not recognizable.");
    }
}

void Geom_file::write_gjf(std::string_view keyword) const {
    if (filename_suffix != ".gjf") throw std::invalid_argument("Error: the suffix should be \".gjf\".");
    fmt::ostream ofile(fmt::output_file(filename_prefix + filename_suffix));

    ofile.print("%chk={0:s}.chk\n{1:s}\n\n{0:s}\n\n{2:d} {3:d}\n", filename_prefix, keyword, 0, 1);
    for (int i = 0; i < elements.size(); ++ i) {
        ofile.print(" {:<2s}    {:13.8f}    {:13.8f}    {:13.8f}\n", elements[i], 
            coordinates(coord_x, i), coordinates(coord_y, i), coordinates(coord_z, i));
    }
    ofile.print("\n");
    ofile.close();
}

void Geom_file::write_xyz() const {
    if (filename_suffix != ".xyz") throw std::invalid_argument("Error: the suffix should be \".xyz\".");
    fmt::ostream ofile(fmt::output_file(filename_prefix + filename_suffix));

    ofile.print("{:5d}\n", natoms);
    ofile.print("{:s}\n", filename_prefix);
    for (int i = 0; i < elements.size(); ++ i) {
        ofile.print(" {:<2s}    {:13.8f}    {:13.8f}    {:13.8f}\n", elements[i], 
            coordinates(coord_x, i), coordinates(coord_y, i), coordinates(coord_z, i));
    }
    ofile.close();
}

