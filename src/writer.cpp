#include "get_backbone.hpp"

void Gjf_file::write_gjf( std::string_view keyword) const {
    fmt::ostream ofile(fmt::output_file(filename_prefix + ".gjf"));

    ofile.print("%chk={0:s}.chk\n{1:s}\n\n{0:s}\n\n{2:d} {3:d}\n", filename_prefix, keyword, 0, 1);
    for (int i = 0; i < elements.size(); ++ i) {
        ofile.print(" {:<2s}    {:13.8f}    {:13.8f}    {:13.8f}\n", elements[i], 
            coordinates(coord_x, i), coordinates(coord_y, i), coordinates(coord_z, i));
    }
    ofile.print("\n");
    ofile.close();
}
