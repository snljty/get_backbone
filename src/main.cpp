#include "get_backbone.hpp"

int main(int argc, const char *argv[]) {
    if (argc - 1 != 1) throw std::invalid_argument(fmt::format("Usage: {:s} xxx.gjf/.xyz", argv[0]));

    Backbone_extracter mol(argv[1]);
    mol.get_backbone();
    mol.write_results();

    return 0;
}

