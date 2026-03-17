#include "get_backbone.hpp"

const PeriodicTable& PeriodicTable::get() {
    static PeriodicTable instance;
    return instance;
}

int PeriodicTable::get_atomic_number(std::string_view symbol) const {
    return elements_map.at(symbol);
}

PeriodicTable::PeriodicTable() {
    elements_map.reserve(max_element_num + 1);
    for (int atomic_number = 0; atomic_number <= max_element_num; ++ atomic_number) {
        elements_map.insert({element_names[atomic_number], atomic_number});
    }
}

