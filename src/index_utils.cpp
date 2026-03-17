#include "get_backbone.hpp"

std::vector<int> indices_str_to_list_from_0(const std::string& indices) {
    std::vector<int> ret;
    std::istringstream indices_iss(indices);
    std::string current_str;
    size_t pos;
    int a, b;
    while (std::getline(indices_iss, current_str, ',')) {
        pos = current_str.find('-');
        if (pos != std::string::npos) {
            current_str[pos] = ' ';
            std::istringstream current_iss(current_str);
            if (!(current_iss >> a >> b)) throw std::invalid_argument("Error: cannot parse string.");
            for (int i = a; i <= b; ++ i) {
                ret.push_back(i - 1);
            }
        } else {
            std::istringstream current_iss(current_str);
            if (!(current_iss >> a)) throw std::invalid_argument("Error: cannot parse string.");
            ret.push_back(a - 1);
        }
    }
    return ret;
}

std::string list_to_indices_str_from_1(const std::vector<int>& nums) {
    if (nums.empty()) return "";
    std::vector<std::string> ret_list;
    int a, b;
    // handle first
    a = nums[0];
    b = a;
    // handle second to the one before last
    for (size_t i = 1; i < nums.size(); ++ i) {
        if (nums[i] == b + 1) {
            ++ b;
        } else {
            if (a == b) {
                ret_list.emplace_back(fmt::format("{:d}", a + 1));
            } else {
                ret_list.emplace_back(fmt::format("{:d}-{:d}", a + 1, b + 1));
            }
            a = nums[i];
            b = a;
        }
    }
    // handle last
    if (a == b) {
        ret_list.emplace_back(fmt::format("{:d}", a + 1));
    } else {
        ret_list.emplace_back(fmt::format("{:d}-{:d}", a + 1, b + 1));
    }
    // join ret_list to ret
    std::string ret;
    if (!ret_list.empty()) ret = ret_list[0];
    for (size_t i = 1; i < ret_list.size(); ++ i) {
        ret += "," + ret_list[i];
    }
    return ret;
}

