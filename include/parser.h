#include <fstream>
#include <regex>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <iterator>
#include <iostream>

namespace CppBOLOS {

struct Collision {
    std::string target;
    std::string kind;
    std::vector<std::vector<double>> data;
    std::string comment;
    double mass_ratio;
    std::string product;
    double threshold;
    double weight_ratio;
};

struct Block {
    std::string target;
    std::string arg;
    std::string comment;
    std::vector<std::vector<double>> data;
};

// std::stringstream clean_file(const std::string& filename);
std::vector<std::string> clean_file(const std::string& filename);

std::vector<std::string> split_string(const std::string& str, const std::string& delim);
std::vector<std::string> read_until_sep(std::istream& file);
std::string join_strings(const std::vector<std::string>& strings, const std::string& delimiter);
std::string remove_carriage_return(std::string line);

Block read_block(std::istream& file, bool has_arg);

Collision read_momentum(std::istream& file);

Collision read_excitation(std::istream& file);

Collision read_attachment(std::istream& file);

std::vector<Collision> parse(std::istream& file);
std::vector<Collision> parse(const std::vector<std::string>& lines);

template <typename... Args>
inline void parserError(Args... args) {
    std::ostringstream oss;
    (oss << ... << args);  // fold expression
    throw std::runtime_error("CppBOLOS::parserError: " + oss.str());
}

}