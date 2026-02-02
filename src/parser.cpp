#include "parser.h"
#include "Logger.h"
#include "utilities.h"

namespace CppBOLOS {

// std::stringstream clean_file(const std::string& filename) {
//     std::ifstream file(filename);
//     // Check if the file was opened successfully
//     if (!file.is_open()) {
//         parserError("CppBOLOS parserError: Could not open file " + filename);
//     }
//
//     // Read each line, remove carriage return, and store the cleaned lines
//     std::stringstream ss;
//     std::string line;
//     while (std::getline(file, line)) {
//         // Remove carriage return and write directly to stringstream
//         if (!line.empty() && line.back() == '\r') {
//             line.pop_back(); // Efficiently removes trailing '\r'
//         }
//         ss << line << "\n";
//     }
//
//     return ss;
// }

std::vector<std::string> clean_file(const std::string& filename) {
    std::ifstream file(filename);

    // Check if the file was opened successfully
    if (!file.is_open()) {
        parserError("CppBOLOS parserError: Could not open file " + filename);
    }

    // Read each line, remove carriage return, and store the cleaned lines
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line)) {
        // Remove carriage return (if present)
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        lines.push_back(line);
    }

    return lines;
}

// This is to mirror split() in python.
// it will handle multiple consecutive delimiters properly by ignoring empty results between them
std::vector<std::string> split_string(const std::string& str, const std::string& delim = "[ \t]+") {
    std::regex re(delim);
    std::sregex_token_iterator first{str.begin(), str.end(), re, -1}, last;
    std::vector<std::string> result;
    for (auto it = first; it != last; ++it) {
        std::string item = it->str();
        if (!item.empty()) {
            result.push_back(item);
        }
    }
    return result;
}

std::vector<std::string> read_until_sep(std::istream& file) {
    std::vector<std::string> lines;
    std::string line;
    std::regex sep("-----+");
    while (std::getline(file, line)) {
        if (std::regex_match(line, sep)) {
            break;
        }
        lines.push_back(line);
    }
    return lines;
}


std::string join_strings(const std::vector<std::string>& strings, const std::string& delimiter) {
    std::ostringstream joined;
    copy(strings.begin(), strings.end(), std::ostream_iterator<std::string>(joined, delimiter.c_str()));
    std::string result = joined.str();
    return result.substr(0, result.size() - delimiter.size());  // Remove the last delimiter
}

std::string remove_carriage_return(std::string line) {
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
    return line;
}

Block read_block(std::istream& file, bool has_arg) {
    Block block;

    // Read the target line
    if (!std::getline(file, block.target)) {
        parserError("Failed to read target line");
    }
    block.target = trim(block.target);
    // If has_arg is true, read the argument line
    if (has_arg) {
        if (!std::getline(file, block.arg)) {
            parserError("Failed to read argument line");
        }
    }

    // Read the comment lines until a separator line
    block.comment = join_strings(read_until_sep(file), "\n");

    // Read the data lines until a separator line
    std::vector<std::string> data_lines = read_until_sep(file);

    for (const std::string& line : data_lines) {
        std::vector<std::string> items = split_string(line);
        std::vector<double> row;
        for (const std::string& item : items) {
            std::stringstream item_ss(item);
            double item_double;
            item_ss >> item_double;
            row.push_back(item_double);
        }
        if (row[0] < 0 or row[1] < 0) parserError("Negative data in ", line);
        // LOG_DEBUG(row[0] <<" \t" << row[1]);
        block.data.push_back(row);
    }
    /*
    // Remove leading zeros
    while (block.data.size() > 1 && block.data.front()[1] == 0.0 && block.data[1][1] == 0.0) {
        LOG_DEBUG("Removing leading zero line: " << block.data.front()[0] << "\t" << block.data.front()[1]);
        block.data.erase(block.data.begin());
    }
    // Remove trailing zeros
    while (block.data.size() > 1 && block.data.back()[1] == 0.0 && block.data[block.data.size() - 2][1] == 0.0) {
        LOG_DEBUG("Removing trailing zero line: " << block.data.back()[0] << "\t" << block.data.back()[1]);
        block.data.pop_back();
    }
    // Remove zeros in the middle
    for (auto it = block.data.begin() + 1; it != block.data.end() - 1;) {
        if ((*it)[1] == 0.0) {
            LOG_DEBUG("Removing middle zero line: " << (*it)[0] << "\t" << (*it)[1]);
            it = block.data.erase(it);
        } else {
            ++it;
        }
    }
    */

    return block;
}

Collision read_momentum(std::istream& file) {
    Block block = read_block(file, true);
    Collision process;
    process.target = trim(split_string(block.target, "->")[0]);
    process.kind = "MOMENTUM";
    process.product = "";
    process.comment = block.comment;
    std::vector<std::string> arg_items = split_string(block.arg);
    process.mass_ratio = std::stod(arg_items[0]);  // Mass ratio is the first item of the argument line
    process.threshold = 0.0;  // No threshold for momentum processes
    process.weight_ratio = -1.0;  // No weight ratio for momentum processes
    process.data = block.data;
    return process;
}

Collision read_excitation(std::istream& file) {
    Block block = read_block(file, true);
    Collision process;
    process.kind = "EXCITATION";
    std::vector<std::string> target_items = split_string(block.target, "<->");
    std::vector<std::string> arg_items = split_string(block.arg);
    if (target_items.size() == 2) {
        process.threshold = std::stod(arg_items[0]);  // Threshold is the first item of the argument line
        process.weight_ratio = std::stod(arg_items[1]);  // Weight ratio is the second item of the argument line
    } else {
        target_items = split_string(block.target, "->");
        process.threshold = std::stod(arg_items[0]);  // Argument line only contains threshold
        process.weight_ratio = 1.0; }
    process.target = trim(target_items[0]);
    process.product = trim(target_items[1]);
    process.mass_ratio = -1.0; // No mass ratio
    process.data = block.data;
    return process;
}

Collision read_attachment(std::istream& file) {
    Block block = read_block(file, false);
    Collision process;
    process.kind = "ATTACHMENT";
    std::regex re("<->|->");
    std::sregex_token_iterator it(block.target.begin(), block.target.end(), re, -1);
    std::vector<std::string> target_items(it, {});
    if (target_items.size() == 2) {
        process.target = trim(target_items[0]);
        process.product = trim(target_items[1]);
    } else {
        process.target = trim(target_items[0]);
        process.product = ""; }
    process.mass_ratio = -1.0; //
    process.weight_ratio = -1.0;  // No weight ratio for attachment processes
    process.threshold = 0.0;
    process.comment = block.comment;
    process.data = block.data;
    return process;
}

// Function type for process-reading functions
using ProcessFunc = Collision (*)(std::istream&);

// KEYWORDS map
std::unordered_map<std::string, ProcessFunc> KEYWORDS = {
        {"MOMENTUM", read_momentum},
        {"ELASTIC", read_momentum},
        {"EFFECTIVE", read_momentum},
        {"EXCITATION", read_excitation},
        {"IONIZATION", read_excitation},
        {"ATTACHMENT", read_attachment}
};

std::vector<Collision> parse(std::istream& file) {
    std::vector<Collision> collisions;
    std::string line;
    while (std::getline(file, line)) {
        auto it = KEYWORDS.find(line);
        if (it != KEYWORDS.end()) {
            LOG_DEBUG("New process of type: '" << line << "'");
            Collision collision = it->second(file);
            collision.kind = line;
            collisions.push_back(collision);
        }
    }
    LOG_INFO("****** Parsing complete with " << collisions.size() << " collisions/processes read. ******\n")
    return collisions;
}

std::vector<Collision> parse(const std::vector<std::string>& lines) {
    // Combine the vector of strings into a single string
    std::string combined;
    for (const auto& line : lines) {
        combined += line + "\n"; // Add each line followed by a newline
    }

    // Wrap the combined string in an istringstream
    std::istringstream stream(combined);

    // Now use the existing stream-based parsing logic
    std::vector<Collision> collisions;
    std::string line;
    while (std::getline(stream, line)) {
        auto it = KEYWORDS.find(line);
        if (it != KEYWORDS.end()) {
            LOG_DEBUG("New process of type: '" << line << "'");
            Collision collision = it->second(stream); // Call the process function
            collision.kind = line;
            collisions.push_back(collision);
        }
    }

    LOG_INFO("****** Parsing complete with " << collisions.size() << " collisions/processes read. ******\n");
    return collisions;
}

}