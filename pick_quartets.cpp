#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <algorithm>
#include <functional>



int most_present_rule(int x0, int x1, int x2, int n_trees){
    // returns the most present index 0..2
    std::vector<int> compare{x0,x1,x2};
    std::vector<int>::iterator result;
    result = std::max_element(compare.begin(), compare.end());
    int count = std::count(compare.begin(), compare.end(), *result);
    // std::cout << "COUNT" << count << std::endl;
    if (count==1){
        // std::cout << "HERE" << std::endl;
        return std::distance(compare.begin(), result);
    }
    return -1;
}

int majority_rule(int x0, int x1, int x2, int n_trees){
    int half = n_trees/2;
    if (x0>half){
        return 0;
    }else if (x1 > half){
        return 1;
    }else if (x2 > half){
        return 2;
    } 
    return -1;
}

std::string create_quartfile(std::function<int(int,int,int,int)> picking_rule, const std::string inputFileName){
    std::ifstream inputFile(inputFileName); // Replace "input.txt" with the actual file name.
    std::string output="";
    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the input file." << std::endl;
        std::exit(1);
    }

    std::string first;
    getline(inputFile,first);
    std::istringstream first_iss(first);
    int n, n_trees;
    first_iss >> n >> n_trees;
    output += std::to_string(n) + "\n";

    // Read and store the data lines
    for (int i = 0; i < n; i++) {
        std::string tmp;
        getline(inputFile, tmp);
        output += tmp + "\n";
    }

    output += "\n";


    std::string line; std::string line2;

    // Read and process the pattern lines
    while (std::getline(inputFile, line)) {
        if (line.empty()) {
            continue; // Skip empty lines
        }
        std::stringstream ss{line};
        std::vector<std::string> v;
        for (int i=0; i<3; i++){
            std::string a;
            std::getline(ss, a, ';');
            const auto strBegin = a.find_first_not_of(" ");
            v.push_back(a.substr(strBegin));
        }

        if (!getline(inputFile, line2) || line.empty()){
            std::cerr << "ERROR: Occurrences not given for last line" << std::endl;
            std::exit(1);
        }

        std::istringstream iss(line2);
        int c,d,e;
        iss >> c >> d >> e;
        int max_idx = picking_rule(c,d,e,n_trees);
        if (max_idx == -1){
            continue;
        }

        std::string pattern = v[max_idx];
        output += pattern + "\n";
    }

    inputFile.close(); // Close the input file
    return output;
}


int main(int argc, char* argv[]) {
    std::string inputFileName="NULL";
    std::string outputFileName="NULL";
    std::string pickingRule = "most-present";
    std::function<int(int,int,int,int)> pick= most_present_rule;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        // Check for specific options
        if (arg == "--input" && i + 1 < argc) {
            inputFileName = argv[i + 1];
            std::cout << "Input file: " << inputFileName << std::endl;
            i++;  // Skip the next argument
        } else if (arg == "--output" && i + 1 < argc) {
            outputFileName = argv[i + 1];
            std::cout << "Output file: " << outputFileName << std::endl;
            i++;  // Skip the next argument
        } else if (arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [--input input_file] [--output output_file]" << std::endl;
            return 0;
        }else if (arg == "--picking-rule" && i + 1 < argc) {
            pickingRule = argv[i + 1];
            std::cout << "Picking Rule: " << pickingRule << std::endl;
            i++;  // Skip the next argument
        } 
        else {
            // Handle unknown arguments or perform other logic
            std::cerr << "Unknown argument: " << arg << std::endl;
            return 1;
        }
    }
    if (inputFileName=="NULL" || outputFileName == "NULL"){
        std::cerr << "Please specify inputFile and outputFile" << std::endl;
        std::exit(1);
    }
    if (pickingRule=="most-present"){
        pick = most_present_rule;
    }else if (pickingRule == "majority"){
        pick = majority_rule;
    }else{
        std::cerr << "Unknown picking rule: " << pickingRule << std::endl;
        return 1;
    }


    std::string output = create_quartfile(pick, inputFileName);
    std::ofstream ofs(outputFileName);
    if (!ofs) {
        std::cerr << "ERROR: Failed to open output file." << std::endl;
        std::exit(1);
    }

    ofs << output;

    return 0;
}
