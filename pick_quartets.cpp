#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <sstream>

using namespace std;

int main() {
    ifstream inputFile("examples/example_quartcount"); // Replace "input.txt" with the actual file name.
    std::string output="";
    if (!inputFile.is_open()) {
        cerr << "Failed to open the input file." << endl;
        return 1;
    }

    std::string n;
    getline(inputFile,n);
    output += n + "\n";

    // Read and store the data lines
    for (int i = 0; i < std::stoi(n); i++) {
        std::string tmp;
        getline(inputFile, tmp);
        output += tmp + "\n";
    }

    output += "\n";


    string line; string line2;

    // Read and process the pattern lines
    while (getline(inputFile, line)) {
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
        int max_idx = -1;
        if (d>e){
            if (c > d){
                max_idx = 0;
            }else if(d>c){
                max_idx = 1;
            }
        }else if (e > d){
            if (c > e){
                max_idx = 0;
            }else if(e>c){
                max_idx = 2;
            }
        }else{
            if (c>d){
                max_idx = 0;
            }
        }
        if (max_idx == -1){
            continue;
        }

        std::string pattern = v[max_idx];
        output += pattern + "\n";
    }

    std::cout << output;

    inputFile.close(); // Close the input file

    return 0;
}
