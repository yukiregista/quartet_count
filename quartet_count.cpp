#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <chrono>
#include <thread>
#include <cmath>
#include <fstream>
#include <array>
#include <cassert>
#include <functional>




int BinomialCoefficient(const int n, const int k) {
  if (k>n){
    return 0;
  }
  if (k==0){
    return 1;
  }
  std::vector<int> aSolutions(k,0);
  aSolutions[0] = n - k + 1;

  for (int i = 1; i < k; ++i) {
    aSolutions[i] = aSolutions[i - 1] * (n - k + 1 + i) / (i + 1);
  }

  return aSolutions[k - 1];
}

int Combinations_to_Unique_Number(std::array<int,4>& combinations,const int& n){
  // IN this version, array combinations should have the number of 1..n (excluding zero)
  std::sort (combinations.begin(), combinations.begin()+4);
  int number = 0;
  for (int i=0; i<4; i++){
    number += BinomialCoefficient(n-combinations[i], 4-i);
  }
  return number;
}

class Taxa {
private:
    int nextId = 1;
    int nextInternalId = -1;
    std::vector<std::string> taxonNames{};
    std::unordered_map<std::string, int> taxonIdMap{};

public:
    const std::vector<std::string>& getTaxonNames() const { return taxonNames; }

    // Method to add a taxon
    void addTaxon(const std::string& name) {
        taxonNames.push_back(name);
        taxonIdMap[name] = nextId++;
    }


    void addInternalTaxon(){
        nextInternalId--;
    }

    int getTaxonId(const std::string& name) const {
        if (taxonIdMap.size()==0){
            return -1;
        }
        auto it = taxonIdMap.find(name);
        return (it != taxonIdMap.end()) ? it->second : -1;
    }
    int getNextId() const {return nextId;}
    int getNextInternalId() const {return nextInternalId;}

};


class Node {
public:
    int getNodeId() const {return NodeId;}
    int getTaxonId() const { return taxonId; }
    double getBranchLength() const { return branchLength; }

    void setParent(Node* par){parent = par;}
    Node* getParent(){return parent;}

    const std::unordered_map<int, Node*>& getChildren() const { 
        return children; 
        }

    // Add existing node as a child
    void addChild(Node* child) {
        child->setParent(this);
        children[child->getNodeId()] = child;
    }

    bool is_leaf() const {return children.size()==0;}

    void set_terminal_taxon_IDs(std::set<int> terminal){this->terminal_taxon_IDs = terminal;}
    void update_terminal_taxon_IDs(std::set<int> terminal){this->terminal_taxon_IDs.insert(terminal.begin(), terminal.end());}
    void clear_terminal_taxon_IDs(){this->terminal_taxon_IDs.clear();}
    std::set<int> getTerminal_taxon_IDs() const {return terminal_taxon_IDs;}
    std::unordered_map<int, Node*> children{};
    std::set<int> terminal_taxon_IDs{};
    double branchLength=std::nan(""); // Stores branch length or NaN if not specified
    

private:
    int taxonId;
    int NodeId;
    Node* parent=nullptr;

    // Constructor with associated Taxon ID and branch length (private)

    // friend class Tree;        // Tree class can access the private constructor
    // friend class NewickParser; // NewickParser class can access the private constructor

public:
    // Factory method to create a Node with associated taxon ID and branch length
    static std::unique_ptr<Node> createNode(int associatedTaxonId = 0, double branchLength = std::nan("")) {
        return std::make_unique<Node>(associatedTaxonId, branchLength);
    }
    Node(int associatedTaxonId, double branchLength) : branchLength(branchLength),taxonId(associatedTaxonId), NodeId(associatedTaxonId){}
};

void set_terminal_taxon_IDs(Node& node) {
    if (node.is_leaf()){
        node.set_terminal_taxon_IDs(std::set<int>{node.getTaxonId()});
    } else{
        for (auto child : node.getChildren()){
            std::set<int> tmp = child.second->getTerminal_taxon_IDs();
            node.update_terminal_taxon_IDs(tmp);
        }
    }
}

void reset_terminal_taxon_IDs(Node& node) {
    node.clear_terminal_taxon_IDs();
    if (node.is_leaf()){
        node.set_terminal_taxon_IDs(std::set<int>{node.getTaxonId()});
    } else{
        for (auto child : node.getChildren()){
            std::set<int> tmp = child.second->getTerminal_taxon_IDs();
            node.update_terminal_taxon_IDs(tmp);
        }
    }
}

void count_quartets_of_the_node(const Node& node, std::vector<int>& quartet_counts, const std::set<int>& all_ids,const int& root_id,const int& n_taxa){
    // count the quartets associated with this node
   std::unordered_map<int, Node*> children = node.getChildren();
   int children_size = children.size();
   if (children_size <=1){
        // either leaf or unification. no need for updateing
        return ;
   }
   //std::cout << "children size bigger" << std::endl;
    std::set<int> terminals = node.getTerminal_taxon_IDs();
    //std::cout << "all_id:";
    // for(auto item: all_ids){std::cout<< item;}
    // std::cout << std::endl;
    // std::cout << "terminals:";
    // for(auto item: terminals){std::cout  << item;}
    // std::cout << std::endl;
    std::set<int> others{}; // will contain other terminals (root is not considered as terminals)
    std::set_difference(all_ids.begin(), all_ids.end(), terminals.begin(), terminals.end(), 
                        std::inserter(others, others.cbegin()));
    // std::cout << "others:";
    int size_of_others = others.size();
    if (size_of_others==0){return;}
    // std::cout << "YES" << std::endl;


    // for (auto child1 = children.begin(); child1 != children.end(); child1++){
    //     std::cout << child1->first << std::endl;
    // }
    // std::cout << "YESS" << std::endl;

    for (auto child1 = children.begin(); child1 != children.end(); child1++) {
        //access all terminal of this child.
        //std::cout << "YES2!" << std::endl;
        if (child1 == children.end()){break;}
        for (int term1 : child1->second->getTerminal_taxon_IDs()){
            // term 1 is picked from child 1
            for (auto child2 = std::next(child1); child2 != children.end(); child2++) {
                for (int term2 : child2->second->getTerminal_taxon_IDs()){
                    // term2 is picked from child 2
                    for (int other_term : others){
                        // (term1, term2) | (root_id, other_term) is present.
                        // First, make an array of these four elements.
                        std::array<int,4> quartet = {term1, term2, other_term, root_id};
                        // std::cout << "term1" << term1 << std::endl;
                        // std::cout << "term2" << term2 << std::endl;
                        // std::cout << "other_term" << other_term << std::endl;
                        // std::cout << "root_id" << root_id << std::endl;
                        //get the unique number of this quartet
                        int quartet_num = 3 * Combinations_to_Unique_Number(quartet, n_taxa);
                        // determine which pattern it is.
                        std::pair<int, int> result1 = std::minmax(term1, term2);
                        std::pair<int, int> result2 = std::minmax(root_id, other_term);
                        int type = 0;
                        if (result1.second < result2.first || result2.second < result1.first){
                            // type 0 (1,2)|(3,4)
                            ;
                        }else if ((result1.first < result2.first && result1.second < result2.second) || (result1.first > result2.first && result1.second > result2.second)){
                            // type 1: (1,3)|(2,4)
                            type = 1;
                        }else{
                            type = 2;
                        }
                        int quartet_topology_id = quartet_num + type;
                        quartet_counts[quartet_topology_id] += 1;
                    }
                }
            }
        }
    }
}


class Tree {
public:
    Node& getRoot() const { return *root; }
    const Taxa& getTaxa() const {return *taxa;}

    // Function to print an ASCII plot of the phylogenetic tree
    void printAsciiTree(const Node& node, int level = 0) const {
        // Print the current node
        for (int i = 0; i < level; ++i) {
            std::cout << "| ";
        }
        std::cout << "+-" << "Taxon ID: " << node.getNodeId() << "\n";

        // Recursively print the children
        for (const auto pair : node.getChildren()) {
            printAsciiTree(*pair.second, level + 1);
        }
    }

    void reroot_with_leaf_ind(int leaf_ind);

    void remove_unificating_root();


    // Function to perform post-order traversal on the tree and apply a function to each node
    template <typename Function>
    void postOrderTraversal_change(Function&& func, Node& node) const {
        for (const auto& child : node.getChildren()) {
            postOrderTraversal_change(std::forward<Function>(func), *child.second);  // Recursively apply to child nodes
        }

        func(node);  // Apply the function to the current node
    }
        // Function to perform post-order traversal on the tree and apply a function to each node
    template <typename Function>
    void postOrderTraversal(Function&& func, const Node& node) const {
        for (const auto& child : node.getChildren()) {
            postOrderTraversal(std::forward<Function>(func), *child.second);  // Recursively apply to child nodes
        }

        func(node);  // Apply the function to the current node
    }

    template <typename Function, typename T>
    void postOrderDFS_modifyData(Function&& func, const Node& node, T& dataContainer) const {
        for (const auto& child : node.getChildren()) {
            postOrderDFS_modifyData(std::forward<Function>(func), *child.second, dataContainer);  // Recursively apply to child nodes
        }
        func(node, dataContainer);  // Apply the function to the current node
    }

    // General function to perform post-order traversal on the tree and apply a function to each node
    template <typename Function, typename... Args>
    void postOrderTraversal(Function&& func,const Node& node, Args&&... args) const {
        for (const auto& child : node.getChildren()) {
            postOrderTraversal(std::forward<Function>(func), *child.second, std::forward<Args>(args)...);
        }

        func(node, std::forward<Args>(args)...);  // Apply the function to the current node
    }

    void postOrderQuartetCounting(const Node& node, std::vector<int>& quartet_counts, const std::set<int>& all_ids,const int& root_id, int& n_taxa){
        for (const auto& child : node.getChildren()) {
            postOrderQuartetCounting(*child.second, quartet_counts, all_ids,root_id,n_taxa);
        }
        count_quartets_of_the_node(node, quartet_counts, all_ids,root_id,n_taxa);
    }

    Node* root=nullptr;
    Taxa* taxa=nullptr;
    std::unordered_map<int,std::unique_ptr<Node>> nodes;
};


void Tree::reroot_with_leaf_ind(int leaf_ind){
    std::unique_ptr<Node>& leaf = nodes[leaf_ind];
    if (root==leaf.get()){
        // do nothing
        return;
    }
    Node* chi = leaf.get();
    Node* par = chi->getParent();
    while (chi!=root){
        Node* nextpar = par->getParent();
        chi->addChild(par); // add par as child of par
        par->setParent(chi); // set chi as parent of par
        par->children.erase(chi->getNodeId()); //delete chi from children of par
        chi = par;
        par = nextpar;
    }
    leaf->setParent(nullptr);
    this->root = leaf.get();
//     std::this_thread::sleep_for(std::chrono::milliseconds(1000));
//     for (auto child :nodes[-1]->children){
//         std::cout << child.second->NodeId << std::endl;
//     }
//     std::this_thread::sleep_for(std::chrono::milliseconds(1000));

}

void Tree::remove_unificating_root(){
    //remove unificating root
    if (root->children.size()!=1){
        std::cout << "children of root is not one. do nothing." << std::endl;
        return;
    }
    Node* newRoot;
    while (root->children.size()==1){
        for (auto child : root->children){
            child.second->setParent(nullptr);
            newRoot = child.second;
            break;
        }
        int rootId = root->getNodeId();
        if (rootId > 0){ nodes.erase(rootId); }
        root = newRoot;
    }
}

class NewickParser {
public:
    std::vector<std::unique_ptr<Tree>> parseNewick(const std::string& newickString,std::unique_ptr<Taxa>& taxa);

private:
    Node* parseInternal(std::string::const_iterator& it, const std::string& newickString,std::unique_ptr<Tree>& tree);

    void skipWhitespace(std::string::const_iterator& it, const std::string& newickString);

    std::string readUntil(std::string::const_iterator& it, const std::set<char>& delimiters, const std::string& newickString);

    double readBranchLength(std::string::const_iterator& it, const std::string& newickString);
};

std::vector<std::unique_ptr<Tree>> NewickParser::parseNewick(const std::string& newickString, std::unique_ptr<Taxa>& taxa) {
    std::string::const_iterator it = newickString.begin();
    std::vector<std::unique_ptr<Tree>> trees;
    while (it != newickString.end()){
        std::unique_ptr<Tree> tree = std::make_unique<Tree>();
        tree->taxa = taxa.get();
        Node* root = parseInternal(it, newickString, tree);
        tree->root = root;
        //taxa = tree->taxa;
        trees.push_back(std::move(tree));
    }
    return trees;
}

Node* NewickParser::parseInternal(std::string::const_iterator& it, const std::string& newickString, std::unique_ptr<Tree>& tree) {
    if (*it == '(') {
        // Internal node
        ++it; // Skip '('
        tree->taxa->addInternalTaxon();
        std::unique_ptr<Node> internal1 = Node::createNode(tree->taxa->getNextInternalId() + 1);
        while (*it != ')') {
            Node* tmp = parseInternal(it, newickString, tree);
            internal1->addChild(tmp);
            skipWhitespace(it, newickString);
            if (*it == ',') {
                ++it; // Skip ','
            }
        }
        // std::cout << *it << std::endl;  
        // std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        ++it; // Skip ')'
        double branchLength;
        if (*it == ':'){
            ++it;
            branchLength = readBranchLength(it, newickString);
        } else if (*it == ','){
            ++it;
            branchLength = std::nan("");
        } else if (*it == ')'){
            // do nothing
            branchLength = std::nan("");
        } else {
            // internal node is labeled -- but we don't store that information.
            std::string taxonName = readUntil(it, { ':', ',', ')' }, newickString);
            if (*it == ':'){
                ++it;
                branchLength = readBranchLength(it, newickString);
            } else if (*it == ','){
                ++it;
                branchLength = std::nan("");
            }else{
                branchLength = std::nan("");
            }
        }
        internal1->branchLength = branchLength;
        skipWhitespace(it, newickString);
        int nodeId = internal1->getNodeId();
        tree->nodes[nodeId] = std::move(internal1);
        return tree->nodes[nodeId].get(); // Internal node has no associated taxon
    } else {
        // Leaf node
        std::string taxonName = readUntil(it, { ':', ',', ')' }, newickString);
        double branchLength;
        if (*it == ':'){
            ++it; //skip the delimiter
            branchLength = readBranchLength(it, newickString);
            if (*it == ','){
                //skip the delimiter
                ++it;
            }
        } else if (*it == ','){
            ++it; //skip the delimiter
            branchLength = std::nan("");
        } else if (*it == ')'){
            branchLength = std::nan("");
        }
        skipWhitespace(it, newickString);
        // add taxon
        int taxonId = tree->taxa->getTaxonId(taxonName);
        if (taxonId == -1){
            tree->taxa->addTaxon(taxonName);
            taxonId = tree->taxa->getNextId()-1;
        }
        tree->nodes[taxonId] = std::make_unique<Node>(taxonId, branchLength);
        return tree->nodes[taxonId].get();
    }
}


void NewickParser::skipWhitespace(std::string::const_iterator& it, const std::string& newickString) {
    while (it != std::end(newickString) && std::isspace(*it)) {
        ++it;
    }
}

std::string NewickParser::readUntil(std::string::const_iterator& it, const std::set<char>& delimiters, const std::string& newickString) {
    std::string result;
    while (it != std::end(newickString) && delimiters.find(*it) == delimiters.end() && !std::isspace(*it)) {
        result += *it;
        ++it;
    }

    return result;
}

double NewickParser::readBranchLength(std::string::const_iterator& it, const std::string& newickString) {
    skipWhitespace(it, newickString);
    std::string branchLengthStr = readUntil(it, {',', ')'}, newickString);
    return std::stod(branchLengthStr);
}

std::vector<std::unique_ptr<Tree>> readNewickFile(const std::string& fileName, std::unique_ptr<Taxa>& taxa){
    std::string newickString;
    {
        std::ifstream ifs(fileName);
        if (!ifs)
        {
            std::cerr << "ERROR: Could not open input file." << std::endl;
            std::exit(1);
            return std::vector<std::unique_ptr<Tree>>{};
        }
        
        std::string buf;
        while (!ifs.eof())
        {
            std::getline(ifs, buf);
            newickString += buf + "\n";
        }
    }
        // parse newick file
        NewickParser newickParser;
        std::vector<std::unique_ptr<Tree>> phylogeneticTree = newickParser.parseNewick(newickString, taxa); //ok
    

    return phylogeneticTree;
}


void printTerminalId(const Node& node) {
    std::cout << "Taxon ID: " << node.getTaxonId() << std::endl;
    for (int id : node.getTerminal_taxon_IDs()){
        std::cout << "Terminal ID: " << id << std::endl;
    }
    std::cout << std::endl;
}


std::vector<int> count_quartets(std::vector<std::unique_ptr<Tree>>& trees){
    // count the number of quartets
    // From c++11, the vector will be returned using move semantics.
    int n_taxa = trees[0]->taxa->getNextId()-1;
    int nC4 = BinomialCoefficient(n_taxa,4);
    std::vector<int> quartet_counts(nC4*3, 0);
    int count = 0;
    for (std::unique_ptr<Tree> & tree : trees){
        std::cout << "Tree" << count++ << std::endl;
        // count quartets in the tree.
        // Iteration over leaf nodes
        for (int reroot_index=1; reroot_index <= n_taxa-3; reroot_index++){
            //std::cout << reroot_index <<std::endl;
            const int root_id = reroot_index;
            // 1. Rerooting the tree with the leaf
            tree->reroot_with_leaf_ind(reroot_index); //rerooting; changing parent-child structure
            // tree->postOrderTraversal(printTerminalId, *tree->root);
            tree->postOrderTraversal_change(reset_terminal_taxon_IDs, *tree->root); // reset terminals
            // tree->postOrderTraversal(printTerminalId, *tree->root);
            // postOrderTraversal to count all quartets inccluding reroot_index
            //tree->postOrderTraversal(printTerminalId, *tree->root);
            std::set<int> all_ids;
            for (auto child: tree->root->children){
                all_ids = child.second->getTerminal_taxon_IDs();
                break;
            }
            tree->postOrderQuartetCounting(*tree->root, quartet_counts,all_ids, root_id, n_taxa);
            // finally remove root from the tree
            tree->remove_unificating_root();
        }   
    }

    // for (int i=0; i<nC4; i++){
    //     int all = quartet_counts[3*i] + quartet_counts[3*i+1] + quartet_counts[3*i+2];
    //     if (all != trees.size()){
    //         std::cout << all << std::endl;
    //     }
    // }//check complete!

    return quartet_counts;
}

int most_present_rule(int x0, int x1, int x2){
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

std::string create_quartfile(std::function<int(int,int,int)> picking_rule, const std::vector<int>& quart_counts, const std::vector<std::string>& taxonNames, const int n_taxa){
    std::string quartfile_str = std::to_string(n_taxa) + "\n";
    // First write index to taxa correspondence.
    for (int i=1; i<= n_taxa; i++){
        quartfile_str += std::to_string(i) + " " + taxonNames[i-1] + "\n";
    }
    quartfile_str += "\n";

    // Write suitable quartets
    int count = 0;
    for (int i=n_taxa-3; i>=1; i--){
        for (int j=n_taxa-2; j>i; j--){
            for (int k=n_taxa-1; k>j; k--){
                for (int l=n_taxa; l>k; l--){
                    // std::cout << "QCOUNTS" << quart_counts[3*count] <<  quart_counts[3*count+1] << quart_counts[3*count+2] << std::endl;
                    int type = picking_rule(quart_counts[3*count], quart_counts[3*count+1],quart_counts[3*count+2]);
                    // std::cout << "TYPE: " << type << std::endl;
                    if (type==0){
                        //type 0 (1,2) | (3,4)
                        quartfile_str += std::to_string(i) + " " + std::to_string(j) + " | " + std::to_string(k) + " " + std::to_string(l) + "\n";
                    }else if(type==1){
                        //type 1 (1,3) | (2,4)
                        quartfile_str += std::to_string(i) + " " + std::to_string(k) + " | " + std::to_string(j) + " " + std::to_string(l) + "\n";
                    }else if(type==2){
                        //type 1 (1,4) | (2,3)
                        quartfile_str += std::to_string(i) + " " + std::to_string(l) + " | " + std::to_string(j) + " " + std::to_string(k) + "\n";
                    }
                    // if other than 0,1,2, do nothing.
                    count++;
                }
            }
        }
    }
    //std::cout << quartfile_str << std::endl;
    return quartfile_str;
}


int check_main(){
// std::unique_ptr<Taxa> taxa = std::make_unique<Taxa>();
    // NewickParser newickParser;
    // std::vector<std::unique_ptr<Tree>> phylogeneticTree = newickParser.parseNewick("((A,B),(C,D),E,F);\n ((A,B),(C,E),(D,F))", taxa); //ok
    // for (int i=0; i<phylogeneticTree.size(); i++){
    //     phylogeneticTree[i]->printAsciiTree(*(phylogeneticTree[i]->root), 0);
    //     std::cout << std::endl;
    // }
    // phylogeneticTree[0]->nodes[-3]->getTerminal_taxon_IDs();
    // phylogeneticTree[0]->postOrderTraversal_change(set_terminal_taxon_IDs, *phylogeneticTree[0]->root);
    // std::cout << "a" << std::endl;

    // phylogeneticTree[0]->postOrderTraversal_change(set_terminal_taxon_IDs, *phylogeneticTree[0]->root);//ok
    // for (auto item : phylogeneticTree[0]->root->getTerminal_taxon_IDs()){
    //     std::cout << item << std::endl;
    // }
    // phylogeneticTree[0]->postOrderTraversal(printTerminalId, *phylogeneticTree[0]->root);//ok

    // // rerooting
    // phylogeneticTree[0]->reroot_with_leaf_ind(1);
    // phylogeneticTree[0]->printAsciiTree(*(phylogeneticTree[0]->root), 0);
    // // reset terminals
    // phylogeneticTree[0]->postOrderTraversal_change(reset_terminal_taxon_IDs, *phylogeneticTree[0]->root);
    // phylogeneticTree[0]->postOrderTraversal(printTerminalId, *phylogeneticTree[0]->root);//ok


    // Check implementation of Combinations_to_unique_number ->OK
    // int n_taxa = 50; int prev = -1;
    // for (int i=n_taxa-3; i>=1; i--){
    //     for (int j=n_taxa-2; j>i; j--){
    //         for (int k=n_taxa-1; k>j; k--){
    //             for (int l=n_taxa; l>k; l--){
    //                 std::array<int,4> check{i,j,k,l};
    //                 int now = Combinations_to_Unique_Number(check, n_taxa);
    //                 if (now-prev !=1){std::cout << "error";}
    //                 prev = now;
    //             }
    //         }
    //     }
    // }

    return 0;
}

int main(int argc, char* argv[]){
    std::string inputFile="NULL";
    std::string outputFile="NULL";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        // Check for specific options
        if (arg == "--input" && i + 1 < argc) {
            inputFile = argv[i + 1];
            std::cout << "Input file: " << inputFile << std::endl;
            i++;  // Skip the next argument
        } else if (arg == "--output" && i + 1 < argc) {
            outputFile = argv[i + 1];
            std::cout << "Output file: " << outputFile << std::endl;
            i++;  // Skip the next argument
        } else if (arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [--input input_file] [--output output_file]" << std::endl;
            return 0;
        } else {
            // Handle unknown arguments or perform other logic
            std::cerr << "Unknown argument: " << arg << std::endl;
            return 1;
        }
    }

    if (inputFile=="NULL" || outputFile == "NULL"){
        std::cerr << "Please specify inputFile and outputFile" << std::endl;
        std::exit(1);
    }

    
    //file IO
    std::unique_ptr<Taxa> taxa = std::make_unique<Taxa>();
    std::vector<std::unique_ptr<Tree>> phylogeneticTree =readNewickFile(inputFile, taxa);
    //phylogeneticTree[0]->printAsciiTree(*phylogeneticTree[0]->root);
    //phylogeneticTree[0]->reroot_with_leaf_ind(1);
    phylogeneticTree[0]->postOrderTraversal_change(reset_terminal_taxon_IDs, *phylogeneticTree[0]->root);
    //phylogeneticTree[0]->postOrderTraversal(printTerminalId, *phylogeneticTree[0]->root);//ok



    std::vector<int> quartet_counts = count_quartets(phylogeneticTree);


    const std::vector<std::string>& taxonNames = taxa->getTaxonNames();
    int n_taxa = taxa->getNextId() - 1;
    std::string quartfile = create_quartfile(most_present_rule, quartet_counts, taxonNames, n_taxa);

    //write quartfile to outputFile

    std::ofstream ofs(outputFile);
    if (!ofs) {
        std::cerr << "ERROR: Failed to open output file." << std::endl;
        std::exit(1);
    }

    ofs << quartfile;
    return 0;
}