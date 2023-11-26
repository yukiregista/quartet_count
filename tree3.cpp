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
    std::vector<std::string> taxonNames;
    std::unordered_map<std::string, int> taxonIdMap;

public:
    const std::vector<std::string>& getTaxaNames() const { return taxonNames; }

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

    std::vector<std::string> getTaxonNames() const {return taxonNames;}
};


class Node {
public:
    int getNodeId() const {return NodeId;}
    int getTaxonId() const { return taxonId; }
    double getBranchLength() const { return branchLength; }

    void setParent(int par){parent = par;}
    int getParent(){return parent;}

    const std::set<int>& getChildren() const { 
        return children; 
        }

    // Add existing node as a child
    void addChild(int childID) {
        tree->nodes[childID]->setParent(this->NodeId);
        children.insert(childID);
    }

    bool is_leaf() const {return children.size()==0;}

    void set_terminal_taxon_IDs(std::set<int> terminal){this->terminal_taxon_IDs = terminal;}
    void update_terminal_taxon_IDs(std::set<int> terminal){this->terminal_taxon_IDs.insert(terminal.cbegin(), terminal.cend());}
    void clear_terminal_taxon_IDs(){this->terminal_taxon_IDs.clear();}
    std::set<int> getTerminal_taxon_IDs() const {return terminal_taxon_IDs;}
    Tree* getTree() const {return tree;}
    

private:
    std::set<int> children{};
    int taxonId;
    int NodeId;
    double branchLength; // Stores branch length or NaN if not specified
    std::set<int> terminal_taxon_IDs{};
    int parent;
    Tree* tree;

    // Constructor with associated Taxon ID and branch length (private)

    friend class Tree;        // Tree class can access the private constructor
    friend class NewickParser; // NewickParser class can access the private constructor

public:
    // Factory method to create a Node with associated taxon ID and branch length and tree
    static std::unique_ptr<Node> createNode(Tree* tree, int associatedTaxonId = 0, double branchLength = std::nan("")) {
        return std::unique_ptr<Node>(new Node(associatedTaxonId, branchLength, tree));
    }
    Node(int associatedTaxonId, double branchLength, Tree* associatedTree) : taxonId(associatedTaxonId), branchLength(branchLength), NodeId(associatedTaxonId), tree(associatedTree) {}
};

void set_terminal_taxon_IDs(Node& node) {
    if (node.is_leaf()){
        node.set_terminal_taxon_IDs(std::set<int>{node.getTaxonId()});
    } else{
        for (int child : node.getChildren()){
            std::set<int> tmp = node.getTree()->nodes[child]->getTerminal_taxon_IDs();
            node.update_terminal_taxon_IDs(tmp);
        }
    }
}

void reset_terminal_taxon_IDs(Node& node) {
    node.clear_terminal_taxon_IDs();
    if (node.is_leaf()){
        node.set_terminal_taxon_IDs(std::set<int>{node.getTaxonId()});
    } else{
        for (int child : node.getChildren()){
            std::set<int> tmp = node.getTree()->nodes[child]->getTerminal_taxon_IDs();
            //for (auto item : tmp){ std::cout << " " << item << " ";}
            node.update_terminal_taxon_IDs(tmp);
            //for (auto item : node.getTerminal_taxon_IDs()){std::cout << item << std::endl;};
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

    Node* root;
    Taxa* taxa;
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
        Node* nextpar = par->parent;
        std::cout << chi->NodeId << std::endl;
        chi->addChild(par); // add par as child of par
        par->setParent(chi); // set chi as parent of par
        par->children.erase(chi->NodeId); //delete chi from children of par
        chi = par;
        par = nextpar;
    }
    leaf->setParent(nullptr);
    root = leaf.get();
//     std::this_thread::sleep_for(std::chrono::milliseconds(1000));
//     for (auto child :nodes[-1]->children){
//         std::cout << child.second->NodeId << std::endl;
//     }
//     std::this_thread::sleep_for(std::chrono::milliseconds(1000));

}

void count_quartets_of_the_node(const Node& node, std::vector<int>& quartet_counts, const std::set<int>& all_ids,const int& root_id,const int& n_taxa){
    // count the quartets associated with this node
   std::set<int> children = node.getChildren();
   Tree* tree = node.getTree();
   int children_size = children.size();
   if (children_size <=1){
        // either leaf or unification. no need for updateing
        return ;
   }
    std::set<int> terminals = node.getTerminal_taxon_IDs();
    std::set<int> others; // will contain other terminals (root is not considered as terminals)
    std::set_difference(all_ids.begin(), all_ids.end(), terminals.begin(), terminals.end(), 
                        std::inserter(others, others.begin()));
    if (others.size()==0){return;}
    //for (auto item : others){std::cout << item << " ";}
    for (auto child1 = children.cbegin(); child1 != std::prev(children.cend()); child1++) {
        //access all terminal of this child.
        for (int term1 : tree->nodes[*child1]->getTerminal_taxon_IDs()){
            // term 1 is picked from child 1
            for (auto child2 = std::next(child1); child2 != children.cend(); child2++) {
                for (int term2 : tree->nodes[*child2]->getTerminal_taxon_IDs()){
                    // term2 is picked from child 2
                    for (int other_term : others){
                        // (term1, term2) | (root_id, other_term) is present.
                        // First, make an array of these four elements.
                        std::array<int,4> quartet = {term1, term2, other_term, root_id};
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
    int count = 0;
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
        std::unique_ptr<Node> internal1 = Node::createNode();
        internal1->NodeId = tree->taxa->getNextInternalId() + 1;
        while (*it != ')') {
            Node* tmp = parseInternal(it, newickString, tree);
            internal1->addChild(tmp->getNodeId());
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
            ;
        } else {
            // internal node is labeled -- but we don't store that information.
            std::string taxonName = readUntil(it, { ':', ',', ')' }, newickString);
            if (*it == ':'){
                ++it;
                branchLength = readBranchLength(it, newickString);
            } else if (*it == ','){
                ++it;
                branchLength = std::nan("");
            }
        }
        skipWhitespace(it, newickString);
        int nodeId = internal1->NodeId;
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
        tree->nodes[taxonId] = std::make_unique<Node>(taxonId, branchLength, tree.get());
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
            std::cout << "ファイルが開けませんでした。" << std::endl;
            std::cin.get();
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
    

    return std::move(phylogeneticTree);
}


void printTerminalId(const Node& node) {
    std::cout << "Taxon ID: " << node.getTaxonId() << std::endl;
    for (int id : node.getTerminal_taxon_IDs()){
        std::cout << "Terminal ID: " << id << std::endl;
    }
    std::cout << std::endl;
}


void count_quartets(std::vector<std::unique_ptr<Tree>>& trees){
    // count the number of quartets
    int n_taxa = trees[0]->taxa->getNextId()-1;
    int nC4 = BinomialCoefficient(n_taxa,4);
    std::vector<int> quartet_counts(nC4*3, 0);
    for (std::unique_ptr<Tree> & tree : trees){
        // count quartets in the tree.

        // Iteration over leaf nodes
        for (int reroot_index=1; reroot_index <= n_taxa; reroot_index++){

            // 1. Rerooting the tree with the leaf
            //tree->reroot_with_leaf_ind(reroot_index); //rerooting; changing parent-child structure
            // tree->postOrderTraversal(printTerminalId, *tree->root);
            tree->postOrderTraversal_change(reset_terminal_taxon_IDs, *tree->root); // reset terminals
            tree->nodes[-3]->getTerminal_taxon_IDs();
            // tree->postOrderTraversal(printTerminalId, *tree->root);
            // postOrderTraversal to count all quartets inccluding reroot_index
            //tree->postOrderTraversal(printTerminalId, *tree->root);
            std::set<int> all_ids;
            for (int child: tree->root->getChildren()){
                all_ids = tree->nodes[child]->getTerminal_taxon_IDs();
                break;
            }
            tree->postOrderQuartetCounting(*tree->root, quartet_counts,all_ids, reroot_index, n_taxa);
            //tree->postOrderTraversal(count_quartets_of_the_node, *tree->root, quartet_counts,all_ids, reroot_index, n_taxa);
        }   
    }
}






int main(){
    std::unique_ptr<Taxa> taxa = std::make_unique<Taxa>();
    NewickParser newickParser;
    std::vector<std::unique_ptr<Tree>> phylogeneticTree = newickParser.parseNewick("((A,B),(C,D),E);\n ((A,C),(B,D),(E,F))", taxa); //ok
    for (int i=0; i<phylogeneticTree.size(); i++){
        phylogeneticTree[i]->printAsciiTree(*(phylogeneticTree[i]->root), 0);
        std::cout << std::endl;
    }
    std::cout << phylogeneticTree[0]->nodes[-3]->getParent() << std::endl;

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


    //file IO
    // std::unique_ptr<Taxa> taxa = std::make_unique<Taxa>();
    // std::vector<std::unique_ptr<Tree>> phylogeneticTree =readNewickFile("true.tre", taxa);
    // phylogeneticTree[0]->printAsciiTree(*phylogeneticTree[0]->root);
    // //phylogeneticTree[0]->reroot_with_leaf_ind(1);
    // phylogeneticTree[0]->postOrderTraversal_change(reset_terminal_taxon_IDs, *phylogeneticTree[0]->root);
    // phylogeneticTree[0]->postOrderTraversal(printTerminalId, *phylogeneticTree[0]->root);//ok
    count_quartets(phylogeneticTree);
    return 0;
}