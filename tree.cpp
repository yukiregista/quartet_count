#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <thread>


class Node;

class Taxa {
private:
    int nextId = 1;
    std::vector<std::string> taxonNames;
    std::unordered_map<std::string, int> taxonIdMap;

public:
    const std::vector<std::string>& getTaxaNames() const { return taxonNames; }

    // Method to add a taxon
    void addTaxon(const std::string& name) {
        taxonNames.push_back(name);
        taxonIdMap[name] = nextId++;
    }

    int getTaxonId(const std::string& name) const {
        auto it = taxonIdMap.find(name);
        return (it != taxonIdMap.end()) ? it->second : -1;
    }
    int getNextId() const {return nextId;}

    std::vector<std::string> getTaxonNames() const {return taxonNames;}
};

#include <cmath> // for NaN


class Node {
public:
    int getTaxonId() const { return taxonId; }
    double getBranchLength() const { return branchLength; }

    void setParent(Node* par){parent = par;}
    Node* getParent(){return parent;}

    const std::vector<std::unique_ptr<Node>>& getChildren() const { return children; }

    // Add existing node as a child
    void addChild(std::unique_ptr<Node> child) {
        child->setParent(this);
        children.push_back(std::move(child));
    }

    // Create a new node with an associated taxon ID and branch length, and add as a child
    void addChild(int associatedTaxonId, double branchLength);

    // Create a new node with an associated taxon name, branch length, and add as a child
    void addChild(const std::string& associatedTaxonName, const Taxa& taxa, double branchLength);

    bool is_leaf() const {return children.size()==0;}

    void set_terminal_taxon_IDs(std::unordered_set<int> terminal){this->terminal_taxon_IDs = terminal;}
    void update_terminal_taxon_IDs(std::unordered_set<int> terminal){this->terminal_taxon_IDs.insert(terminal.cbegin(), terminal.cend());}
    std::unordered_set<int> getTerminal_taxon_IDs() const {return terminal_taxon_IDs;}
    std::vector<std::unique_ptr<Node>> children;

private:
    int taxonId;
    double branchLength; // Stores branch length or NaN if not specified
    std::unordered_set<int> terminal_taxon_IDs;
    Node* parent;

    // Constructor with associated Taxon ID and branch length (private)

    friend class Tree;        // Tree class can access the private constructor
    friend class NewickParser; // NewickParser class can access the private constructor

public:
    // Factory method to create a Node with associated taxon ID and branch length
    static std::unique_ptr<Node> createNode(int associatedTaxonId = 0, double branchLength = std::nan("")) {
        return std::unique_ptr<Node>(new Node(associatedTaxonId, branchLength));
    }
    Node(int associatedTaxonId, double branchLength) : taxonId(associatedTaxonId), branchLength(branchLength) {}
};

bool is_leaf(const Node& node){
    bool a=node.getChildren().size()==0;
    return a;
}

void Node::addChild(int associatedTaxonId, double branchLength) {
    std::unique_ptr<Node> newNode = Node::createNode(associatedTaxonId, branchLength);
    children.push_back(std::move(newNode));
    newNode->setParent(this);
}

void Node::addChild(const std::string& associatedTaxonName, const Taxa& taxa, double branchLength) {
    int associatedTaxonId = taxa.getTaxonId(associatedTaxonName);

    if (associatedTaxonId != -1) {
        std::unique_ptr<Node> newNode = Node::createNode(associatedTaxonId, branchLength);
        children.push_back(std::move(newNode));
        newNode->setParent(this);
    } else {
        // Handle error, print a message, or throw an exception
        std::cerr << "Error: Taxon not found for name " << associatedTaxonName << std::endl;
    }
}

class Edge {
public:
    const Node* getParent() const { return parent; }
    const Node* getChild() const { return child; }
    double getBranchLength() const { return branchLength; }

    Edge(const Node* parent, const Node* child, double branchLength)
        : parent(parent), child(child), branchLength(branchLength) {}

private:
    const Node* parent;
    const Node* child;
    double branchLength;
};




class Tree {
public:
    Node& getRoot() const { return *root; }
    const Taxa getTaxa() const {return taxa;}

    Tree(std::unique_ptr<Node> root, Taxa taxa) : root(std::move(root)), taxa(taxa) {
        postOrderTraversal_change([this](Node& node){
            if (node.is_leaf()){
                this->add_leafnodes(&node);
            }
        }, *this->root);
    }

    void add_leafnodes(Node* newNode){
        leaf_nodes.push_back(newNode);
    }
    
    std::vector<Node*> const getLeafNodes(){
        return leaf_nodes;
    }
    

    // Function to set the root of the tree
    void setRoot(std::unique_ptr<Node> newRoot) {
        root = std::move(newRoot);
    }

    // Function to print an ASCII plot of the phylogenetic tree
    void printAsciiTree(const Node& node, int level = 0) const {
        // Print the current node
        for (int i = 0; i < level; ++i) {
            std::cout << "| ";
        }
        std::cout << "+-" << "Taxon ID: " << node.getTaxonId() << "\n";

        // Recursively print the children
        for (const auto& child : node.getChildren()) {
            printAsciiTree(*child, level + 1);
        }
    }

    // Function to get a vector of edges in the tree
    std::vector<Edge> getEdges(const Node& node, const Taxa& taxa, double branchLength = 1.0) const {
        std::vector<Edge> edges;

        for (const auto& child : node.getChildren()) {
            const Node* parentNode = &node;
            const Node* childNode = child.get();
            edges.emplace_back(parentNode, childNode, branchLength);

            // Recursively get edges for child nodes
            auto childEdges = getEdges(*child, taxa, branchLength);
            edges.insert(edges.end(), childEdges.begin(), childEdges.end());
        }

        return edges;
    }

    // Function to print edges in the tree
    void printEdges(const Taxa& taxa) const {
        auto edges = getEdges(*root, taxa);

        std::cout << "Edges in the tree:\n";
        for (const auto& edge : edges) {
            std::cout << "Parent Taxon ID: " << edge.getParent()->getTaxonId()
                      << ", Child Taxon ID: " << edge.getChild()->getTaxonId()
                      << ", Branch Length: " << edge.getBranchLength() << "\n";
        }
    }

    // // Function to perform DFS on the tree and apply a function to each node
    // template <typename Function>
    // void depthFirstSearch(Function&& func, const Node& node) const {
    //     func(node);  // Apply the function to the current node

    //     for (const auto& child : node.getChildren()) {
    //         depthFirstSearch(std::forward<Function>(func), *child);  // Recursively apply to child nodes
    //     }
    // }

    // Function to perform post-order traversal on the tree and apply a function to each node
    template <typename Function>
    void postOrderTraversal_change(Function&& func, Node& node) const {
        for (const auto& child : node.getChildren()) {
            postOrderTraversal_change(std::forward<Function>(func), *child);  // Recursively apply to child nodes
        }

        func(node);  // Apply the function to the current node
    }
        // Function to perform post-order traversal on the tree and apply a function to each node
    template <typename Function>
    void postOrderTraversal(Function&& func, const Node& node) const {
        for (const auto& child : node.getChildren()) {
            postOrderTraversal(std::forward<Function>(func), *child);  // Recursively apply to child nodes
        }

        func(node);  // Apply the function to the current node
    }

private:
    std::unique_ptr<Node> root;
    Taxa taxa;
    std::vector<Node*> leaf_nodes;
};


class NewickParser {
public:
    Tree parseNewick(const std::string& newickString);

private:
    std::unique_ptr<Node> parseInternal(std::string::const_iterator& it, const std::string& newickString, Taxa & taxa);

    void skipWhitespace(std::string::const_iterator& it, const std::string& newickString);

    std::string readUntil(std::string::const_iterator& it, const std::unordered_set<char>& delimiters, const std::string& newickString);

    double readBranchLength(std::string::const_iterator& it, const std::string& newickString);
};


Tree NewickParser::parseNewick(const std::string& newickString) {
    std::string::const_iterator it = newickString.begin();
    Taxa taxa;
    std::unique_ptr<Node> a = parseInternal(it, newickString, taxa);
    a->setParent(nullptr);
    Tree tree(std::move(a), taxa);
    return tree;
}

std::unique_ptr<Node> NewickParser::parseInternal(std::string::const_iterator& it, const std::string& newickString, Taxa & taxa) {
    if (*it == '(') {
        // Internal node
        ++it; // Skip '('

        std::unique_ptr<Node> internal1 = Node::createNode();
        while (*it != ')') {
            std::unique_ptr<Node> tmp = std::move(parseInternal(it, newickString, taxa));
            internal1->addChild(std::move(tmp));
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
        return internal1; // Internal node has no associated taxon
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
        taxa.addTaxon(taxonName);
        return std::make_unique<Node>(taxa.getNextId()-1, branchLength);
    }
}


void NewickParser::skipWhitespace(std::string::const_iterator& it, const std::string& newickString) {
    while (it != std::end(newickString) && std::isspace(*it)) {
        ++it;
    }
}

std::string NewickParser::readUntil(std::string::const_iterator& it, const std::unordered_set<char>& delimiters, const std::string& newickString) {
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

void printTaxonId(const Node& node) {
    std::cout << "Taxon ID: " << node.getTaxonId() << std::endl;
}

void set_terminal_taxon_IDs(Node& node) {
    if (is_leaf(node)){
        node.set_terminal_taxon_IDs(std::unordered_set<int>{node.getTaxonId()});
    } else{
        for (const std::unique_ptr<Node>& child : node.getChildren()){
            std::unordered_set<int> tmp = child->getTerminal_taxon_IDs();
            node.update_terminal_taxon_IDs(tmp);
        }
    }
}

void printTerminalId(const Node& node) {
    std::cout << "Taxon ID: " << node.getTaxonId() << std::endl;
    for (int id : node.getTerminal_taxon_IDs()){
        std::cout << "Terminal ID: " << id << std::endl;
    }
    std::cout << std::endl;
}


void reroot_tree_at_taxon_id(Tree & tree, int taxon_id){
    // reroot tree with taxon 
    const Taxa& taxa = tree.getTaxa();
    if (taxa.getNextId() <= taxon_id){
        // taxon does not exist. Do nothing.
        return;
    }
    
}

void reroot_tree_at_leaf_index(Tree & tree, int leaf_ind){
    // reroot tree with leaf
    // for rerooting one at a time.
     const Taxa& taxa = tree.getTaxa();
     Node* leaf = tree.getLeafNodes()[leaf_ind];
     if (&tree.getRoot()==leaf){
        // do nothing
        return;
     }
     Node* chi = leaf;
     while (chi!=&tree.getRoot()){
        Node* par = chi->getParent();
        chi->addChild(std::move(*par));
        auto it = std::remove_if(par->children.begin(), par->children.end(),
                             [chi](const std::unique_ptr<Node>& ptr) {
                                 return ptr.get() == chi;
                             });
        par->children.erase(it, par->children.end()); // reduce physical size
     }
}


int main() {
    // Example usage
    // Taxa taxa;
    // taxa.addTaxon("Species1");
    // std::unique_ptr<Node> leaf1 = Node::createNode(taxa.getTaxonId("Species1"));
    // taxa.addTaxon("Species2");
    // std::unique_ptr<Node> leaf2 = Node::createNode(taxa.getTaxonId("Species2"));
    // taxa.addTaxon("Species3");
    // std::unique_ptr<Node> leaf3 = Node::createNode(taxa.getTaxonId("Species3"));
    // taxa.addTaxon("Species4");
    // std::unique_ptr<Node> leaf4 = Node::createNode(taxa.getTaxonId("Species4"));
    // taxa.addTaxon("Species5");
    // std::unique_ptr<Node> leaf5 = Node::createNode(taxa.getTaxonId("Species5"));

    
    // std::unique_ptr<Node> internal1 = Node::createNode();
    // std::unique_ptr<Node> internal2 = Node::createNode();
    

    // internal1->addChild(std::move(leaf1));
    // internal1->addChild(std::move(leaf2));
    // internal2->addChild(std::move(leaf3));
    // internal2->addChild(std::move(leaf4));
    // internal2->addChild(std::move(leaf5));

    // std::unique_ptr<Node> rootNode = Node::createNode();
    // rootNode->addChild(std::move(internal1));
    // rootNode->addChild(std::move(internal2));


    // Tree phylogeneticTree(std::move(rootNode));

    // // Print the ASCII tree
    // phylogeneticTree.printAsciiTree(phylogeneticTree.getRoot());

    // // Print edges in the tree
    // phylogeneticTree.printEdges(taxa);


    //std::string newickString = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5)E:0.6;";
    std::string newickString = "(A,B,(C,D), (G,(E,F)));";
    NewickParser newickParser;
    Tree phylogeneticTree = newickParser.parseNewick(newickString);
    std::vector<std::string> a = phylogeneticTree.getTaxa().getTaxonNames();
    for (std::string b : a){
        std::cout << b << std::endl;
    }

    // Now you can use the parsed tree as needed, for example, print the ASCII tree
    phylogeneticTree.printAsciiTree(phylogeneticTree.getRoot());

    reroot_tree_at_leaf_index(phylogeneticTree, 0);
    phylogeneticTree.printAsciiTree(phylogeneticTree.getRoot());

    phylogeneticTree.postOrderTraversal(printTaxonId, phylogeneticTree.getRoot());
    phylogeneticTree.postOrderTraversal(is_leaf, phylogeneticTree.getRoot());

    phylogeneticTree.postOrderTraversal_change(set_terminal_taxon_IDs, phylogeneticTree.getRoot());
    phylogeneticTree.postOrderTraversal(printTerminalId, phylogeneticTree.getRoot());

    std::vector<Node*> leaves = phylogeneticTree.getLeafNodes();
    for (Node* leaf : leaves){
        std::cout << leaf->getTaxonId() << std::endl;
    }

    return 0;
}
