////
//// Created by fritsche on 17/05/2021.
////
//
#include "Taxonomy.h"
#include "Utils.h"
#include <iostream>
#include <fstream>
#include <unordered_set>


void protal::taxonomy::StdTaxonomy::LoadNodes(std::string path) {
    std::ifstream is(path.c_str(), std::ios::in);

    bool header = false;
    std::string line;
    std::vector<std::string> tokens;
    while (std::getline(is, line)) {
        Utils::split(tokens, line, "\t");

        if (header) {
            header = false;
            continue;
        }


        int id = stoi(tokens[0]);
        int parent_id = stoi(tokens[2]);
        std::string rank = tokens[4];

        Node* node = GetNode(id);
        node->rank = rank;

        Node* pnode = nullptr;

        // Root is defined in the file that node id and parent node id are identical
        if (id == parent_id) {
            root_id_ = id;
        } else {
            pnode = GetNode(parent_id);
        }

        // Set parent, is nullptr if root
        node->parent_node = pnode;
        if (pnode) {
            pnode->children.emplace_back(node);
        }
    }
    Level();
}

protal::taxonomy::Node* protal::taxonomy::StdTaxonomy::GetNode(int taxid) {
    Node *node = nullptr;
    if (node_from_id_.find(taxid) == node_from_id_.end()) {
        node = new Node();
        node->id = taxid;
        node_from_id_.insert({ taxid, node });
        nodes_.emplace_back(node);
    } else {
        node = node_from_id_[taxid];
    }
    return node;
}

void protal::taxonomy::StdTaxonomy::DeleteNode(Node* node) {
    // Safety belt!
    // Node cannot be in the node_from_id_ map if it is gonna be deleted from memory
    for (NodeFromId ::const_iterator it = node_from_id_.begin(); it != node_from_id_.end(); ++it) {
        if (it->second == node) {
            std::cerr << "Trying to delete node that is still in the map. " << node->id << std::endl;
            exit(13);
        }
    }

    nodes_.erase(std::remove(nodes_.begin(), nodes_.end(), node), nodes_.end());
    delete node;
}

protal::taxonomy::StdTaxonomy::StdTaxonomy() {
    InitRankMaps();
}

protal::taxonomy::StdTaxonomy::~StdTaxonomy() {
    while (!nodes_.empty()) {
        Node* tmp = nodes_.back();
        nodes_.pop_back();
        delete tmp;
    }
}

std::vector<protal::taxonomy::Node *> &protal::taxonomy::StdTaxonomy::Nodes() {
    return nodes_;
}


void protal::taxonomy::StdTaxonomy::ClearNodes() {
//    std::vector<Node*> tmp;
//    for (auto e : nodes_) {
//        if (!node_from_id_.contains(e->id)) {
//            tmp.insert(e);
//        }
//    }
}

void protal::taxonomy::StdTaxonomy::Level() {
    int level = 0;
    std::cout << "rootid: " << root_id_ << std::endl;
    std::cout << "level: " << node_from_id_[root_id_]->id << std::endl;
    Level(node_from_id_[root_id_], level);
}

void protal::taxonomy::StdTaxonomy::Level(Node* node, int level) {
    node->level = level++;
    for (auto child : node->children) {
        Level(child, level);
    }
}

void protal::taxonomy::StdTaxonomy::InitRankMaps() {
    rank_to_id_.insert({ "no rank", 0 });
    rank_to_id_.insert({ "domain", 1 });
    rank_to_id_.insert({ "phylum", 2 });
    rank_to_id_.insert({ "class", 3 });
    rank_to_id_.insert({ "order", 4 });
    rank_to_id_.insert({ "family", 5 });
    rank_to_id_.insert({ "genus", 6 });
    rank_to_id_.insert({ "species", 7 });
    rank_to_id_.insert({ "subspecies", 8 });

    for (auto pair : rank_to_id_) {
        rank_from_id_.insert({ pair.second, pair.first });
    }
}

void protal::taxonomy::StdTaxonomy::LoadNames(std::string path) {

    std::ifstream is(path.c_str(), std::ios::in);

    bool header = false;
    std::string line;
    std::vector<std::string> tokens;
    while (std::getline(is, line)) {
        Utils::split(tokens, line, "\t");

        if (header) {
            header = false;
            continue;
        }

        int id = stoi(tokens[0]);


        if (node_from_id_.find(id) == node_from_id_.end()) {
            std::cerr << "Warning: there is no node with id: " << id << std::endl;
            continue;
        }

        if (tokens.size() < 7) {
            std::cerr << "corrupt line in names" << std::endl;
            std::cerr << line << std::endl;
            exit(2);
        }
        std::string name = tokens[2];
        std::string name_unique = tokens[4];
        std::string name_class = tokens[6];


        auto node = node_from_id_[id];

        if (name_class.compare("scientific name") == 0 || !node->name) {
            if (!node->name) {
                node->name = new Name();
            }

            auto name_obj = node->name;
            name_obj->name = name;
            name_obj->name_class = name_class;
            name_obj->name_unique = name_unique;
            name_obj->id = id;
        }

    }
}


void protal::taxonomy::StdTaxonomy::Export(std::string path) {
    std::cout << "Export to " << path << std::endl;
    std::ofstream ofs(path);
    if (!ofs)
        std::cerr << "unable to open outfile " << path << std::endl;

    std::string main_name = "scientific name";
    // The old id (e.g. ncbi) is stored in name
    std::string oid_name = "old id";

    ofs << "id\tparent_id\texternal_id\tname\trank\tlevel" << std::endl;

    for (auto node : nodes_) {
        if (node->external_id == -1) exit(8);
        auto name = node->name;

        ofs << node->id << separator;
        ofs << (node->parent_node ? node->parent_node->id : root_id_) << separator;
        ofs << node->external_id << separator;
        ofs << name->name << separator;
        ofs << node->rank << separator;
        ofs << node->level << separator;
        ofs << node->rep_genome << std::endl;
    }

    ofs.close();
}

void protal::taxonomy::StdTaxonomy::SubsetAndRelabel(std::string path, std::string collapse_level) {
    std::cout << "Delete intermediate nodes" << std::endl;

    bool collapse_if_only_child = true;

    std::cout << "#nodes: " << nodes_.size() << std::endl;
    std::cout << "sanity check" << std::endl;


    auto r = GetNode(1);
    std::cout << "before root1 furst: " << r->id << " " << r->children.size() << std::endl;
    std::for_each(r->children.begin(), r->children.begin() + std::min(20lu, r->children.size()), [](Node* const node) {
        std::cout << node->id << ",";
    });
    std::cout << std::endl;

    std::cout << "RS_GCF_004361695.1 (51724) exists? " << nodes_.at(51724)->name->name << std::endl;

    SanityCheck();

    KeepRanksAndLeaves(nullptr, GetNode(root_id_));

    std::cout << "RS_GCF_004361695.1 (51724) exists? " << nodes_.at(51724)->name->name << std::endl;

    std::cout << "#nodes: " << nodes_.size() << std::endl;
    std::cout << "before root1a second: " << r->id << " " << r->children.size() << std::endl;

    std::cout << "before root1 furst: " << r->id << " " << r->children.size() << std::endl;
    std::for_each(r->children.begin(), r->children.begin() + std::min(20lu, r->children.size()), [](Node* const node) {
        std::cout << node->id << ",";
    });
    std::cout << std::endl;

    // Subset according to list with species int pairs
    std::cout << "subset file " << path << std::endl;
    if (!Utils::exists(path)) {
        std::cerr << path << " does not exist." << std::endl;
        exit(8);
    }

    std::ifstream is(path.c_str(), std::ios::in);

//    auto str = csv::load_file(path.c_str());
//    auto parser = csv::make_parser( str , '\t');

    bool header = false;

    tsl::sparse_map<std::string, TaxId> species2id;
    tsl::sparse_set<int> blocked_ids;

    std::cout << "RS_GCF_004361695.1 (51724) exists? " << (species2id.find("RS_GCF_004361695.1") != species2id.end()) << std::endl;

    std::cout << "before roota: " << GetNode(1)->id << " " << GetNode(1)->children.size() << std::endl;

    std::string line;
    std::vector<std::string> tokens;
    while (std::getline(is, line)) {
        Utils::split(tokens, line, "\t");
//    for (auto&& row : parser ) {
        if (header) {
            header = false;
            continue;
        }

        std::string species = tokens.at(0);
        int id = stoi(tokens.at(1));

        species2id.insert({ species, id });
        blocked_ids.insert(id);
    }
    std::cout << "before rootb: " << GetNode(1)->id << " " << GetNode(1)->children.size() << std::endl;

    std::cout << "species2id: " << species2id.size() << std::endl;

    std::cout << "before root1 furst: " << r->id << " " << r->children.size() << std::endl;
//    for (auto& pair : species2id) {
//        std::cout << pair.first << " " << pair.second << std::endl;
//    }

    std::cout << "RS_GCF_004361695.1 (51724) exists? " << nodes_.at(51724)->name->name << std::endl;

    std::cout << " len : " << species2id.size() << std::endl;

    tsl::sparse_set<Node*> nodes_to_keep;

    node_from_id_.clear();

    std::cout << "collect nodes to keep.. (" << nodes_.size() << ")" << std::endl;



    for (auto node : nodes_) {
        std::string scientific_name = node->name->name;

        if (node->IsLeaf())
            std::cout << node->name->name << " -> " << (species2id.find(scientific_name) != species2id.end()) << std::endl;
        // If node is in nodes to keep
        if (species2id.find(scientific_name) != species2id.end()) {

            Node* parent = node;

            // Node is flagged
            nodes_to_keep.insert(node);

            if (node->children.size() == 1) {
                auto& child = node->children[0];
                if (child->name) {
                    node->rep_genome = child->name->name;
//                    std::cout << "name: " << child->name->name << std::endl;
                }
            }


            while (!parent->IsRoot()) {
                parent = parent->parent_node;
                nodes_to_keep.insert(parent);
            }
        } else {
            //std::cout << scientific_name << std::endl;
        }
    }
    std::cout << "before root  (9): " << r->id << " " << r->children.size() << std::endl;

    std::cout << "species2id: " << species2id.size() << std::endl;
    std::cout << "keep: " << nodes_to_keep.size() << std::endl;

    if (nodes_to_keep.empty()) {
        std::cout << "no nodes to keep..." << std::endl;
        exit(9);
    }

    std::cout << "Old size " << nodes_.size() << std::endl;


    // New ids (start with 1)
    int running_id = 1;


    std::vector<Node*> nodes2remove;
//    auto itr =  nodes_.begin();
    size_t counter = 0;
//    while (itr != nodes_.end())

    std::cout << "before root  (10): " << r->id << " " << r->children.size() << std::endl;
    for (auto node : r->children) {
        if (node->rank == "phylum") {
            std::cout << node->id << " " << node->name->name << std::endl;
        }
    }
//    exit(9);
    node_from_id_.clear();
    for (auto& node : nodes_)
    {
        if (!nodes_to_keep.contains(node)) {
            nodes2remove.push_back(node);
        } else {
            // We have a new id for this node
            if (species2id.find(node->name->name) != species2id.end()) {
//                std::cout << "known: " << node->name->name << " " << species2id[node->name->name] << std::endl;
                node->ChangeId(species2id[node->name->name]);
            } else {
                while (blocked_ids.contains(running_id))
                    running_id++;

                if (!node->parent_node) root_id_ = running_id;

                if (collapse_if_only_child && node->children.size() == 1 && species2id.find(node->children[0]->name->name) != species2id.end()) {
                    node->ChangeId(species2id[node->children[0]->name->name]);
                    node->rep_genome = node->children[0]->name->name;
                    nodes2remove.push_back(node->children[0]);
//                    std::cout << "remove_leaf: " << node->external_id << " new: " << node->id << std::endl;
                    continue;
                }

                node->ChangeId(running_id);
//                std::cout << "old: " << node->external_id << " new: " << node->id << std::endl;

                running_id++;
            }
            node_from_id_[node->id] = node;

//            ++itr;
        }
    }

    std::cout << "keep: " << nodes_to_keep.size() << std::endl;
    std::cout << "Remove and unlink " << nodes2remove.size() << "/" << nodes_.size() << std::endl;

    for (auto n : nodes2remove)
        n->RemoveRelink();


    nodes_.clear();
    std::cout << "RepopulateNodes" << std::endl;
    RepopulateNodes();
    node_from_id_.clear();
    std::cout << "Reinsert in map" << std::endl;
    for (auto n : nodes_)
        node_from_id_.insert( { n->id, n } );
    std::cout << "Delete unused nodes" << std::endl;
    for (auto n : nodes2remove)
        delete n;

    std::cout << "New size " << nodes_.size() << std::endl;
}

void protal::taxonomy::StdTaxonomy::SanityCheck(Node* node) {
    if (node == nullptr)
        node = GetNode(1);

    for (auto child : node->children) {
        if (!node_from_id_.contains(child->id)) {
            std::cerr << "failed for " << node->id << " and child " << child->id << std::endl;
            exit(23);
        }

        SanityCheck(child);
    }
}

void protal::taxonomy::StdTaxonomy::KeepRanksAndLeaves(tsl::sparse_set<std::string> &ranks, std::vector<Node*> &remove_list, Node* parent_node, Node* node) {
//    std::cout << "KeepRanksAndLeaves " << " " << node->id << std::endl;
    bool init = false;
    if (parent_node == nullptr) {
        init = true;
        parent_node = node;
    } else if ( (node->IsLeaf() || ranks.contains(node->rank)) ) {
        parent_node = node;
    } else {
        if (node->id == 1) {
            std::cerr << node->id << std::endl;
            exit(17);
        }
        remove_list.push_back(node);
    }
//    std::cout << "iterate children (" << node->children.size() << ") of " << node->id << std::endl;
    for (auto child : node->children) {
        if (!node_from_id_.contains(child->id)) {
            std::cerr << "Failure child: " << child->id << " " << node->id << " #children: " << node->children.size() << std::endl;
            exit(9);
        }
        KeepRanksAndLeaves(ranks, remove_list, parent_node, child);
    }
    if (init) {
        std::cout << "Remove and unlink" << std::endl;
        for (auto n : remove_list) {
//            node_from_id_.erase(n->id);
            n->RemoveRelink();
//            DeleteNode(n);
        }
        nodes_.clear();
        std::cout << "RepopulateNodes" << std::endl;
        RepopulateNodes();
        node_from_id_.clear();
        std::cout << "Reinsert in map" << std::endl;
        for (auto n : nodes_)
            node_from_id_.insert( { n->id, n } );
        std::cout << "Delete unused nodes" << std::endl;
        for (auto n : remove_list)
            delete n;

    }
}

void protal::taxonomy::StdTaxonomy::KeepRanksAndLeaves(Node *parent_node, protal::taxonomy::Node *node) {
    tsl::sparse_set<std::string> keep_ranks({
        "phylum", "class", "order", "family", "genus", "species"
    });
    std::vector<Node*> remove_list;
    KeepRanksAndLeaves(keep_ranks, remove_list, parent_node, node);
}


void protal::taxonomy::StdTaxonomy::RepopulateNodes(Node* node) {
    if (node == nullptr)
        node = GetNode(root_id_);

    nodes_.push_back(node);

    for (auto child : node->children)
        RepopulateNodes(child);
}

void protal::taxonomy::StdTaxonomy::Add(std::string path) {
    std::ifstream is(path.c_str(), std::ios::in);

    bool header = false;
    std::string line;
    std::vector<std::string> tokens;

    while (std::getline(is, line)) {
        Utils::split(tokens, line, "\t");
        if (header) {
            header = false;
            continue;
        }

        std::string accession = tokens.at(0);
        int external_id = stoi(tokens.at(1));
    }
}



void protal::taxonomy::IntTaxonomy::Load(std::string path) {
    std::ifstream is(path.c_str(), std::ios::in);

    bool header = true;
    std::string line;
    std::vector<std::string> tokens;

    while (std::getline(is, line)) {
        Utils::split(tokens, line, "\t");

        if (header) {
            header = false;
            continue;
        }

        int id = stoi(tokens.at(0));

        int parent_id = stoi(tokens.at(1));
        int external_id = stoi(tokens.at(2));
        std::string name = tokens.at(3);
        std::string rank = tokens.at(4);
        int level = stoi(tokens.at(5));
        std::string rep_genome = tokens.at(6);

        if (map.find(id) == map.end())
            map.insert( { id, IntNode() } );

        if (map.find(parent_id) == map.end())
            map.insert( { parent_id, IntNode() } );

        if (id == parent_id) {
            root_id = id;
        }

        map.at(id).id = id;
        map.at(id).parent_id = parent_id;
        map.at(id).external_id = external_id;
        map.at(id).scientific_name = name;
        map.at(id).rank = rank;
        map.at(id).level = level;
        map.at(id).rep_genome = rep_genome;

        if (rank.empty()) {
            std::cout << "?: " << id << " " << parent_id << " " << name << " " << rank << std::endl;
        }
        map.at(parent_id).children.emplace_back(id);
    }
}

protal::taxonomy::IntTaxonomy::IntTaxonomy(std::string path) {
    Load(path);
}

int protal::taxonomy::IntTaxonomy::LCA(int t1, int t2) {
    if (map.find(t1) == map.end() || map.find(t2) == map.end()) {
        std::cout << "t1: " << t1 << std::endl;
        std::cout << "t2: " << t2 << std::endl;
    }
    assert(map.find(t1) != map.end() && map.find(t2) != map.end());
    int l1 = map.at(t1).level;
    int l2 = map.at(t2).level;


//    std::cout << t1 << " " << l1 << "   " << t2 << " " << l2 << std::endl;
    while (l1 > l2) {
//        auto n = map[map[t1].parent_id];
        auto n = map.at(map.at(t1).parent_id);
        l1 = n.level;
        t1 = n.id;
    }
    while (l2 > l1) {
//        auto n = map[map[t2].parent_id];
        auto n = map.at(map.at(t2).parent_id);
        l2 = n.level;
        t2 = n.id;
    }


//    std::cout << "lca3" << std::endl;
    auto iter = 0;
    while (t1 != t2) {
//        auto n1 = map[map[t1].parent_id];
//        auto n2 = map[map[t2].parent_id];
        auto n1 = map.at(map.at(t1).parent_id);
        auto n2 = map.at(map.at(t2).parent_id);

        l1 = n1.level;
        l2 = n2.level;
        t1 = n1.id;
        t2 = n2.id;


    }

    return t1;
}

protal::taxonomy::IntNode &protal::taxonomy::IntTaxonomy::Get(int t) {
    return map.at(t);
}

protal::taxonomy::IntNode &protal::taxonomy::IntTaxonomy::GetParent(int t) {
    return map.at(map.at(t).parent_id);
}

size_t protal::taxonomy::IntTaxonomy::GetSubtreeSize(int t) {
    auto& node = Get(t);
    size_t num = 0;
    if (node.children.size() == 0) {
        return 1;
    } else {
        for (auto& child_id : node.children) {
            num += GetSubtreeSize(child_id);
        }
    }
    return num;
}

size_t protal::taxonomy::IntTaxonomy::MaxTaxid() {
    return std::max_element(map.begin(), map.end(), [] (auto e1, auto e2) { return e1.first < e2.first; })->first;
}
size_t protal::taxonomy::IntTaxonomy::MaxLeafId() {
    return std::max_element(map.begin(), map.end(), [] (auto e1, auto e2) {
        auto e1_id = e1.second.children.empty() ? e1.first : 0;
        auto e2_id = e2.second.children.empty() ? e2.first : 0;
        return e1_id < e2_id; })->first;
}

bool protal::taxonomy::IntTaxonomy::IsNodeAncestor(int node_id, int leaf_id) {
    auto* leaf = &map.at(leaf_id);
    auto* node = &map.at(node_id);

    while (leaf->level > node->level) {
        leaf = &GetParent(leaf->id);
    }

    return node->id == leaf->id;
}

std::string protal::taxonomy::IntTaxonomy::LineageStr(int t, const std::vector<std::string> ranks, std::string divider)  {
    auto node = map.at(t);

    std::vector<std::string> result(ranks.size());

    if (node.IsLeaf()) {
        auto it = std::find(ranks.begin(), ranks.end(), "genome");
        if (it != ranks.end()) {
            int index = it - ranks.begin();
            result[index] = node.scientific_name;
        }
    }

    int nid = node.id;
    int pid = node.parent_id;

    while (nid != pid) {
        auto node = map.at(nid);

        auto it = std::find(ranks.begin(), ranks.end(), node.rank);
        if (it != ranks.end()) {
            int index = it - ranks.begin();
            result[index] = node.scientific_name;
        }
        nid = pid;
        pid = map.at(nid).parent_id;
    }

//    return lineage;
    return Utils::join(result, "|");
}

std::string protal::taxonomy::IntTaxonomy::LineageExternalIds(int t, const std::vector<std::string> ranks, std::string divider)  {
    auto node = map.at(t);

    std::vector<std::string> result(ranks.size());

    if (node.IsLeaf()) {
        auto it = std::find(ranks.begin(), ranks.end(), "genome");
        if (it != ranks.end()) {
            int index = it - ranks.begin();
            result[index] = std::to_string(node.external_id);
        }
    }

    int nid = node.id;
    int pid = node.parent_id;

    while (nid != pid) {
        auto node = map.at(nid);

        auto it = std::find(ranks.begin(), ranks.end(), node.rank);
        if (it != ranks.end()) {
            int index = it - ranks.begin();
            result[index] = std::to_string(node.external_id);
        }
        nid = pid;
        pid = map.at(nid).parent_id;
    }

//    return lineage;
    return Utils::join(result, "|");
}

std::string protal::taxonomy::IntTaxonomy::Lineage(int t, std::string divider)  {
    auto init = map.at(t);

    if (init.IsLeaf() && !init.rep_genome.empty())
        init = map.at(init.parent_id);
    if (init.rank == "genome")
        init = map.at(init.parent_id);

    int nid = init.id;
    int pid = init.parent_id;

    if (nid == pid) return divider;

    std::string tmp = init.scientific_name;
    std::string lineage = init.scientific_name;

    nid = pid;
    pid = map.at(nid).parent_id;

    while (nid != pid) {
        auto node = map.at(nid);

        lineage = node.scientific_name + divider + lineage;

        nid = pid;
        pid = map.at(nid).parent_id;
    }

    return lineage;
}

std::string protal::taxonomy::IntTaxonomy::LineageExternal(int t, std::string divider)  {
    auto init = map.at(t);

    if (init.IsLeaf() && !init.rep_genome.empty())
        init = map.at(init.parent_id);
    if (init.rank == "genome")
        init = map.at(init.parent_id);

    int nid = init.id;
    int pid = init.parent_id;

    if (nid == pid) return divider;

    std::string tmp = std::to_string(init.external_id);
    std::string lineage = std::to_string(init.external_id);

    nid = pid;
    pid = map.at(nid).parent_id;

    while (nid != pid) {
        auto node = map.at(nid);

        lineage = std::to_string(node.external_id) + divider + lineage;

        nid = pid;
        pid = map.at(nid).parent_id;
    }

    return lineage;
}

bool protal::taxonomy::IntTaxonomy::IsRoot(int t) {
    return t == root_id;
}

bool protal::taxonomy::IntTaxonomy::IsLeaf(size_t taxid) {
    return Get(taxid).children.empty();
}

bool protal::taxonomy::IntTaxonomy::HasNode(int taxid) {
    return map.contains(taxid);
}