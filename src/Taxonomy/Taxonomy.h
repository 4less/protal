//
// Created by fritsche on 17/05/2021.
//

#pragma once

#include "sparse_map.h"
#include "sparse_set.h"

namespace protal::taxonomy {
    struct Node;
    struct IntNode;
    struct Name;

    typedef int32_t TaxId;
    typedef int32_t RankId;

    typedef tsl::sparse_map<std::string, RankId> RankToId;
    typedef tsl::sparse_map<RankId, std::string> RankFromId;

    typedef tsl::sparse_map<TaxId, Node*> NodeFromId;
    typedef tsl::sparse_map<TaxId, IntNode> IntNodeMap;
    typedef tsl::sparse_map<std::string, TaxId> StringIdMap;
//    typedef tsl::sparse_map<TaxId, tsl::sparse_map<std::string, Name>> NamesFromId;/
    typedef std::vector<Name> Names;
    typedef std::vector<Node> Nodes;

    struct Name {
        int id = -1;
        std::string name = "";
        std::string name_unique = "";
        std::string name_class = "";

        std::string ToString(std::string sep) {
            return std::to_string(id) + sep + name + sep + name_unique + sep + name_class;
        }
    };

    struct Node {
        TaxId id = -1;
        TaxId external_id = -1;
        TaxId ncbi_id = -1;

        Node* parent_node = nullptr;
        std::vector<Node*> children;
        std::string rank = "";
        std::string rep_genome = "";

        uint32_t level = -1;
        uint32_t inorder_id = 0;

        Name* name = nullptr;

        ~Node() {
            delete name;
        };

        void RemoveFromTree() {
            // delete from parents children
            if (parent_node)
                parent_node->RemoveChild(this);
            for (auto child : children) {
                child->parent_node = nullptr;
            }
        }

        void RemoveRelink() {
            for (auto child : children) {
                parent_node->children.push_back(child);
                child->parent_node = parent_node;
            }
            parent_node->RemoveChild(this);
            parent_node = nullptr;
            children.clear();
        }

        void RemoveChild(Node* node) {
            children.erase(std::remove(children.begin(), children.end(), node), children.end());
        }

        bool IsRoot() {
            return !parent_node;
        }

        bool IsLeaf() {
            return children.empty();
        }

        void ChangeId(int new_id) {
            external_id = id;
            id = new_id;
        }
    };


    class StdTaxonomy {
        RankToId rank_to_id_;
        RankFromId rank_from_id_;

        NodeFromId node_from_id_;

        // This is the master list of nodes. If it gets removed from this vector
        // the node needs to be deleted from memory
        std::vector<Node*> nodes_;
        Names names_;

        // Root id
        TaxId root_id_ = -1;
        std::string separator = "\t";

        // Initialize rank maps
        void InitRankMaps();
        void InitRankMaps(std::string path);

        // helper functions
        Node *GetNode(int taxid);

        void DeleteNode(Node *node);
        void Level();
        void Level(Node *, int level);

    public:
        StdTaxonomy();
        ~StdTaxonomy();

        Node* LCA(Node* t1, Node* t2);
        Node* LCA(TaxId t1, TaxId t2);

        std::vector<Node*>& Nodes();

        void ClearNodes();
        void KeepRanksAndLeaves(tsl::sparse_set<std::string> &ranks, std::vector<Node*> &remove_list, Node *parent_node, Node *node);
        void KeepRanksAndLeaves(Node *parent_node, Node *node);
        void RepopulateNodes(Node* node=nullptr);

        void LoadNodes(std::string path);
        void LoadNames(std::string path);
        void Add(std::string path);
        void SaveCustom(std::string path);
        void Export(std::string path);

        // Subset and relabel some taxa based on a tsv file species to id
        void SubsetAndRelabel(std::string path, std::string collapse_level="");

        void SanityCheck(Node *node=nullptr);
    };

    struct IntNode {
        int id;
        int parent_id = -1;
        int external_id;
        int level;
        std::string scientific_name = "";
        std::string rank;
        std::string rep_genome = "";
        std::vector<int> children;

        IntNode& operator = (const IntNode &t)
        {
            id = t.id;
            parent_id = t.parent_id;
            external_id = t.external_id;
            level = t.level;
            scientific_name = t.scientific_name;
            rank = t.rank;
            children = t.children;
            rep_genome = t.rep_genome;
            return *this;
        }

        bool IsLeaf() {
            return children.empty();
        }

        std::string ToString() const {
            std::string s;

            s += std::to_string(id);
            s += '\t';
            s += std::to_string(parent_id);
            s += '\t';
            s += std::to_string(external_id);
            s += '\t';
            s += scientific_name;
            s += '\t';
            s += rank;
            s += '\t';
            s += std::to_string(level);

            return s;
        }

    };

    class IntTaxonomy {
        void Load(std::string path);

    public:
        StringIdMap string_to_id;
        IntNodeMap map;
        TaxId root_id = 0;

        IntTaxonomy(std::string path);
        int LCA(int t1, int t2);
        IntNode &Get(int t1);
        IntNode &Get(int t1, std::string rank);
        TaxId Get(std::string taxon);
        std::string Lineage(int t, std::string divider=";");
        bool IsRoot(int t);

        bool IsLeaf(size_t taxid);

        IntNode &GetParent(int t);

        size_t MaxTaxid();
        size_t MaxLeafId();

        bool IsNodeAncestor(int node_id, int leaf_id);

        size_t GetSubtreeSize(int t);

        std::string LineageExternal(int t, std::string divider=";");

        std::string LineageStr(int t, const std::vector<std::string> ranks={ "phylum", "class", "order", "family", "genus", "species" }, std::string divider=";");

        std::string LineageExternalIds(int t, const std::vector<std::string> ranks, std::string divider);

        static const int Rank2Id(std::string rank) {
            if (rank == "domain" || rank == "superkingdom") {
                return 0;
            } else if (rank == "phylum") {
                return 1;
            } else if (rank == "class") {
                return 2;
            } else if (rank == "order") {
                return 3;
            } else if (rank == "family") {
                return 4;
            } else if (rank == "genus") {
                return 5;
            } else if (rank == "species") {
                return 6;
            } else if (rank == "subspecies") {
                return 7;
            } else if (rank == "genome") {
                return 8;
            }
            return -1;
        }

        bool HasNode(int taxid);
    };

}