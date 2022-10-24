//
// Created by fritsche on 15/07/22.
//

#pragma once
#include <string>
#include <sparse_map.h>
#include <fstream>
#include <err.h>
#include "Utils.h"
#include <sysexits.h>

namespace protal {
    class Gene {
    private:
        size_t m_id = 0;
        std::string m_sequence = "";
        size_t m_start_byte;
        size_t m_length;
        std::ifstream* m_is = nullptr;

        void Load(std::string& into, size_t start_byte, size_t length) {
            if(m_is && m_is->is_open())
            {
                m_is->seekg(start_byte);
                into.resize(length);
                m_is->read(&into[0], length);
            }
        };

    public:
        Gene(){};

        void Set(size_t id, size_t start_byte, size_t length, std::ifstream* is) {
            m_id = id;
            m_start_byte = start_byte;
            m_length = length;
            m_is = is;
        }

        bool IsSet() const {
            return m_id != 0 || m_length != 0;
        }

        bool IsLoaded() const {
            return !m_sequence.empty();
        }
        void Load() {
            Load(m_sequence, m_start_byte, m_length);
        };
        const std::string& Sequence() const {
            return m_sequence;
        }
        const size_t GetId() const {
            return m_id;
        }
    };


    class Genome {
    public:
        using GeneKey = size_t;
        using GenomeKey = size_t;
        using GeneList = std::vector<Gene>;

    private:
        GeneList m_genes;
        GenomeKey m_key;
        bool m_is_loaded = false;


        const size_t GeneKeyToIndex(GeneKey const& key) const {
            return key - 1;
        }
    public:

//        Genome() {};
        explicit Genome(GenomeKey key) : m_key(key) {};

        GenomeKey GetKey() {
            return m_key;
        }

        void AddGene(GeneKey key, size_t id, size_t start_byte, size_t length, std::ifstream* is) {
            auto index = GeneKeyToIndex(key);
            if (index >= m_genes.size()) {
                m_genes.resize(index+1);
            }
            m_genes[index].Set(key, start_byte, length, is);
        }

        void LoadGene(GeneKey key) {
            m_genes[GeneKeyToIndex(key)].Load();
        };

        Gene& GetGene(GeneKey key) {
            return m_genes.at(GeneKeyToIndex(key));
        }

        const GeneList GetGeneList() const {
            return m_genes;
        }

        size_t GeneNum() const {
            return std::count_if(m_genes.begin(), m_genes.end(), [](Gene const& gene) {
                return gene.IsSet();
            });
        }

        void LoadGenome() {
            std::for_each(m_genes.begin(), m_genes.end(), [](Gene& gene){ gene.Load(); });
            m_is_loaded = true;
        };

        bool IsLoaded() const {
            return m_is_loaded;
//            return std::any_of(m_genes.begin(), m_genes.end(), [](Gene const& gene){ return gene.IsLoaded(); });
        }

        void LoadGenomeOMP() {
#pragma omp critical(genome_loader)
            {
                if (!IsLoaded()) {
                    std::for_each(m_genes.begin(), m_genes.end(), [](Gene &gene) { gene.Load(); });
                    m_is_loaded = true;
                }
            }
        };

        Gene& GetGeneOMP(GeneKey key) {
            if (!IsLoaded()) {
                LoadGenomeOMP();
            }
            return m_genes.at(GeneKeyToIndex(key));
        };
    };


    class GenomeLoader {
        using GenomeKey = Genome::GenomeKey;
        using GenomeMap = tsl::sparse_map<GenomeKey, Genome>;
        using GeneKey = Genome::GeneKey;


        std::string m_path;
        std::ifstream m_is;
        GenomeMap m_genomes;

        Genome& AddOrGetGenome(GenomeKey const& key) {
            if (!m_genomes.contains(key)) {
                m_genomes.insert( { key, Genome(key) } );
            }
            return m_genomes.at(key);
        }

    public:
        GenomeLoader(std::string genome_path, std::string genome_map) :
                m_path(genome_path),
                m_is(genome_path, std::ios::in) {
            LoadPositionMap(genome_map);
        };

        ~GenomeLoader() {
            m_is.close();
        }

        size_t GetLoadedGenomeCount() const {
            return std::count_if(m_genomes.begin(), m_genomes.end(), [](auto const& pair) { return pair.second.IsLoaded(); });
        }

        const Genome& GetGenome(GenomeKey const& key) const {
            return m_genomes.at(key);
        }

        GenomeMap& GetGenomeMap() {
            return m_genomes;
        }

        void LoadAllGenomes() {
            std::vector<GenomeKey> keys;
            for (auto& pair : m_genomes) {
                keys.emplace_back(pair.first);
            }
            std::sort(keys.begin(), keys.end());

            for (auto& key : keys) {
                m_genomes.at(key).LoadGenome();
            }
        }

        Genome& GetGenome(GenomeKey const& key) {
            assert(m_genomes.contains(key));
            return m_genomes.at(key);
        }

        void LoadPositionMap(std::string file_path) {
            std::ifstream is(file_path, std::ios::in);
            std::string line;

            while (std::getline(is, line)) {
                auto tokens = Utils::split(line, "\t");


                if (tokens.size() != 4) {
                    errx(EX_DATAERR, "Position map file needs to have exactly 4 columns (%s line has %zul).", line.c_str(), tokens.size());
                }

                GenomeKey genome_id = std::stoull(tokens[0]);
                GeneKey gene_key = std::stoull(tokens[1]);
                size_t start = std::stoull(tokens[2]);
                size_t length = std::stoull(tokens[3]) - start;

                auto& genome = AddOrGetGenome(genome_id);

                genome.AddGene(gene_key, gene_key, start, length, &m_is);
            }
            is.close();
        }
    };
}