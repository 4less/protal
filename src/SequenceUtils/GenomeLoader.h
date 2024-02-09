//
// Created by fritsche on 15/07/22.
//

#pragma once
#include <sparse_map.h>
#include <string>
#include <fstream>
#include <err.h>
#include <KmerUtils.h>
#include <ranges>
#include <sparse_set.h>

#include "Utils.h"
#include <sysexits.h>

#include "Benchmark.h"

namespace protal {
    class Gene {
        static const size_t DEFAULT = UINT64_MAX;
    private:
        size_t m_id = 0;
        std::string m_sequence = "";
        size_t m_start_byte;
        size_t m_length = DEFAULT;
        std::ifstream* m_is = nullptr;

        size_t m_short_unique = 0;
        size_t m_long_unique = 0;
        size_t m_long_super_unique = 0;
        size_t m_total_kmers = 0;

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
            return (m_id != 0 || m_length > 0) && m_length != DEFAULT;
        }

        bool IsLoaded() const {
            return !m_sequence.empty();
        }

        void Load() {
            if (!IsLoaded()) {
                Load(m_sequence, m_start_byte, m_length);
            }
        };

        void LoadOMP() {
#pragma omp critical(genome_loader)
            if (!IsLoaded()) {
                Load(m_sequence, m_start_byte, m_length);
            }
        };

        size_t GetLength() const {
            return m_length;
        }

        const std::string& Sequence() const {
            return m_sequence;
        }
        const size_t GetId() const {
            return m_id;
        }

        void SetUniqueValues(size_t short_unique, size_t long_unique, size_t long_super_unique, size_t total_kmers) {
            m_short_unique = short_unique;
            m_long_unique = long_unique;
            m_long_super_unique = long_super_unique;
            m_total_kmers = total_kmers;
        }

        [[nodiscard]] std::tuple<size_t, size_t, size_t, size_t> GetUniqueKmerCounts() const {
            return { m_short_unique, m_long_unique, m_long_super_unique, m_total_kmers };
        }

        double UniqueRate() const {
            return m_long_unique/static_cast<double>(m_total_kmers);
        }

        double SuperUniqueRate() const {
            return m_long_super_unique/static_cast<double>(m_total_kmers);
        }
    };


    class Genome {
    public:
        using GeneKey = size_t;
        using GenomeKey = size_t;
        using GeneList = std::vector<Gene>;

    private:
        using GeneID = uint16_t;
        GeneList m_genes;
        GenomeKey m_key;
        tsl::sparse_set<GeneID> m_hittable_genes;
        bool m_is_loaded = false;

        size_t m_short_unique = 0;
        size_t m_long_unique = 0;
        size_t m_long_super_unique = 0;
        size_t m_total_kmers = 0;



        const size_t GeneKeyToIndex(GeneKey const& key) const {
            return key - 1;
        }
    public:

//        Genome() {};
        explicit Genome(GenomeKey key) : m_key(key) {};

        GenomeKey GetKey() {
            return m_key;
        }

        void SetUniqueValues() {
            for (auto& gene : m_genes) {
                auto [su, lu, lsu, total] = gene.GetUniqueKmerCounts();
                m_short_unique += su;
                m_long_unique += lu;
                m_long_super_unique += lsu;
                m_total_kmers += total;
            }
        }

        [[nodiscard]] std::tuple<size_t, size_t, size_t, size_t> GetUniqueKmerCounts() const {
            return { m_short_unique, m_long_unique, m_long_super_unique, m_total_kmers };
        }

        void AddGene(GeneKey key, size_t id, size_t start_byte, size_t length, std::ifstream* is) {
            auto index = GeneKeyToIndex(key);
            if (index >= m_genes.size()) {
                m_genes.resize(index+1);
            }
            m_genes[index].Set(key, start_byte, length, is);
        }

        void AddHittableGene(GeneID geneid) {
            m_hittable_genes.insert(geneid);
        }

        bool IsGeneHittable(GeneID geneid) const {
            return m_hittable_genes.empty() ? true : m_hittable_genes.contains(geneid);
        }

        std::vector<uint32_t> GetHittableGenes() {
            std::vector<uint32_t> genes;
            for (auto i = 1; i <= 120; i++) {
                if (IsGeneHittable(i)) genes.emplace_back(i);
            }
            return genes;
        }

        void LoadGene(GeneKey key) {
            m_genes[GeneKeyToIndex(key)].Load();
        };

        bool ValidGene(GeneKey key) {
            return GeneKeyToIndex(key) < m_genes.size();
        }

        Gene& GetGene(GeneKey key) {
            return m_genes.at(GeneKeyToIndex(key));
        }

        const GeneList GetGeneList() const {
            return m_genes;
        }

        size_t GeneNum() const {
            return m_hittable_genes.empty() ? std::count_if(m_genes.begin(), m_genes.end(), [](Gene const& gene) {
                return gene.IsSet();
            }) : m_hittable_genes.size();
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
#pragma omp critical(genome_loader)
            {
                if (!IsLoaded()) {
                    LoadGenome();
                }
            }
            return m_genes.at(GeneKeyToIndex(key));
        };
    };


    class GenomeLoader {
        using GenomeKey = Genome::GenomeKey;
        using GenomeMap = tsl::sparse_map<GenomeKey, Genome>;
        using GeneKey = Genome::GeneKey;


        std::string m_path;
        std::string m_genome_map;
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
                m_genome_map(genome_map),
                m_is(genome_path, std::ios::in) {
            LoadPositionMap(genome_map);
        };

        GenomeLoader(const GenomeLoader& other) :
                m_path(other.m_path),
                m_genome_map(other.m_genome_map),
                m_is(other.m_path, std::ios::in) {
            LoadPositionMap(other.m_genome_map);
        }

        ~GenomeLoader() {
            m_is.close();
        }

        void PrintHittableGenes() {
            for (auto& [gid, _] : m_genomes) {
                auto& genome = m_genomes.at(gid);
                auto hg = genome.GetHittableGenes();
                std::string hgstr = "";
                for (auto gene : hg) {
                    hgstr += std::to_string(gene) + ',';
                }
                std::cout << "Hittable\t" << gid << '\t' << hg.size() << '\t' << hgstr << std::endl;
            }
        }

        void LoadUniqueKmers(std::string const& file) {
            std::ifstream is(file, std::ios::in);

            std::vector<std::string> tokens;
            std::string line;
            while (std::getline(is, line)) {
                Utils::split(tokens, line, "\t");


                auto taxid = std::stoull(tokens[0]);
                auto geneid = std::stoull(tokens[1]);

                auto short_unique = std::stoull(tokens[2]);
                auto short_unique_rate = std::stod(tokens[3]);
                auto long_unique = std::stoull(tokens[4]);
                auto long_unique_rate = std::stod(tokens[5]);
                auto long_super_unique = std::stoull(tokens[6]);
                auto long_super_unique_rate = std::stod(tokens[7]);
                auto total_kmers = std::stoull(tokens[8]);

                auto& taxon = m_genomes.at(taxid);
                if (short_unique + long_unique > 0) {
                    taxon.AddHittableGene(geneid);
                }
                auto& gene = taxon.GetGene(geneid);
                gene.SetUniqueValues(short_unique, long_unique, long_super_unique, total_kmers);
            }
            is.close();

            for (auto& tid : views::keys(m_genomes)) {
                m_genomes.at(tid).SetUniqueValues();
            }

            // PrintHittableGenes();
        }

        void LoadHittableGenes(std::string const& file) {
            std::ifstream is(file, std::ios::in);

            std::vector<std::string> tokens;
            std::string line;
            while (std::getline(is, line)) {
                Utils::split(tokens, line, "\t");
                auto taxid = std::stoull(tokens[0]);
                auto gene_str = tokens[2];

                auto& taxon = m_genomes.at(taxid);

                Utils::split(tokens, gene_str, ",");
                for (auto const& g : tokens) {
                    auto gid = std::stoul(g);
                    taxon.AddHittableGene(gid);
                }
            }
            is.close();
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

        void WriteSamHeader(std::ostream& os=std::cout) {
            os << "@HD\tVN:1.6\n";
            for (auto& [key, genome] : m_genomes) {
                auto& genes = genome.GetGeneList();
                for (auto i = 0; i < genes.size(); i++) {
                    if (genes[i].IsSet()) {

                        os << "@SQ\tSN:" << key << '_' << genes[i].GetId() << '\t' << "LN:" << genes[i].GetLength() << '\n';
                    }
                }
            }
        }

        void LoadAllGenomes() {
            std::vector<GenomeKey> keys;
            for (auto& pair : m_genomes) {
                keys.emplace_back(pair.first);
            }
            std::sort(keys.begin(), keys.end());

            for (auto& key : keys) {
                if (!m_genomes.at(key).IsLoaded()){
                    m_genomes.at(key).LoadGenomeOMP();
                }
            }
        }

        bool AllGenomesLoaded() {
            for (auto& [id, genome] : m_genomes) {
                if (!genome.IsLoaded()){
                    return false;
                }
            }
            return true;
        }

        Genome& GetGenome(GenomeKey const& key) {
            assert(m_genomes.contains(key));
            if (!m_genomes.contains(key)) {
                std::cout << "Genomes Key: " << key << std::endl;
                exit(10);
            }
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
