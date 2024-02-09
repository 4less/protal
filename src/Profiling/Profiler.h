//
// Created by fritsche on 07/10/22.
//

#pragma once

#include "Strain.h"
#include "AlignmentUtils.h"
#include <cstdint>
#include <cstddef>
#include <vector>
#include <sparse_map.h>
#include <numeric>
#include <unordered_set>
#include <Profiler/ProfilerDefinitions.h>

#include "LineSplitter.h"
#include "Taxonomy.h"
#include "InternalReadAlignment.h"
#include "Constants.h"
#include "SNPUtils.h"
#include "ProfilerOptions.h"
#include "ScoreAlignments.h"
#include "gzstream/gzstream.h"
#include "cPMML.h"
#include "sparse_map.h"
#include "Benchmark.h"

namespace protal {
    bool IsDigit(std::string &test) {
        return false;
    }

    using TruthSet = tsl::robin_set<uint32_t>;
    TruthSet GetTruth(std::string const& file_path, taxonomy::IntTaxonomy& taxonomy) {
        std::ifstream is(file_path, std::ios::in);
        TruthSet truths;
        std::vector<std::string> tokens;
        std::string delim = "\t";
        int lineage_column = -1;


        std::string line;
        while (std::getline(is, line)) {
            LineSplitter::Split(line, delim, tokens);
            if (lineage_column == -1) {
                for (auto i = 0; i < tokens.size(); i++) {
                    if (tokens[i].starts_with("d__") || tokens[i].starts_with("s__")) {
                        lineage_column = i;
                        break;
                    }
                }
            }


            auto lineage_str = tokens[lineage_column];
            LineSplitter::Split(lineage_str, ";", tokens);
            std::string species = "";
            for (auto& e : tokens) {
                if (e.starts_with("s__")) {
                    species = e;
                }
            }
            if (species == "s__") {
                std::cerr << "No Species classification in Gold Profile: " << line << std::endl;
                continue;
            }
            truths.insert(taxonomy.Get(species));
        }
        return truths;
    }

    namespace profiler {
        // Forward declaration
        class TaxonFilterForest;
        class TaxonFilter;

        struct MicrobialProfile;
        struct Gene;
        struct Taxon;

        using GeneMap = tsl::sparse_map<uint32_t, Gene>;
        using TaxonMap = tsl::sparse_map<uint32_t, Taxon>;
        using GeneRef = protal::Gene;

        using TaxonFilterObj = TaxonFilterForest;

        // Data structures

        struct AlleleCounts {
            uint32_t filtered = 0;
            uint32_t mono = 0;
            uint32_t bi = 0;
            uint32_t tri = 0;
            uint32_t tetra = 0;

            uint32_t Multi() {
                return bi + tri + tetra;
            }

            uint32_t Noisy() {
                return Multi() + filtered;
            }

            uint32_t AllValid() {
                return Multi() + mono;
            }

            uint32_t All() {
                return Multi() + mono;
            }
        };

        class Gene {
        public:
//            using SNPs = std::vector<SNP>;

            size_t m_mapped_reads = 0;
            size_t m_mapped_length = 0;
            size_t m_total_kmers = 0;
            size_t m_mapq_sum = 0;
            size_t m_unique_mers = 0;
            size_t m_unique_two_mers = 0;
            size_t m_unique_mer_reads = 0;
            size_t m_unique_two_mers_reads = 0;
            size_t m_unique_two_mer_reads = 0;
            double m_ani_sum = 0;

            size_t m_gene_length = 0;
            std::vector<SamEntry*> m_sams;
            StrainLevelContainer m_strain_level;


            GeneRef* m_gene_ref = nullptr;

            Gene(protal::Gene& gene_ref) :
                    m_gene_ref(&gene_ref),
                    m_strain_level(gene_ref) {
            };

            Gene& operator=(const Gene& other) {
                m_gene_ref = other.m_gene_ref;
                return *this;
            }

            size_t Coverage(size_t above=0) {
                auto cov_vec = GetStrainLevel().GetSequenceRangeHandler().CalculateCoverageVector2();
                auto cov = ranges::count_if(cov_vec, [above](const uint16_t e){ return e > above; });
                return cov;
            }

            AlleleCounts AlleleSNPCounts(uint32_t min_cov, size_t min_qual_sum) {
                std::vector<size_t> alleles(5, 0);
                GetAlleles(alleles, min_cov, min_qual_sum);
                AlleleCounts counts;
                counts.filtered = alleles[0];
                counts.mono = alleles[1];
                counts.bi = alleles[2];
                counts.tri = alleles[3];
                counts.tetra = alleles[4];

                return counts;
            }

            void AddRead(GenePos pos) {
                m_mapped_reads++;
            }

            void ClearSams() {
                m_sams.clear();
            }

            void GetAlleles(std::vector<size_t>& allele_counts, uint32_t min_cov=0, uint32_t min_qual_sum=0) const {
                auto& vars = m_strain_level.GetVariantHandler().GetVariants();

                for (const auto& [pos, _] : vars) {
                    auto const& var = vars.at(pos);
                    size_t total_cov = std::accumulate(var.begin(), var.end(), 0, [](size_t acc, Variant const& var) {
                        return acc + var.Observations();
                    });

                    size_t valid_alleles = 0;
                    for (auto& v : var) {
                        valid_alleles += (v.Observations() >= min_cov & v.QualitySum() >= min_qual_sum);
                    }
                    if (valid_alleles >= allele_counts.size()) {
                        allele_counts.resize(valid_alleles+1, 0);
                    }
                    allele_counts[valid_alleles]++;
                }
            }

            const StrainLevelContainer& GetStrainLevel() const {
                return m_strain_level;
            }

            StrainLevelContainer& GetStrainLevel() {
                return m_strain_level;
            }

            std::string GetStatisticsString() {
                std::string s = "";
                auto cov_vec = m_strain_level.GetSequenceRangeHandler().GetCoverageVector();
                auto cov_sum = std::accumulate(cov_vec.begin(), cov_vec.end(), 0);
                auto cov = m_strain_level.GetSequenceRangeHandler().CoveredPortion();
                s += std::to_string(m_mapped_reads) + '\t';
                s += std::to_string(cov_sum) + '\t';
                s += std::to_string(m_gene_length) + '\t';
                s += std::to_string(m_ani_sum) + '\t';
                s += std::to_string(m_mapq_sum) + '\t';
                s += std::to_string(m_ani_sum/static_cast<double>(m_mapped_reads)) + '\t';
                s += std::to_string(m_mapq_sum/static_cast<double>(m_mapped_reads)) + '\t';
                s += std::to_string(m_strain_level.GetSequenceRangeHandler().CoveredPortion()) + '\t';
                s += std::to_string(m_strain_level.GetSequenceRangeHandler().CoveredPortion()/static_cast<double>(m_gene_length));
//                s += std::to_string(m_strain_level.GetSequenceRangeHandler().CoveredPortion(5)) + '\t';
//                s += std::to_string(m_strain_level.GetSequenceRangeHandler().CoveredPortion(5)/static_cast<double>(m_gene_length));
                return s;
            }

            bool AddSam(SamEntry const& sam, size_t read_id, double ani=0.0, bool store_sam=true, bool no_strain=true) {
                if (!no_strain) {
                    auto successful = m_strain_level.AddSam(sam, read_id);
                    if (!successful) {
                        return false;;
                    }
                }

                if (store_sam) {
                    m_sams.emplace_back(const_cast<SamEntry*>(&sam));
                }


                // if (extract_snps) {
                //     ExtractSNPs(sam, m_gene_ref->Sequence(), m_snps, 0, 0, read_id);
                // }

                size_t length = AlignmentLengthRef(sam.m_cigar);

                m_mapped_reads++;
                m_mapped_length += length;
                m_total_kmers += (length - 30) * 0.2;
                m_mapq_sum += sam.m_mapq;
                m_ani_sum += ani;

                m_unique_mers += sam.m_uniques;
                m_unique_two_mers += sam.m_uniques_two;

                m_unique_mer_reads += sam.m_uniques > 0;
                m_unique_two_mer_reads += sam.m_uniques_two > 0;

                return true;
            }

            void SetLength(size_t length) {
                m_gene_length = length;
            }

            size_t UniqueKmers() const {
                return m_unique_mers;
            }

            size_t UniqueKmersReads() const {
                return m_unique_mer_reads;
            }

            size_t UniqueTwoKmers() const {
                return m_unique_two_mers;
            }

            size_t UniqueTwoKmersReads() const {
                return m_unique_two_mers_reads;
            }

            size_t TotalKmers() const {
                return m_total_kmers;
            }

            double UniqueRate() const {
                return m_unique_mers / static_cast<double>(m_total_kmers);
            }

            double SuperUniqueRate() const {
                return m_unique_two_mers / static_cast<double>(m_total_kmers);
            }


            double VerticalCoverage() const {


                if (m_mapped_length/10 < m_sams.size()) {


                    std::cout << "VerticalCoverage()" << std::endl;
                    std::cout << static_cast<double>(m_mapped_length) << "/" << static_cast<double>(m_gene_length) << std::endl;
                    std::cout << "Read num: " << m_sams.size() << std::endl;
                    std::cout << "Mapped reads: " << m_mapped_reads << std::endl;
                    std::cout << "Gene length: " << m_gene_ref->Sequence().length() << std::endl;

                    for (auto sam : m_sams) {
                        std::cout << sam->ToString() << std::endl;
                    }

                    Utils::Input();
                }
                return static_cast<double>(m_mapped_length)/static_cast<double>(m_gene_length);
            }

//            SNPs& GetSNPs() {
//                return m_snps;
//            }
//
//            const SNPs& GetSNPs() const {
//                return m_snps;
//            }

            std::string ToString() const {
                std::string str = "{";
                str += std::to_string(m_mapped_reads) + ",";
                str += std::to_string(static_cast<double>(m_mapq_sum)/m_mapped_reads) + ",";
                str += std::to_string(m_ani_sum/m_mapped_reads) + "}";


                return str;
            }

            std::string ToString2() const {
                std::string str;
                str += std::to_string(m_mapped_reads) + "\t";
                str += std::to_string(static_cast<double>(m_mapq_sum)/m_mapped_reads) + "\t";
                str += std::to_string(m_ani_sum/m_mapped_reads) + "\t";
//                str += std::to_string(m_snps.size());
                return str;
            }
        };

        class Taxon {
        private:
            size_t m_id;
            std::string m_name;
            size_t m_total_hits = 0;
            size_t m_total_kmers = 0;
            size_t m_unique_mers = 0;
            size_t m_unique_mer_reads = 0;
            size_t m_unique_hits = 0;
            double m_ani_sum = 0;
            size_t m_mapq_sum = 0;
            double m_vcov = -1;
            GeneMap m_genes;

            Genome m_genome;
            size_t m_genome_gene_count = 0;

        public:

            Taxon(Genome& genome) : m_genome(genome), m_genome_gene_count(genome.GeneNum()) {}

            void AddHit(GeneId geneid, GenePos genepos, double ani, bool unique) {
                if (!m_genes.contains(geneid)) {
                    auto& g = m_genome.GetGene(geneid);
                    g.LoadOMP();
                    m_genes.insert( { geneid, profiler::Gene(g) } );
                    m_genes.at(geneid).SetLength(m_genome.GetGene(geneid).Sequence().length());
                }

                m_genes.at(geneid).AddRead(genepos);
                m_total_hits++;
                m_ani_sum += ani;
                m_unique_hits += unique;
            }

            std::vector<size_t> GetAlleles(size_t min_cov=0, size_t min_qual_sum=0) const {
                std::vector<size_t> alleles(5, 0);
                for (auto& [id, _] : GetGenes()) {
                    auto& gene = GetGenes().at(id);
                    gene.GetAlleles(alleles, min_cov, min_qual_sum);
                }
                alleles.resize(5,0);
                return alleles;
            }

            const GeneMap& GetGenes() const {
                return m_genes;
            }

            GeneMap& GetGenes() {
                return m_genes;
            }

            void ClearSams() {
                for (auto& [id, _] : m_genes) m_genes.at(id).ClearSams();
            }

            bool AddSam(GeneId geneid, SamEntry const& sam, double score, bool unique, size_t read_id, bool no_strain=true) {
                if (!m_genes.contains(geneid)) {
                    auto& g = m_genome.GetGene(geneid);
                    g.LoadOMP();
                    m_genes.insert( { geneid, profiler::Gene(g) } );
                    m_genes.at(geneid).SetLength(m_genome.GetGene(geneid).Sequence().length());
                }



                bool success = m_genes.at(geneid).AddSam(sam, read_id, score, true, no_strain);
                if (!success) return false;

                m_unique_mers += sam.m_uniques;
                m_unique_mer_reads += sam.m_uniques > 0;
                m_total_hits++;
                m_total_kmers += (sam.m_seq.length() - 30) * 0.2; //TODO: store that info in sam
                m_ani_sum += score;
                m_mapq_sum += sam.m_mapq;
                m_unique_hits += unique;
                return true;
            }

            size_t GetGenomeGeneNumber() const {
                return m_genome.GeneNum();
            }

            size_t UniqueKmers() const {
                return std::accumulate(m_genes.begin(), m_genes.end(), 0, [](auto acc, auto const& pair) { return acc + pair.second.UniqueKmers(); });
            }

            size_t UniqueTwoKmers() const {
                return std::accumulate(m_genes.begin(), m_genes.end(), 0, [](auto acc, auto const& pair) { return acc + pair.second.UniqueTwoKmers(); });
            }

            size_t GenesWithUniqueKmers() const {
                return std::count_if(m_genes.begin(), m_genes.end(), [](auto const& pair) { return pair.second.UniqueKmers() > 0; });
            }

            size_t GenesWithUniqueTwoKmers() const {
                return std::count_if(m_genes.begin(), m_genes.end(), [](auto const& pair) { return pair.second.UniqueTwoKmers() > 0; });
            }

            Genome& GetGenome() {
                return m_genome;
            }

            size_t PresentGenes() const {
                return m_genes.size();
            }

            double Uniqueness() const {
                return static_cast<double>(m_unique_hits)/static_cast<double>(m_total_hits);
            }

            double GetMeanANI() const {
                return m_ani_sum/m_total_hits;
            }

            double GetMeanMAPQ() const {
                return m_mapq_sum/m_total_hits;
            }

            size_t TotalHits() const {
                return m_total_hits;
            }

            size_t TotalLength() const {
                return std::accumulate(m_genes.begin(), m_genes.end(), 0, [](size_t acc, std::pair<uint32_t, Gene> const& pair){
                    return acc + pair.second.m_mapped_length;
                });
            }

            double GetGeneVariance(size_t min_gene_hits = 1) const {
                std::vector<double> vec = VerticalCoverageVector();
                std::vector<double> vec_filtered;
                std::copy_if(vec.begin(), vec.end(), std::back_inserter(vec_filtered),
                 [min_gene_hits](double value) { return value >= min_gene_hits; });

                auto var = Variance(vec_filtered);
                return var;
            }

            size_t UniqueHits() const {
                return m_unique_hits;
            }

            void SetId(size_t id) {
                m_id = id;
            }

            void SetName(std::string name) {
                m_name = name;
            }

            std::vector<double> VerticalCoverageVector() const {
                std::vector<double> vcovs;
                for (auto& [geneid, gene] : m_genes) {
                    vcovs.emplace_back(gene.VerticalCoverage());
                }
                return vcovs;
            }



            double Variance(std::vector<double> v) const {
                if (v.size() < 2) return 500.0f;
                double sum = std::accumulate(v.begin(), v.end(), 0.0);
                double mean = sum / v.size();

                std::vector<double> diff(v.size());
                std::transform(v.begin(), v.end(), diff.begin(),
                               std::bind2nd(std::minus<double>(), mean));
                double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
                double stdev = std::sqrt(sq_sum / v.size());
                return stdev;
            }

            double StandardDeviation(std::vector<double> v) const {
                return std::sqrt(Variance(v));
            }

            double VCovStdDev() const {
                auto v = VerticalCoverageVector();
                double sum = std::accumulate(v.begin(), v.end(), 0.0);
                double mean = sum / v.size();

                std::vector<double> diff(v.size());
                std::transform(v.begin(), v.end(), diff.begin(),
                               std::bind2nd(std::minus<double>(), mean));
                double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
                double stdev = std::sqrt(sq_sum / v.size());
                return stdev;
            }

            static double Median(std::vector<double> const& v) {
                if (v.empty()) return 0;
                size_t mid_index = v.size()/2;
                return v.size() % 2 == 1 ?
                    v.at(mid_index) :
                    (v.at(mid_index - 1) + v.at(mid_index)) / 2;
            }

            double VerticalCoverage(bool force=false) {
                if (m_vcov == -1 || force) {
                    std::vector<double> vcovs;
                    for (auto& [geneid, gene] : m_genes) {
                        vcovs.emplace_back(gene.VerticalCoverage());
                    }
                    std::sort(vcovs.begin(), vcovs.end());

                    if (vcovs.empty()) return 0.0;

                    m_vcov = Median(vcovs);
                }
                return m_vcov;
            }

            double VerticalCoverage() const {
                return m_vcov;
            }

            Gene& GetGene(int geneid) {
                if (!m_genes.contains(geneid)) {
                    std::cerr << geneid << " not in genes " << std::endl;
                    exit(9);
                }
                return m_genes.at(geneid);
            }

            double GetAbundance(double total_vertical_coverage) {
                return VerticalCoverage()/total_vertical_coverage;
            }


            std::string ToString() const {
                std::string str;

                str += std::to_string(m_id);
                str += "\t{ ";
                str += "Genes: " + std::to_string(m_genes.size()) + ", ";
                str += "Total Hits: " + std::to_string(m_total_hits) + ", ";
                str += "MeanMAPQ: " + std::to_string(GetMeanMAPQ()) + ", ";
                str += "VCOV: " + std::to_string(VerticalCoverage()) + ", ";
                str += " }";

                return str;
            }

            std::string ToString() {
                std::string str;

                str += std::to_string(m_id);
                str += "\t{ ";
                str += "Genes: " + std::to_string(m_genes.size()) + ", ";
                str += "Total Hits: " + std::to_string(m_total_hits) + ", ";
                str += "MeanMAPQ: " + std::to_string(GetMeanMAPQ()) + ", ";
                str += "VCOV: " + std::to_string(VerticalCoverage(true)) + ", ";
                str += " }";

                return str;
            }

            std::string& GetName() {
                return m_name;
            }

            const std::string& GetName() const {
                return m_name;
            }

            std::string ToString(taxonomy::IntTaxonomy& taxonomy) const {
                std::string str;
                auto& genes = m_genome.GetGeneList();
                auto total_gene_length = std::accumulate(genes.begin(), genes.end(), 0, [](size_t acc, protal::Gene const& gene){
                    return acc + (gene.IsSet() ? gene.GetLength() : 0);
                });

                str += taxonomy.Get(m_id).scientific_name;
                str += "\t" + std::to_string(m_id) + "\t";
                str += "{ ";
                str += "Genes: " + std::to_string(m_genes.size()) + ", ";
                str += "HittableGenes: " + std::to_string(m_genome.GeneNum()) + ", ";
                str += "Total Hits: " + std::to_string(m_total_hits) + ", ";
                str += "Total MGenome: " + std::to_string(total_gene_length) + ", ";
                str += "MeanANI: " + std::to_string(m_ani_sum/(double)m_total_hits) + ", ";
                str += "MeanMAPQ: " + std::to_string(GetMeanMAPQ()) + ", ";
                str += "VCOV: " + std::to_string(VerticalCoverage()) + ", ";
                str += " }";

//                for (auto& [id, gene] : m_genes) {
//                    auto& variant_handler = gene.GetStrainLevel().GetVariantHandler();
//                    str += "\tGene\t" + std::to_string(id) + "\t" + gene.ToString2() + '\n';
//                    for (auto& [key, bin] : variant_handler.GetVariants()) {
//                        str += std::to_string(key) + '\t' + variant_handler.VariantBinToString(bin) + '\n';
////                        str += "\t\t" + snp.ToString() + '\n';
//                    }
//                }

                return str;
            }
        };



        class TaxonFilter {
            double m_min_mean_ani = 0.94;
            double m_min_gene_presence = 0.7;
            size_t m_min_reads = 50;
            size_t m_min_mean_mapq = 5;
            double m_min_uniqueness = 0.1;
        public:
            TaxonFilter() {}
            TaxonFilter(double min_mean_ani, double min_gene_presence, double min_reads, size_t min_mean_mapq) :
                    m_min_mean_ani(min_mean_ani),
                    m_min_gene_presence(min_gene_presence),
                    m_min_reads(min_reads),
                    m_min_mean_mapq(min_mean_mapq) {
            }

            static double ExpectedGenes(size_t marker_genes_count, size_t mapped_reads) {
                return (1 - pow((double) (marker_genes_count - 1)/marker_genes_count, mapped_reads)) * marker_genes_count;
            }

            static double ExpectedGenePresence(Taxon const& taxon) {
                return ExpectedGenes(taxon.GetGenomeGeneNumber(), taxon.TotalHits());
            }

            static double ExpectedGenePresenceRatio(Taxon taxon) {
                return static_cast<double>(taxon.PresentGenes()) / ExpectedGenes(taxon.GetGenomeGeneNumber(), taxon.TotalHits());
            }

            bool Formula1(Taxon const& taxon) const {
                bool pass = // Currently mix of expected gene presence, min ani and min gene presence
                        taxon.GetMeanANI() >= m_min_mean_ani &&
                        taxon.TotalHits() >= m_min_reads &&
                        ExpectedGenePresenceRatio(taxon) >= m_min_gene_presence &&
                        taxon.GetMeanMAPQ() >= m_min_mean_mapq;
                return pass;
            }

            void PrintViolatingValues(Taxon const& taxon) const {
                if (taxon.GetMeanANI() < m_min_mean_ani) {
                    std::cout << "MeanANI:      " << taxon.GetMeanANI() << " < " << m_min_mean_ani << " (min)" << std::endl;
                }
                if (taxon.TotalHits() < m_min_reads) {
                    std::cout << "TotalHits:    " << taxon.TotalHits() << " < " << m_min_reads << " (min)" << std::endl;
                }
                if (ExpectedGenePresenceRatio(taxon) < m_min_gene_presence) {
                    std::cout << "GenePresence: " << ExpectedGenePresenceRatio(taxon) << " < " << m_min_gene_presence << " (min)" << std::endl;
                }
                if (taxon.GetMeanMAPQ() < m_min_mean_mapq) {
                    std::cout << "MeanMAPQ: " << taxon.GetMeanMAPQ() << " < " << m_min_mean_mapq << " (min)" << std::endl;
                }
            }

            bool Pass(Taxon const& taxon) const {
                return Formula1(taxon);
            }
        };

        class HittableGeneMask {
        private:
            tsl::sparse_map<Profiler::TaxonID, tsl::sparse_set<Profiler::GeneID>> m_hittable_map;

        public:
            void Load(std::string const& file) {
                std::ifstream is(file, std::ios::in);

                std::vector<std::string> tokens;
                std::string line;
                while (std::getline(is, line)) {
                    Utils::split(tokens, line, "\t");
                    auto taxid = std::stoull(tokens[0]);
                    auto gene_str = tokens[2];

                    Utils::split(tokens, gene_str, ",");
                    for (auto const& g : tokens) {
                        auto gid = std::stoul(g);
                        m_hittable_map[taxid].insert(gid);
                    }
                }
            }

            size_t HittableGenes(Profiler::TaxonID taxid) {
                return m_hittable_map.at(taxid).size();
            }

            bool IsHittable(Profiler::TaxonID taxid, Profiler::GeneID geneid) {
                return m_hittable_map[taxid].contains(geneid);
            }
        };

        class TaxonFilterForest {
            cpmml::Model m_model{};
            double m_threshold = 0.5;

            mutable std::unordered_map<std::string, std::string> m_sample;


        public:

            TaxonFilterForest(const std::string& path) : m_model(path) {}

            TaxonFilterForest(const TaxonFilterForest& other) :
                    m_model(other.m_model), m_threshold(other.m_threshold), m_sample() {
            }

            TaxonFilterForest(const TaxonFilterForest&& other) :
                    m_model(other.m_model),
                    m_threshold(other.m_threshold),
                    m_sample(other.m_sample) {
            }

            bool Pass(Taxon const& taxon) const {
                auto a = taxon.GetAlleles();
                auto af = taxon.GetAlleles(2, 60);
                // a.resize(5, 0);
                // af.resize(5, 0);

                // "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence"
                m_sample["present_genes"] = std::to_string(taxon.PresentGenes());
                m_sample["total_hits"] = std::to_string(taxon.TotalHits());
                m_sample["unique_hits"] = std::to_string(taxon.UniqueHits());
                m_sample["expected_gene_presence"] = std::to_string(TaxonFilter::ExpectedGenePresence(taxon));
                m_sample["mean_ani"] = std::to_string(taxon.GetMeanANI());
                //m_sample["mean_mapq"] = std::to_string(taxon.GetMeanMAPQ());

                m_sample["variance1"] = std::to_string(taxon.GetGeneVariance(1));
                m_sample["stddev"] = std::to_string(taxon.VCovStdDev());


                m_sample["A1"] = std::to_string(a[1]);
                m_sample["A2"] = std::to_string(a[2]);
                m_sample["AF0"] = std::to_string(af[0]);
                m_sample["AF1"] = std::to_string(af[1]);
                auto prediction_str = m_model.predict(m_sample);
                return prediction_str == "TRUE";
            }
        };

        class MicrobialProfile {
        public:
            Benchmark bm_add_sam{"Add Sam profile"};
            MicrobialProfile(GenomeLoader& genome_loader) : m_genome_loader(genome_loader) {}

            void AddRead(InternalReadAlignment const& ira, bool unique=true) {
                if (!m_taxa.contains(ira.taxid)) {
                    auto& genome = m_genome_loader.GetGenome(ira.taxid);
                    m_taxa.insert( { ira.taxid, Taxon(genome) } );
                    m_taxa.at(ira.taxid).SetName(std::to_string(ira.taxid));
                }
                auto& taxon = m_taxa.at(ira.taxid);
                taxon.AddHit(ira.geneid, ira.genepos, ira.alignment_ani, unique);
            }

            bool AddSam(int taxid, int geneid, SamEntry const& sam, double score, bool unique=true, int read_id=0, bool no_strain=true) {
                if (!m_genome_loader.GetGenome(taxid).IsGeneHittable(geneid)) {
                    return false;
                }

                if (!m_taxa.contains(taxid)) {
                    auto &genome = m_genome_loader.GetGenome(taxid);

                    m_taxa.insert( { taxid, Taxon(genome) } );
                    m_taxa.at(taxid).SetId(taxid);
                    m_taxa.at(taxid).SetName(std::to_string(taxid));
                }

                // Debug
                if (!m_taxa.contains(taxid)) {
                    std::cout << "Error: " << taxid << std::endl;
                    std::cout << sam.ToString() << std::endl;
                }
                auto& taxon = m_taxa.at(taxid);
                unique = sam.m_mapq > 20;
                bm_add_sam.Start();
                bool success = taxon.AddSam(geneid, sam, score, unique, read_id, no_strain);
                bm_add_sam.Stop();
                return success;
            }

            void SetName(std::string name) {
                m_name = name;
            }

            const std::string& GetName() const {
                return m_name;
            }

            std::string& GetName() {
                return m_name;
            }

            void PostProcessSNPs() {
                size_t min_observations = 3;
                size_t min_observations_fwdrev = 2;
                double min_frequency = 0.2;
                size_t min_avg_quality = 10;


                for (auto& [tid, _] : m_taxa) {
                    auto& taxon = m_taxa.at(tid);
                    for (auto& [gid, __] : taxon.GetGenes()) {
                        auto& gene = taxon.GetGenes().at(gid);
                        auto& strain_handler = gene.GetStrainLevel();
                        strain_handler.PostProcess(min_observations, min_observations_fwdrev, min_frequency, min_avg_quality);
                    }
                }
            }

            std::unordered_set<size_t> GetKeySet(std::optional<TaxonFilter> filter={}) {
                std::unordered_set<size_t> keys;
                for (auto& [id, taxon] : m_taxa) {
                    if (filter.has_value() && !filter->Pass(taxon)) continue;
                    keys.insert(id);
                }
                return keys;
            }

            std::string ToString(taxonomy::IntTaxonomy& taxonomy, std::optional<TaxonFilterObj> filter={}) {
                std::string str;

                str += m_name + '\t';
                str += std::to_string(m_taxa.size());
                if (filter.has_value()) {
                    str += " (Filtered)";
                }
                str += '\n';
                for (auto [id, _] : m_taxa) {
                    auto& taxon = m_taxa.at(id);
                    if (filter.has_value() && !filter->Pass(taxon)) continue;
                    taxon.VerticalCoverage(true);
                    str += taxon.ToString(taxonomy) + '\n';
                }

                return str;
            }



            std::string ToString() {
                std::string str;

                str += m_name + '\t';
                str += std::to_string(m_taxa.size());
                for (auto [id, taxon] : m_taxa) {
                    str += taxon.ToString() + '\n';
                }

                return str;
            }

            TaxonMap& GetTaxa() {
                return m_taxa;
            }

            const TaxonMap& GetTaxa() const {
                return m_taxa;
            }

            void ClearSams() {
                for (auto& [id, _] : m_taxa) m_taxa.at(id).ClearSams();
            }

            void AnnotateWithTruth(TruthSet const& set, TaxonFilterObj& filter, std::string& output) {
                std::ofstream os(output, std::ios::out);

//                std::cout << "Truth: ________________" << std::endl;

                for (auto& [key, taxon] : m_taxa) {
                    bool positive = set.contains(key);
                    bool prediction = filter.Pass(taxon);

//                    if (!positive && !prediction) continue;
                    std::vector<size_t> alleles = taxon.GetAlleles();
                    std::vector<size_t> filtered_alleles = taxon.GetAlleles(2, 60);
                    alleles.resize(5, 0);
                    filtered_alleles.resize(5, 0);


                    os << positive << "\t"; //1
                    os << prediction << "\t";
                    os << key << "\t";
                    os << taxon.PresentGenes() << "\t";
                    os << taxon.TotalHits() << "\t";
                    os << taxon.UniqueHits() << "\t";
                    os << taxon.GetMeanANI() << "\t";
                    os << TaxonFilter::ExpectedGenePresence(taxon) << '\t';
                    os << TaxonFilter::ExpectedGenePresenceRatio(taxon) << '\t';
                    os << taxon.Uniqueness() << '\t'; //10
                    os << taxon.GetMeanMAPQ() << '\t';
                    os << taxon.GetGeneVariance(1) << '\t';
                    os << taxon.GetGeneVariance(5) << '\t';
                    for (auto i = 0; i < 5; i++) {
                        os << alleles[i] << '\t';
                    }
                    for (auto i = 0; i < 5; i++) {
                        os << filtered_alleles[i] << '\t';
                    }
                    os << taxon.VCovStdDev() << '\t'; // 24
                    os << taxon.GetGenomeGeneNumber() << '\t';
                    os << taxon.UniqueKmers() << '\t'; //26
                    os << taxon.GenesWithUniqueKmers() << '\t';
                    os << taxon.UniqueTwoKmers() << '\t';
                    os << taxon.GenesWithUniqueTwoKmers() << '\t'; // 29

                    auto [su, lu, lsu, all] = m_genome_loader.GetGenome(key).GetUniqueKmerCounts();

                    // Defined in index for each species (genome)
                    os << su << '\t';
                    os << lu << '\t';
                    os << lsu << '\t';
                    os << all << '\t';

                    os << lu/static_cast<double>(taxon.UniqueKmers()) << '\t';
                    os << lsu/static_cast<double>(taxon.GenesWithUniqueTwoKmers());
                    os << std::endl;

//
//                    if (prediction && taxon.GetMeanANI() < 0.95) {
//                        exit(19);
//                    }
                }

                os.close();
            }

            void WriteSparseProfile(taxonomy::IntTaxonomy& taxonomy, TaxonFilterObj const& filter, std::ostream &os_filtered=std::cout, std::ostream* os_total=nullptr, std::ostream* os_dismissed=nullptr) {
                bool one_pass = false;

                std::string gene_covs_str = "";
                std::string gene_cov_ratios_str = "";
                std::vector<size_t> gene_covs(120, 0);
                std::vector<double> gene_cov_ratios(120, 0.0);


                double total_vcov = 0;
                for (auto& pair : m_taxa) {
                    auto& taxon = m_taxa.at(pair.first);

                    pair.second.TotalHits();
                    bool prediction = filter.Pass(taxon);
                    if (!prediction) continue;
                    total_vcov += taxon.VerticalCoverage();
                }

                for (auto& [key, _] : m_taxa) {
                    // this is necessary as taxon cannot be constant
                    auto& taxon = m_taxa.at(key);
                    bool prediction = filter.Pass(taxon);
                    one_pass |= prediction;

                    gene_covs_str.clear();
                    gene_cov_ratios_str.clear();
                    std::fill(gene_covs.begin(), gene_covs.end(), 0);
                    std::fill(gene_cov_ratios.begin(), gene_cov_ratios.end(), 0.0);

                    std::vector<size_t> alleles = taxon.GetAlleles();
                    std::vector<size_t> filtered_alleles = taxon.GetAlleles(2, 60);
                    alleles.resize(5, 0);
                    filtered_alleles.resize(5, 0);

                    for (auto& [id, _] : taxon.GetGenes()) {
                        auto& gene = taxon.GetGenes().at(id);
                        size_t cov = gene.GetStrainLevel().GetSequenceRangeHandler().CoveredPortion();
                        double ratio = static_cast<double>(cov) / gene.m_gene_length;
                        gene_covs[id] = cov;
                        gene_cov_ratios[id] = ratio;

                        auto stats_line = gene.GetStatisticsString();
                        *os_dismissed << key << '\t';
                        *os_dismissed << taxonomy.LineageStr(key) << '\t';
                        *os_dismissed << taxon.GetName() << '\t';
                        *os_dismissed << id << '\t';
                        *os_dismissed << stats_line << '\n';
                    }


                    gene_covs_str = std::accumulate(gene_covs.begin(), gene_covs.end(), std::string{}, [](std::string acc, size_t val) {
                        return(acc + "\t" + std::to_string(val));
                    });
                    gene_cov_ratios_str = std::accumulate(gene_cov_ratios.begin(), gene_cov_ratios.end(), std::string{}, [](std::string acc, double val) {
                        return(acc + "\t" + std::to_string(val));
                    });

                    double mean_gene_covs = std::accumulate(gene_covs.begin(), gene_covs.end(), 0) /
                                            static_cast<double>(taxon.GetGenes().size());
                    double mean_gene_cov_ratios = std::accumulate(gene_cov_ratios.begin(), gene_cov_ratios.end(), 0) /
                                            static_cast<double>(taxon.GetGenes().size());

                    auto node = taxonomy.Get(key);

                    if (prediction) {
                        os_filtered << node.rep_genome << '\t' << "d__Bacteria|" << taxonomy.LineageStr(key) << '\t' << taxon.GetAbundance(total_vcov) << std::endl;
                    }


                    // Print all outputs
                    if (os_total) {
                        *os_total << (prediction ? "1" : "0") << "\t" << node.rep_genome << '\t' << "d__Bacteria|" << taxonomy.LineageStr(key) << '\t' << (prediction ? taxon.GetAbundance(total_vcov) : 0);
                        *os_total << '\t' << taxon.VCovStdDev();
                        *os_total << '\t' << taxon.GetGeneVariance();
                        *os_total << '\t' << taxon.GetGeneVariance(5);
                        *os_total << '\t' << taxon.ToString(taxonomy);
                        *os_total << '\t' << taxon.VerticalCoverage();
                        *os_total << '\t' << mean_gene_covs << '\t' << mean_gene_cov_ratios << '\t' << gene_covs_str;
                        *os_total << '\t' << gene_cov_ratios_str << std::endl;
                    }
                }


                if (!one_pass) {
                    std::cout << "No taxon passes criteria" << std::endl;
                }
            }


        private:
            std::string m_name;
            mutable TaxonMap m_taxa;
            GenomeLoader &m_genome_loader;
        };


        using SamPairList = std::vector<std::vector<AlignmentPair>>;
        using SamPairs = std::vector<AlignmentPair>;
        using SamPairList_ptr = std::vector<std::vector<AlignmentPair>*>;
        using SamPairs_ptr = std::vector<AlignmentPair*>;

        template<typename ScoringSystem=ScoreAlignments>
        class Profiler {

        public:
            Profiler(GenomeLoader& genome_loader) :
                    m_genome_loader(genome_loader) {}

        private:
            GenomeLoader& m_genome_loader;


            SamPairs_ptr m_pairs_unique_ptr;
            SamPairList_ptr m_pairs_nonunique_ptr;

            ScoringSystem m_score;

            CigarInfo m_info1;
            CigarInfo m_info2;

            size_t m_min_alignment_length = 50;
            size_t m_min_mapq = 4;

            bool m_no_strain = true;

        public:

            SamPairs m_pairs_unique;
            SamPairList m_pairs_nonunique;
            SamPairs m_pairs_nonunique_best;
            Benchmark m_post_process_bm{"Post-processing"};

            void SetNoStrain(bool val) {
                m_no_strain = val;
            }

            static std::pair<double, int> ScorePairedAlignment(OptIRA const& a, OptIRA const& b, double divide_penalty=2.2, double alone_penalty=2) {
                int score = 0;
                double ani = 0;
                if (a.has_value() && b.has_value()) {
                    score = (a.value().alignment_score + b.value().alignment_score) / divide_penalty;
                    ani = (a.value().alignment_ani + b.value().alignment_ani) / 2;
                } else if (a.has_value()) {
                    score = a.value().alignment_score - alone_penalty;
                    ani = a.value().alignment_ani;
                } else if (b.has_value()) {
                    score = b.value().alignment_score - alone_penalty;
                    ani = b.value().alignment_ani;
                }
                return { ani, score };
            }

            static std::pair<double, int> ScorePairedAlignment(OptIRAPair const& pair, double divide_penalty=2.2, double alone_penalty=2) {
                return ScorePairedAlignment(pair.first, pair.second);
            }

            static bool SameReadId(OptIRA const &ira1, OptIRA const &ira2, OptIRAPair const& pair) {
                assert(pair.first.has_value() || pair.second.has_value());
                assert(ira1.has_value() || ira2.has_value());

                if (ira1.has_value() && ira2.has_value() && ira1.value().readid != ira2.value().readid) {
                    exit(12);
                }

                auto& rid = ira1.has_value() ? ira1.value().readid : ira2.value().readid;
                auto& rid2 = pair.first.has_value() ? pair.first.value().readid : pair.second.value().readid;

                return rid == rid2;
            }

            static bool SameRead(AlignmentPair& pair, AlignmentPair& other) {
                std::string& qname1 = pair.Any().m_qname;
                std::string& qname2 = other.Any().m_qname;

                std::string_view pairless_header1(qname1.c_str(), qname1.length() - 1);
                std::string_view pairless_header2(qname2.c_str(), qname2.length() - 1);

                return pairless_header1 == pairless_header2;
            }

//            void Process(std::optional<InternalReadAlignment> const &ira,
//                         std::optional<InternalReadAlignment> const &ira2,
//                         bool unique) {
//                assert(ira.has_value() || ira2.has_value());
//                if (unique) {
//                    ProcessUnique(ira, ira2);
//                } else {
////                    ProcessNonUnique(ira, ira2);
//                }
//            }

//            void Process(size_t read_id, SamEntry &sam1, SamEntry &sam2, bool has_sam1, bool has_sam2, bool unique) {
//                CigarInfo info;
//                size_t min_alignment_length = 70;
//                if (has_sam1) {
//                    CompressedCigarInfo(sam1.m_cigar, info);
//                    if (info.clipped_alignment_length < min_alignment_length) has_sam1 = false;
//                }
//                if (has_sam2) {
//                    CompressedCigarInfo(sam2.m_cigar, info);
//                    if (info.clipped_alignment_length < min_alignment_length) has_sam2 = false;
//                }
//                if (!has_sam1 && !has_sam2) return;
//                OptIRA oira1 = has_sam1 ? OptIRA(InternalReadAlignment(read_id, sam1)) : OptIRA();
//                OptIRA oira2 = has_sam2 ? OptIRA(InternalReadAlignment(read_id, sam2)) : OptIRA();
//                Process(oira1, oira2, unique);
//            }

            void PrintStats() {
                std::cout << "Unique Sams: " << m_pairs_unique.size() << std::endl;
                std::cout << "NonUnique Sams: " << m_pairs_nonunique.size() << std::endl;
                std::cout << "NonUnique Best Sams: " << m_pairs_nonunique_best.size() << std::endl;
            }

            bool HasReads() const {
                return !(m_pairs_unique.empty() && m_pairs_nonunique.empty() && m_pairs_nonunique_best.empty());
            }

            void FromSam(std::string file_path, bool truth_in_header=false) {
                m_pairs_unique.clear();
                m_pairs_nonunique.clear();
                m_pairs_nonunique_best.clear();
//                std::ifstream file(file_path, std::ios::in);
                igzstream file{ file_path.c_str() };

                std::string delim = "\t";
                std::string line;
                std::vector<std::string> tokens;
                SamEntry current;
                SamEntry current_other;
                SamEntry next;
                SamEntry next_other;

                bool has_current1 = false, has_current2 = false, has_next1 = false, has_next2 = false;

                std::string current_qname = "__";
                std::string next_qname = "__";
                std::string last_qname = "__";

                if (!std::getline(file, line)) {
                    std::cerr << "sam file is empty " << file_path << std::endl;
                    exit(9);
                }

                auto count = 0;
                while (line[0] == '@' && !file.eof()) {
                    std::getline(file, line);
                }

                if (file.eof()) {
                    std::cerr << "sam file contains no sequences (header only)\n\t" << file_path << std::endl;
                    return;
                }

                GetSamPair(file, line, tokens, current, current_other, has_current1, has_current2, true);

                std::vector<AlignmentPair> pair_list;

                size_t read_id = 0;
                size_t count_read_lines = 0;
                while (GetSamPair(file, line, tokens, next, next_other, has_next1, has_next2)) {
                    current_qname = has_current1 ? current.m_qname : current_other.m_qname;
                    next_qname = has_next1 ? next.m_qname : next_other.m_qname;
                    bool unique = current_qname != last_qname && current_qname != next_qname;
                    last_qname = current_qname;

                    AlignmentPair new_pair = AlignmentPair(
                        has_current1 ? std::optional<SamEntry>{current} : std::optional<SamEntry>{},
                        has_current2 ? std::optional<SamEntry>{current_other} : std::optional<SamEntry>{}
                    );

                    if (unique) {
                        m_pairs_unique.emplace_back(std::move(new_pair));
                    } else {
                        if (!pair_list.empty() && !SameRead(new_pair, pair_list.front())) {
                            m_pairs_nonunique_best.emplace_back(pair_list.front());
                            m_pairs_nonunique.emplace_back(std::move(pair_list));
                            pair_list.clear();
                        }
                        pair_list.emplace_back(new_pair);
                    }

                    std::swap(current, next);
                    std::swap(current_other, next_other);
                    has_current1 = has_next1;
                    has_current2 = has_next2;
                    count_read_lines++;
                }
                if (count_read_lines == 0) {
                    std::cerr << "File: " << file_path << " does not contain any alignments." << std::endl;
                }


                AlignmentPair new_pair = AlignmentPair(
                        has_current1 ? std::optional<SamEntry>{current} : std::optional<SamEntry>{},
                        has_current2 ? std::optional<SamEntry>{current_other} : std::optional<SamEntry>{}
                );
                current_qname = has_current1 ? current.m_qname : current_other.m_qname;

                bool unique = current_qname != last_qname;
                if (unique) {
                    m_pairs_unique.emplace_back(std::move(new_pair));
                } else {
                    if (!pair_list.empty() && !SameRead(new_pair, pair_list.front())) {
                        m_pairs_nonunique_best.emplace_back(pair_list.front());
                        m_pairs_nonunique.emplace_back(std::move(pair_list));
                        pair_list.clear();
                    }
                    pair_list.emplace_back(new_pair);
                }

//                std::cout << "Uniques:          " << m_pairs_unique.size() << std::endl;
//                std::cout << "Non-Uniques:      " << m_pairs_nonunique.size() << std::endl;
//                std::cout << "Non-Uniques best: " << m_pairs_nonunique_best.size() << std::endl;

                if (truth_in_header) {
                    OutputErrorData(m_pairs_unique, m_pairs_nonunique);
                }
            }


            std::string ErrorLine(int id, int length1, int score1, bool forward1, int length2, int score2, bool forward2, bool possibly_true, bool same_gene, int distance, bool paired, int ref_length) {
                std::string line;
                line += std::to_string(id) + '\t';
                line += std::to_string(length1) + '\t';
                line += std::to_string(score1) + '\t';
                line += std::to_string(forward1) + '\t';
                line += std::to_string(length2) + '\t';
                line += std::to_string(score2) + '\t';
                line += std::to_string(forward2) + '\t';
                line += std::to_string(possibly_true) + '\t';
                line += std::to_string(same_gene) + '\t';
                line += std::to_string(distance) + '\t';
                line += std::to_string(paired) = '\t';
                line += std::to_string(ref_length);
                return line;
            }

            void TestSNPUtils(std::vector<AlignmentPair> &pairs) {
//                std::cout << "TestSNPUtils" << std::endl;
                SNPList snps;
                size_t sum = 0;
                Benchmark bm("SNPs");
                for (auto &pair: pairs) {
                    auto &any = pair.Any();
                    auto &[tid, gid] = ExtractTaxidGeneid(any.m_rname);
                    auto &gene = m_genome_loader.GetGenome(tid).GetGeneOMP(gid);

                    gene.LoadOMP();
                    auto &ref = gene.Sequence();
//                        std::cout << "Extract: " << any.m_qname << " ---> " <<  any.m_rname << std::endl;

                    PrintAlignment(any, ref);
                    ExtractSNPs(any, ref, snps, tid, gid);
//                    for (auto &snp: snps) {
//                        std::cout << snp.ToString() << std::endl;
//                    }
//                    sum += snps.size();
//                    std::cout << ".........>" << std::endl;
//                    Utils::Input();
                }
                bm.PrintResults();
                std::cout << "Sum: " << sum << std::endl;

            }

            void TestSNPUtils(std::vector<std::vector<AlignmentPair>> &pair_lists) {
//                std::cout << "TestSNPUtils" << std::endl;
                SNPList snps;
                size_t sum = 0;
                Benchmark bm("SNPs");
                for (auto& pairlist : pair_lists) {
                    for (auto &pair: pairlist) {
                        auto& any = pair.Any();
                        auto &[ tid, gid ] = ExtractTaxidGeneid(any.m_rname);
                        auto& gene = m_genome_loader.GetGenome(tid).GetGene(gid);
                        gene.LoadOMP();
                        auto& ref = gene.Sequence();
//                        std::cout << "Extract: " << any.m_qname << " ---> " <<  any.m_rname << std::endl;

                        PrintAlignment(any, ref);
                        ExtractSNPs(any, ref, snps, tid, gid);
                        for (auto& snp : snps) {
                            std::cout << snp.ToString() << std::endl;
                        }
                        sum += snps.size();
                    }
                }
                bm.PrintResults();
                std::cout << "Sum: " << sum << std::endl;

            }

            void OutputErrorData(std::vector<std::vector<AlignmentPair>> &pair_lists) {
                int pairlist_id = 0;
                CigarInfo info1;
                CigarInfo info2;

                size_t count_truth = 0;
                size_t count_false = 0;
                for (auto& pairlist : pair_lists) {
                    for (auto& pair : pairlist) {
                        auto &[ tid, gid ] = ExtractTaxidGeneid(pair.Any().m_rname);
                        auto &[ atid, agid ] = ExtractTaxidGeneid(pair.Any().m_qname);
//                        std::cerr << pair.Any().m_qname << " " << atid << " " << agid << std::endl;
//                        std::cerr << "CAPTURE\t";
                        if (pair.IsPair()) {
                            auto &[ tid1, gid1 ] = ExtractTaxidGeneid(pair.First().m_rname);
                            auto &[ tid2, gid2 ] = ExtractTaxidGeneid(pair.Second().m_rname);
                            auto &[ ttid1, tgid1 ] = ExtractTaxidGeneid(pair.First().m_qname);
                            auto &[ ttid2, tgid2 ] = ExtractTaxidGeneid(pair.Second().m_qname);
                            CompressedCigarInfo(pair.First().m_cigar, info1);
                            CompressedCigarInfo(pair.Second().m_cigar, info2);
                            bool same_gene = gid1 == gid2;
                            int insert = std::max(pair.First().m_pos, pair.Second().m_pos) - std::min(pair.First().m_pos, pair.Second().m_pos);
                            count_truth += tid == ttid1;
                            count_false += tid != ttid1;
//                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag),  tid == ttid1, same_gene, insert, true, m_genome_loader.GetGenome(tid1).GetGene(gid1).GetLength()) << std::endl;
                        } else if (pair.HasFirst()) {
                            CompressedCigarInfo(pair.First().m_cigar, info1);
                            auto &[ ttid, tgid ] = ExtractTaxidGeneid(pair.First().m_qname);
                            count_truth += tid == ttid;
                            count_false += tid != ttid;
//                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), -1, -1, false, tid == ttid, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
                        } else if (pair.HasSecond()) {
                            CompressedCigarInfo(pair.Second().m_cigar, info2);
                            auto &[ ttid, tgid ] = ExtractTaxidGeneid(pair.Second().m_qname);
                            count_truth += tid == ttid;
                            count_false += tid != ttid;
//                            std::cerr << ErrorLine(pairlist_id, -1, -1, false, info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag), tid == ttid, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
                        }
                    }

                    pairlist_id++;
                }
                std::cout << "Truth: " << count_truth << std::endl;
                std::cout << "False: " << count_false << std::endl;
            }

//            void OutputErrorData(TruthSet const& set, std::vector<std::vector<AlignmentPair>> &pair_lists) {
//                int pairlist_id = 0;
//                CigarInfo info1;
//                CigarInfo info2;
//                for (auto& pairlist : pair_lists) {
//                    for (auto& pair : pairlist) {
//                        auto &[ tid, gid ] = ExtractTaxidGeneid(pair.Any().m_rname);
//                        bool maybe_true = set.contains(tid);
//                        std::cerr << "CAPTURE\t";
//                        if (pair.IsPair()) {
//                            auto &[ tid1, gid1 ] = ExtractTaxidGeneid(pair.First().m_rname);
//                            auto &[ tid2, gid2 ] = ExtractTaxidGeneid(pair.Second().m_rname);
//                            CompressedCigarInfo(pair.First().m_cigar, info1);
//                            CompressedCigarInfo(pair.Second().m_cigar, info2);
//                            bool same_gene = gid1 == gid2;
//                            int insert = std::max(pair.First().m_pos, pair.Second().m_pos) - std::min(pair.First().m_pos, pair.Second().m_pos);
//
//                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag),  maybe_true, same_gene, insert, true, m_genome_loader.GetGenome(tid1).GetGene(gid1).GetLength()) << std::endl;
//                        } else if (pair.HasFirst()) {
//                            CompressedCigarInfo(pair.First().m_cigar, info1);
//                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), -1, -1, false, maybe_true, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
//                        } else if (pair.HasSecond()) {
//                            CompressedCigarInfo(pair.Second().m_cigar, info2);
//                            std::cerr << ErrorLine(pairlist_id, -1, -1, false, info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag), maybe_true, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
//                        }
//                    }
//
//                    pairlist_id++;
//                }
//            }


            std::pair<bool,bool> IsSamCorrect(SamEntry &entry) {
                auto &[ rtid, rgid ] = ExtractTaxidGeneid(entry.m_rname);
                auto &[ qtid, qgid ] = ExtractTaxidGeneid(entry.m_rname);
                return { rtid == qtid, rtid == qtid && rgid == qgid };
            }

            std::pair<bool,bool> IsCorrect(SamEntry &entry, TruthSet* set=nullptr) {
                if (set != nullptr) {
                    auto &[ tid, gid ] = ExtractTaxidGeneid(entry.m_rname);
                    return { set->contains(tid), false };
                } else {
                    return IsSamCorrect(entry);
                }
            }

            void OutputErrorData(SamPairs &pairs, SamPairList &pair_lists, TruthSet* set=nullptr) {
                int pairlist_id = 0;
                CigarInfo info1;
                CigarInfo info2;

                size_t count_truth = 0;
                size_t count_false = 0;
                for (auto& pairlist : pair_lists) {
                    for (auto& pair : pairlist) {
//                        std::cerr << "CAPTURE\t";
                        if (pair.IsPair()) {
                            auto &[ tid1, gid1 ] = ExtractTaxidGeneid(pair.First().m_rname);
                            auto &[ tid2, gid2 ] = ExtractTaxidGeneid(pair.Second().m_rname);
                            auto [true_taxid1, true_geneid1] = IsCorrect(pair.First(), set);
                            auto [true_taxid2, true_geneid2] = IsCorrect(pair.Second(), set);
                            auto &[read_truth_tid, read_truth_gid] = ExtractTaxidGeneid(pair.First().m_qname);
                            bool same_gene = gid1 == gid2;
                            CompressedCigarInfo(pair.First().m_cigar, info1);
                            CompressedCigarInfo(pair.Second().m_cigar, info2);

                            count_truth += read_truth_tid == tid1;
                            count_false += read_truth_tid != tid1;

                            int insert = std::max(pair.First().m_pos, pair.Second().m_pos) - std::min(pair.First().m_pos, pair.Second().m_pos);
//                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag), true_taxid1, same_gene, insert, true, m_genome_loader.GetGenome(tid1).GetGene(gid1).GetLength()) << std::endl;
                        } else if (pair.HasFirst()) {
                            auto &[ tid, gid ] = ExtractTaxidGeneid(pair.First().m_rname);
                            auto [true_taxid, true_geneid] = IsCorrect(pair.First(), set);
                            auto &[read_truth_tid, read_truth_gid] = ExtractTaxidGeneid(pair.First().m_qname);
                            CompressedCigarInfo(pair.First().m_cigar, info1);
                            count_truth += read_truth_tid == tid;
                            count_false += read_truth_tid != tid;
//                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), -1, -1, false, true_taxid, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
                        } else if (pair.HasSecond()) {
                            auto &[ tid, gid ] = ExtractTaxidGeneid(pair.Second().m_rname);
                            auto [true_taxid, true_geneid] = IsCorrect(pair.Second(), set);
                            auto &[read_truth_tid, read_truth_gid] = ExtractTaxidGeneid(pair.Second().m_qname);
                            CompressedCigarInfo(pair.Second().m_cigar, info2);
                            count_truth += read_truth_tid == tid;
                            count_false += read_truth_tid != tid;

//                            std::cerr << ErrorLine(pairlist_id, -1, -1, false, info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag), true_taxid, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
                        }
                        break;
                    }

                    pairlist_id++;
                }
                for (auto &pair: pairs) {
//                    std::cerr << "CAPTURE\t";
                    if (pair.IsPair()) {
                        auto &[tid1, gid1] = ExtractTaxidGeneid(pair.First().m_rname);
                        auto &[tid2, gid2] = ExtractTaxidGeneid(pair.Second().m_rname);
                        auto &[read_truth_tid, read_truth_gid] = ExtractTaxidGeneid(pair.First().m_qname);
                        auto [true_taxid1, true_geneid1] = IsCorrect(pair.First(), set);
                        auto [true_taxid2, true_geneid2] = IsCorrect(pair.Second(), set);
                        bool same_gene = gid1 == gid2;

                        std::cout << tid1 << " " << read_truth_tid << std::endl;

                        CompressedCigarInfo(pair.First().m_cigar, info1);
                        CompressedCigarInfo(pair.Second().m_cigar, info2);
                        count_truth += read_truth_tid == tid1;
                        count_false += read_truth_tid != tid1;
                        int insert = std::max(pair.First().m_pos, pair.Second().m_pos) -
                                     std::min(pair.First().m_pos, pair.Second().m_pos);
//                        std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(),
//                                               Flag::IsRead1ReverseComplement(pair.First().m_flag),
//                                               info2.clipped_alignment_length, info2.Score(),
//                                               Flag::IsRead2ReverseComplement(pair.Second().m_flag), true_taxid1,
//                                               same_gene, insert, true,
//                                               m_genome_loader.GetGenome(tid1).GetGene(gid1).GetLength()) << std::endl;
                    } else if (pair.HasFirst()) {
                        auto [true_taxid, true_geneid] = IsCorrect(pair.First(), set);
                        auto &[ tid, gid ] = ExtractTaxidGeneid(pair.First().m_rname);
                        auto &[read_truth_tid, read_truth_gid] = ExtractTaxidGeneid(pair.First().m_qname);
                        CompressedCigarInfo(pair.First().m_cigar, info1);
                        count_truth += read_truth_tid == tid;
                        count_false += read_truth_tid != tid;
//                        std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(),
//                                               Flag::IsRead1ReverseComplement(pair.First().m_flag), -1, -1, false,
//                                               true_taxid, false, -1, false,
//                                               m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
                    } else if (pair.HasSecond()) {
                        auto &[ tid, gid ] = ExtractTaxidGeneid(pair.Second().m_rname);
                        auto [true_taxid, true_geneid] = IsCorrect(pair.Second(), set);
                        auto &[read_truth_tid, read_truth_gid] = ExtractTaxidGeneid(pair.Second().m_qname);
                        CompressedCigarInfo(pair.Second().m_cigar, info2);
                        count_truth += read_truth_tid == tid;
                        count_false += read_truth_tid != tid;
//                        std::cerr
//                                << ErrorLine(pairlist_id, -1, -1, false, info2.clipped_alignment_length, info2.Score(),
//                                             Flag::IsRead2ReverseComplement(pair.Second().m_flag), true_taxid, false,
//                                             -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength())
//                                << std::endl;
                    }
                }


                std::cout << "Truth: " << count_truth << std::endl;
                std::cout << "False: " << count_false << std::endl;

                pairlist_id++;
            }

            Benchmark m_add_sam{"Add sam.."};
            Benchmark m_cigar_info{"Compressed cigar info"};


            bool ProcessMAPQ(MicrobialProfile& profile, AlignmentPair& ap, int read_id=0) {
                bool valid_sam = true;

                bool take_first = false;
                bool take_second = false;
                m_cigar_info.Start();
                if (ap.HasFirst()) CompressedCigarInfo(ap.First().m_cigar, m_info1);
                if (ap.HasSecond()) CompressedCigarInfo(ap.Second().m_cigar, m_info2);
                m_cigar_info.Stop();

                if ((ap.HasFirst() && ap.First().m_mapq < m_min_mapq) || (ap.HasSecond() && ap.Second().m_mapq < m_min_mapq)) {
                    return true;
                }

                if (ap.HasFirst() && m_info1.clipped_alignment_length > m_min_alignment_length) {
                    take_first = true;
                }
                if (ap.HasSecond() && m_info2.clipped_alignment_length > m_min_alignment_length) {
                    take_second = true;
                }
                if (!take_first && !take_second) return true;


                if (take_first && take_second) {
                    auto [tid1, geneid1] = ExtractTaxidGeneid(ap.First().m_rname);
                    auto [tid2, geneid2] = ExtractTaxidGeneid(ap.Second().m_rname);

                    double score = m_score.Score(ap.First(), ap.Second());
                    double score1 = m_score.Score(ap.First());
                    double score2 = m_score.Score(ap.Second());
                    double score_diff = std::abs(score1 - score2);

                    m_add_sam.Start();
                    valid_sam &= profile.AddSam(tid1, geneid1, ap.First(), m_info1.Ani(), true, read_id, m_no_strain);
                    valid_sam &= profile.AddSam(tid2, geneid2, ap.Second(), m_info2.Ani(), true, read_id, m_no_strain);
                    m_add_sam.Stop();

                } else if (take_first) {
                    auto [tid, geneid] = ExtractTaxidGeneid(ap.First().m_rname);
                    double score = m_score.Score(ap.First());

                    m_add_sam.Start();
                    valid_sam &=profile.AddSam(tid, geneid, ap.First(), m_info1.Ani(), true, read_id);
                    m_add_sam.Stop();



                } else if (take_second) {
                    auto [tid, geneid] = ExtractTaxidGeneid(ap.Second().m_rname);
                    double score = m_score.Score(ap.Second());

                    m_add_sam.Start();
                    valid_sam &=profile.AddSam(tid, geneid, ap.Second(), m_info2.Ani(), true, read_id);
                    m_add_sam.Stop();
                }
                return valid_sam;
            }

            void ProcessUnique(MicrobialProfile& profile, AlignmentPair& ap) {
                bool take_first = false;
                bool take_second = false;
                if (ap.HasFirst()) CompressedCigarInfo(ap.First().m_cigar, m_info1);
                if (ap.HasSecond()) CompressedCigarInfo(ap.Second().m_cigar, m_info2);

                if (ap.HasFirst() && ap.HasSecond() && (m_info1.clipped_alignment_length + m_info2.clipped_alignment_length) > m_min_alignment_length) {
                    take_first = true;
                    take_second = true;
                } else if (ap.HasFirst() && m_info1.clipped_alignment_length > m_min_alignment_length) {
                    take_first = true;
                } else if (ap.HasSecond() && m_info2.clipped_alignment_length > m_min_alignment_length) {
                    take_second = true;
                }

                if (!take_first && !take_second) return;

                if (take_first && take_second) {
                    auto [tid1, geneid1] = ExtractTaxidGeneid(ap.First().m_rname);
                    auto [tid2, geneid2] = ExtractTaxidGeneid(ap.Second().m_rname);

                    double score = m_score.Score(ap.First(), ap.Second());
                    double score1 = m_score.Score(ap.First());
                    double score2 = m_score.Score(ap.Second());
                    double score_diff = std::abs(score1 - score2);

                    profile.AddSam(tid1, geneid1, ap.First(), score, true);
                    profile.AddSam(tid2, geneid2, ap.Second(), score, true);

                } else if (take_first) {
                    auto [tid, geneid] = ExtractTaxidGeneid(ap.First().m_rname);
                    double score = m_score.Score(ap.First());

                    profile.AddSam(tid, geneid, ap.First(), score, true);
                } else if (take_second) {
                    auto [tid, geneid] = ExtractTaxidGeneid(ap.Second().m_rname);
                    double score = m_score.Score(ap.Second());

                    profile.AddSam(tid, geneid, ap.Second(), score, true);
                }
            }

            using OptionalRefOstream = std::optional<std::reference_wrapper<std::ostream>>;
            MicrobialProfile Profile(OptionalRefOstream erroneous_sam_out={}) {
                MicrobialProfile profile(m_genome_loader);

                size_t read_id = 0;

                Benchmark bm_process{"Process sams"};
                bm_process.Start();
                for (auto& sam_pair : m_pairs_unique) {
                    bool valid = ProcessMAPQ(profile, sam_pair, read_id++);

                    if (!valid & erroneous_sam_out.has_value()) {
                        auto& os = erroneous_sam_out.value().get();
#pragma omp critical (err_sam)
                        {
                            if (sam_pair.first.has_value()) {
                                os << sam_pair.first.value().ToString() << std::endl;
                            }
                            if (sam_pair.second.has_value()) {
                                os << sam_pair.second.value().ToString() << std::endl;
                            }
                        }
                    }
                }
                bm_process.Stop();
                bm_process.PrintResults();

                // m_add_sam.PrintResults();
                // m_cigar_info.PrintResults();
                // std::cout << "Processed all sams" << std::endl;

                for (auto sam_pair : m_pairs_nonunique_best) {
                    ProcessMAPQ(profile, sam_pair, read_id++);
                }

                m_post_process_bm.Start();
                profile.PostProcessSNPs();
                m_post_process_bm.Stop();


                profile.ClearSams();


                return profile;
            }
        };
    }
}