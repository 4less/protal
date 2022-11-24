//
// Created by fritsche on 07/10/22.
//

#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <sparse_map.h>
#include "LineSplitter.h"
#include "Taxonomy.h"
#include "InternalReadAlignment.h"
#include "Constants.h"
#include "AlignmentUtils.h"
#include "SNPUtils.h"

namespace protal {


    using TruthSet = tsl::robin_set<uint32_t>;
    TruthSet GetTruth(std::string const& file_path) {
        std::ifstream is(file_path, std::ios::in);
        TruthSet truths;
        std::vector<std::string> tokens;
        std::string delim = "\t";

        std::string line;
        while (std::getline(is, line)) {
            LineSplitter::Split(line, delim, tokens);
            uint32_t truth = std::stoul(tokens[0]);
            truths.insert(truth);
        }
        return truths;
    }

    namespace profiler {
        struct MicrobialProfile;
        struct Gene;
        struct Taxon;

        using GeneMap = tsl::sparse_map<uint32_t, Gene>;
        using TaxonMap = tsl::sparse_map<uint32_t, Taxon>;

        // Data structures

        class Gene {
        public:
            size_t m_mapped_reads = 0;
            size_t m_gene_length = 0;

            void AddRead(GenePos pos) {
                m_mapped_reads++;
            }

            void SetLength(size_t length) {
                m_gene_length = length;
            }

            double VerticalCoverage() const {
//                std::cout << static_cast<double>(m_mapped_reads * 150) << "/" << static_cast<double>(m_gene_length) << std::endl;
                return static_cast<double>(m_mapped_reads * 150)/static_cast<double>(m_gene_length);
            }
        };

        class Taxon {
        private:
            std::string m_name;
            size_t m_total_hits = 0;
            size_t m_unique_hits = 0;
            double m_ani_sum = 0;
            double m_vcov = -1;
            GeneMap m_genes;


            Genome m_genome;
            size_t m_genome_gene_count = 0;

        public:
            Taxon(Genome& genome) : m_genome(genome), m_genome_gene_count(genome.GeneNum()) {}

            void AddHit(GeneId geneid, GenePos genepos, double ani, bool unique) {
                if (!m_genes.contains(geneid)) {
                    if (!m_genome.GetGene(geneid).IsLoaded()) {
                        m_genome.GetGene(geneid).Load();
                    }
                    m_genes[geneid] = Gene();
                    m_genes[geneid].SetLength(m_genome.GetGene(geneid).Sequence().length());
                }

                m_genes[geneid].AddRead(genepos);
                m_total_hits++;
                m_ani_sum += ani;
                m_unique_hits += unique;
            }

            size_t GetGenomeGeneNumber() const {
                return m_genome.GeneNum();
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

            size_t TotalHits() const {
                return m_total_hits;
            }

            size_t UniqueHits() const {
                return m_unique_hits;
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

            double VerticalCoverage(bool force=false) {
                if (m_vcov == -1 || force) {
                    std::vector<double> vcovs;
                    for (auto& [geneid, gene] : m_genes) {
                        vcovs.emplace_back(gene.VerticalCoverage());
                    }
                    std::sort(vcovs.begin(), vcovs.end());

                    if (vcovs.empty()) return 0.0;

                    size_t mid_index = vcovs.size()/2;
                    m_vcov = vcovs.size() % 2 == 1 ?
                                  vcovs.at(mid_index) :
                                  (vcovs.at(mid_index - 1) + vcovs.at(mid_index)) / 2;
                }


                return m_vcov;
            }

            double VerticalCoverage() const {
                return m_vcov;
            }

            double GetAbundance(double total_vertical_coverage) {
                return VerticalCoverage()/total_vertical_coverage;
            }

            std::string ToString() const {
                std::string str;

                str += m_name + "{ ";
                str += "Genes: " + std::to_string(m_genes.size()) + ", ";
                str += "Total Hits: " + std::to_string(m_total_hits) + ", ";
                str += "MeanANI: " + std::to_string(m_ani_sum/m_total_hits);
                str += " }";

                return str;
            }
        };

        class TaxonFilter {
            double m_min_mean_ani = 0.94;
            double m_min_gene_presence = 0.7;
            size_t m_min_reads = 50;
            double m_min_uniqueness = 0.1;
        public:
            TaxonFilter() {}
            TaxonFilter(double min_mean_ani, double min_gene_presence, double min_reads) :
                    m_min_mean_ani(min_mean_ani),
                    m_min_gene_presence(min_gene_presence),
                    m_min_reads(min_reads) {
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
                        taxon.GetMeanANI() > m_min_mean_ani &&
                        taxon.TotalHits() > m_min_reads &&
                        ExpectedGenePresenceRatio(taxon) > m_min_gene_presence;
                return pass;
            }

            bool Pass(Taxon const& taxon) const {
                return Formula1(taxon);
            }
        };


        class MicrobialProfile {
        public:
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

            void SetName(std::string name) {
                m_name = name;
            }

            std::string ToString() {
                std::string str;

                str += m_name + '\t';
                str += std::to_string(m_taxa.size());

                return str;
            }

            TaxonMap& GetTaxa() {
                return m_taxa;
            }


            void AnnotateWithTruth(TruthSet const& set) {
//                std::cout << "Truth: ________________" << std::endl;
                TaxonFilter filter(0.95, 0.7, 70);
                for (auto& [key, taxon] : m_taxa) {
                    bool positive = set.contains(key);
                    bool prediction = filter.Pass(taxon);

//                    if (!positive && !prediction) continue;


                    std::cerr << "ANNOTATED\t" << positive << "\t";
                    std::cerr << key << "\t";
                    std::cerr << taxon.PresentGenes() << "\t";
                    std::cerr << taxon.TotalHits() << "\t";
                    std::cerr << taxon.UniqueHits() << "\t";
                    std::cerr << taxon.GetMeanANI() << "\t";
                    std::cerr << TaxonFilter::ExpectedGenePresence(taxon) << '\t';
                    std::cerr << TaxonFilter::ExpectedGenePresenceRatio(taxon) << '\t';
                    std::cerr << taxon.Uniqueness() << '\t';
                    std::cerr << taxon.VCovStdDev() << '\t';
                    std::cerr << prediction << std::endl;

                    if (prediction && taxon.GetMeanANI() < 0.95) {
                        exit(19);
                    }
                }
            }

            void WriteSparseProfile(taxonomy::IntTaxonomy& taxonomy, TaxonFilter const& filter, std::ostream &os=std::cout) {

                bool one_pass = false;

                double total_vcov = 0;
                for (auto& pair : m_taxa) {
                    auto& taxon = m_taxa.at(pair.first);

                    pair.second.TotalHits();
                    bool prediction = filter.Pass(taxon);
                    if (!prediction) continue;
//                    std::cout << "Taxon vertical coverage:  " << taxon.VerticalCoverage() << std::endl;
                    total_vcov += taxon.VerticalCoverage();
                }

//                std::cout << "total_vcov: " << total_vcov << std::endl;

                for (auto& [key, _] : m_taxa) {
                    auto& taxon = m_taxa.at(key);
                    bool prediction = filter.Pass(taxon);
                    if (!prediction) continue;

                    one_pass = true;

                    auto node = taxonomy.Get(key);
                    os << node.rep_genome << '\t' << "d__Bacteria|" << taxonomy.LineageStr(key) << '\t' << taxon.GetAbundance(total_vcov) << std::endl;
                }

                if (!one_pass) {
                    for (auto &[key, taxon] : m_taxa) {
                        bool prediction = filter.Pass(taxon);

                        std::cerr << key << "\t";
                        std::cerr << taxon.PresentGenes() << "\t";
                        std::cerr << taxon.TotalHits() << "\t";
                        std::cerr << taxon.UniqueHits() << "\t";
                        std::cerr << taxon.GetMeanANI() << "\t";
                        std::cerr << TaxonFilter::ExpectedGenePresence(taxon) << '\t';
                        std::cerr << TaxonFilter::ExpectedGenePresenceRatio(taxon) << '\t';
                        std::cerr << taxon.Uniqueness() << '\t';
                        std::cerr << taxon.VerticalCoverage() << '\t';
                        std::cerr << prediction << std::endl;
                    }
                }
            }


        private:
            std::string m_name;
            mutable TaxonMap m_taxa;
            GenomeLoader &m_genome_loader;
        };



        class Profiler {

        public:
            Profiler(GenomeLoader& genome_loader) :
                    m_genome_loader(genome_loader) {}

            enum ReadEvidenceQuality{
                EXCELLENT,
                GOOD,
                MEDIUM,
                BAD,
                SIZE
            };
            static const ReadEvidenceQuality quality_enum_list[];
            static const std::string quality_map[];


            ReadEvidenceQuality GetReadEvidenceQuality(double ani, int alignment_score) {
                if (ani >= 0.99) return ReadEvidenceQuality::EXCELLENT;
                if (ani >= 0.97) return ReadEvidenceQuality::GOOD;
                if (ani >= 0.95) return ReadEvidenceQuality::MEDIUM;
                return ReadEvidenceQuality::BAD;
            }

            class AlignmentLists {
            private:
                std::vector<InternalReadAlignmentList> m_unique_alignments;
                std::vector<InternalNonUniqueReadAlignmentListList> m_non_unique_alignments;
                std::vector<InternalNonUniqueReadOptIRAPairListList> m_non_unique_paired_alignments;
                size_t m_size = ReadEvidenceQuality::SIZE;

            public:
                AlignmentLists() {
                    for (int i = 0; i < m_size; i++)
                        m_unique_alignments.emplace_back(InternalReadAlignmentList{});
                    for (int i = 0; i < m_size; i++)
                        m_non_unique_alignments.emplace_back(InternalNonUniqueReadAlignmentListList{});
                    for (int i = 0; i < m_size; i++)
                        m_non_unique_paired_alignments.emplace_back(InternalNonUniqueReadOptIRAPairListList{});
                };

                InternalReadAlignmentList &GetUniqueList(ReadEvidenceQuality quality) {
                    return m_unique_alignments[quality];
                }
                InternalNonUniqueReadAlignmentListList &GetNonUniqueList(ReadEvidenceQuality quality) {
                    return m_non_unique_alignments[quality];
                }
                InternalNonUniqueReadOptIRAPairListList &GetNonUniquePairList(ReadEvidenceQuality quality) {
                    return m_non_unique_paired_alignments[quality];
                }
            };

        private:
            AlignmentLists m_alignments;
            InternalNonUniqueReadAlignmentList m_tmp;
            InternalNonUniqueReadOptIRAPairList m_tmp_pair;
            GenomeLoader& m_genome_loader;



        public:


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


            void
            ProcessUnique(std::optional<InternalReadAlignment> const &ira, std::optional<InternalReadAlignment> const &ira2, ReadEvidenceQuality min_qual = ReadEvidenceQuality::GOOD) {
                auto [ani, score] = ScorePairedAlignment(ira, ira2);
                auto qual = GetReadEvidenceQuality(ani, score);
                static int good_qual_counter = 0;
                if (ira.has_value()) m_alignments.GetUniqueList(qual).emplace_back(ira.value());
                if (ira2.has_value()) m_alignments.GetUniqueList(qual).emplace_back(ira2.value());
            }

            void ProcessNonUniqueWorker(std::vector<InternalReadAlignment> &non_uniques) {
                if (non_uniques.empty()) return;
                // Assume scores are
//                std::sort(non_uniques.begin(), non_uniques.end(),
//                          [](InternalReadAlignment const &a, InternalReadAlignment const &b) {
//                              return a.alignment_score > b.alignment_score;
//                          });
                auto &best = non_uniques.front();
                auto qual = GetReadEvidenceQuality(best.alignment_ani, best.alignment_score);

                m_alignments.GetNonUniqueList(qual).emplace_back(std::move(non_uniques));
                non_uniques.clear();
            }


            void ProcessNonUniqueWorker(InternalNonUniqueReadOptIRAPairList &non_uniques) {
                if (non_uniques.empty()) return;
                std::sort(non_uniques.begin(), non_uniques.end(), [](OptIRAPair const &a, OptIRAPair const &b) {
                    auto status_a = PairedStatus(a);
                    auto status_b = PairedStatus(b);

                    bool one_paired_one_single =
                            status_a == PAIRED && (status_b == FIRST_ONLY || status_b == SECOND_ONLY) ||
                            status_b == PAIRED && (status_a == FIRST_ONLY || status_a == SECOND_ONLY);
                    if (one_paired_one_single) {
                        auto pair = status_a == PAIRED ? a : b;
                        auto single = status_a == PAIRED ? b : a;

                        if (std::max(pair.first.value().alignment_score, pair.second.value().alignment_score) == (single.first.has_value() ? single.first.value().alignment_score : single.second.value().alignment_score)) {
                            return status_a == PAIRED;
                        }
                    }
                    auto score_a1 = a.first.has_value() ? a.first.value().alignment_ani : 0;
                    auto score_a2 = a.second.has_value() ? a.second.value().alignment_ani : 0;
                    auto score_a = std::max(score_a1, score_a2);

                    auto score_b1 = b.first.has_value() ? b.first.value().alignment_ani : 0;
                    auto score_b2 = b.second.has_value() ? b.second.value().alignment_ani : 0;
                    auto score_b = std::max(score_b1, score_b2);


                    return score_a > score_b;
                });

                auto &best = non_uniques.front();

                auto [ani, score] = ScorePairedAlignment(best);

                auto qual = GetReadEvidenceQuality(ani, score);

                for (auto& nu : non_uniques) {
                    if (nu.first.has_value() && nu.second.has_value() && nu.first->readid != nu.second->readid) {
                        std::cout << "Error" << std::endl;
                        exit(9);
                    }
                }

                if (!IsAmbiguous(non_uniques)) {
                    m_alignments.GetNonUniquePairList(qual).emplace_back(std::move(non_uniques));
                }

                non_uniques.clear();
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

//                std::cout << pairless_header1 << "\n" << pairless_header2 << "\n--" << std::endl;

                return pairless_header1 == pairless_header2;
            }

            void ProcessNonUnique(std::optional<InternalReadAlignment> const &ira, std::optional<InternalReadAlignment> const &ira2) {
                static int good_qual_counter = 0;
                static int non_unique_same_id_counter = 0;
                if (!m_tmp_pair.empty() && !SameReadId(ira, ira2, m_tmp_pair.front())) {
                    ProcessNonUniqueWorker(m_tmp_pair);
                    non_unique_same_id_counter++;
                }
//                std::cout << "Emplace back " << static_cast<int64_t>(ira.has_value() ? ira->readid : -1) << " " << static_cast<int64_t>(ira2.has_value() ? ira2->readid : -1) << std::endl;
                if (ira.has_value() && ira2.has_value() && ira->readid != ira2->readid) {
                    std::cout << "Error" << std::endl;
                    exit(9);
                }
                m_tmp_pair.emplace_back(OptIRAPair { ira, ira2 });
            }

            void Process(std::optional<InternalReadAlignment> const &ira,
                         std::optional<InternalReadAlignment> const &ira2,
                         bool unique) {
                assert(ira.has_value() || ira2.has_value());
                if (unique) {
//                    std::cout << "Unique: " << std::endl;
                    ProcessUnique(ira, ira2);
//                    std::cout << (ira.has_value() ? ira.value().ToString() : "") << std::endl;
//                    std::cout <<  (ira2.has_value() ? ira2.value().ToString() : "") << std::endl;
                } else {
//                    std::cout << "Non-Unique: " << std::endl;
//                    std::cout << (ira.has_value() ? ira.value().ToString() : "") << std::endl;
//                    std::cout <<  (ira2.has_value() ? ira2.value().ToString() : "") << std::endl;
                    ProcessNonUnique(ira, ira2);
                }
//                Utils::Input();
            }


//            void RegisterRead(size_t read_id, TaxId taxid, GeneId geneid, GenePos genepos,
//                              bool forward, size_t pairid, int alignment_score,
//                              double alignment_ani, bool unique = true) {
//
//                Process(std::optional<InternalReadAlignment>{
//                    InternalReadAlignment(
//                        read_id, pairid, taxid, geneid, genepos, forward,
//                        alignment_score, alignment_ani)
//                    }, unique);
//            }

            void Process(size_t read_id, SamEntry &sam1, SamEntry &sam2, bool has_sam1, bool has_sam2, bool unique) {
                CigarInfo info;
                size_t min_alignment_length = 70;
                if (has_sam1) {
                    CompressedCigarInfo(sam1.m_cigar, info);
                    if (info.clipped_alignment_length < min_alignment_length) has_sam1 = false;
                }
                if (has_sam2) {
                    CompressedCigarInfo(sam2.m_cigar, info);
                    if (info.clipped_alignment_length < min_alignment_length) has_sam2 = false;
                }
                if (!has_sam1 && !has_sam2) return;
                OptIRA oira1 = has_sam1 ? OptIRA(InternalReadAlignment(read_id, sam1)) : OptIRA();
                OptIRA oira2 = has_sam2 ? OptIRA(InternalReadAlignment(read_id, sam2)) : OptIRA();
                Process(oira1, oira2, unique);
            }


            void FromSam(std::string &file_path) {
                std::ifstream file(file_path, std::ios::in);

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
                    std::cerr << "Profile is empty" << std::endl;
                    exit(9);
                }

                while (line[0] == '@') {
                    std::getline(file, line);
                }

                GetSamPair(file, line, tokens, current, current_other, has_current1, has_current2, true);
                //last_qname = has_current1 ? current.m_qname.substr(0, current.m_qname.length()-2) : current_other.m_qname.substr(0, current_other.m_qname.length()-2);

//                std::cout << "Has current1: " << has_current1 << " " << current.m_rname << std::endl;
//                std::cout << "Has current2: " << has_current2 << " " << current_other.m_rname << std::endl;
                size_t read_id = 0;
                while (GetSamPair(file, line, tokens, next, next_other, has_next1, has_next2)) {

                    current_qname = has_current1 ? current.m_qname.substr(0, current.m_qname.length()-2) : current_other.m_qname.substr(0, current_other.m_qname.length()-2);
                    next_qname = has_next1 ? next.m_qname.substr(0, next.m_qname.length()-2) : next_other.m_qname.substr(0, next_other.m_qname.length()-2);
                    bool unique = current_qname != last_qname && current_qname != next_qname;

//                    if (!unique) {
//                        std::cout << "Has current1: " << has_current1 << " " << current.m_rname << std::endl;
//                        std::cout << "Has current2: " << has_current2 << " " << current_other.m_rname << std::endl;
//                        std::cout << "Has next1: " << has_next1 << " " << next.m_rname << std::endl;
//                        std::cout << "Has next2: " << has_next2 << " " << next_other.m_rname << std::endl;
//                    }
//                    std::cout << "current_qname: " << current_qname << " <<<<" << std::endl;

                    read_id += current_qname != last_qname;
                    last_qname = current_qname;

//                    std::cout << "------------------Pair" << std::endl;
//                    std::cout << std::string(79, '-') << " unique? " << unique << std::endl;
//                    if (has_current1) {
//                        std::cout << read_id << " " << current.ToString() << std::endl;
//                    }
//                    if (has_current2) {
//                        std::cout << read_id << " " << current_other.ToString() << std::endl;
//                    }
//                    if (has_next1) {
//                        std::cout << next.ToString() << std::endl;
//                    }
//                    if (has_next2) {
//                        std::cout << next_other.ToString() << std::endl;
//                    }
//                    Utils::Input();

                    Process(read_id, current, current_other, has_current1, has_current2, unique);
//                    std::cout << "After Process" << std::endl;

                    std::swap(current, next);
                    std::swap(current_other, next_other);
                    has_current1 = has_next1;
                    has_current2 = has_next2;
                }
                ProcessNonUniqueWorker(m_tmp);

                current_qname = has_current1 ? current.m_qname : current_other.m_qname;
                bool unique = current_qname != last_qname;
                read_id += unique;
                Process(read_id, current, current_other, has_current1, has_current2, unique);


//                return Profile();
            }


            void FromSam(std::string &file_path, std::vector<std::vector<AlignmentPair>> &pairs, std::vector<AlignmentPair> &unique_pairs) {
                std::ifstream file(file_path, std::ios::in);

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
                    std::cerr << "Profile is empty" << std::endl;
                    exit(9);
                }

                while (line[0] == '@') {
                    std::getline(file, line);
                }

                GetSamPair(file, line, tokens, current, current_other, has_current1, has_current2, true);

                std::vector<AlignmentPair> pair_list;

                size_t read_id = 0;
                while (GetSamPair(file, line, tokens, next, next_other, has_next1, has_next2)) {

                    current_qname = has_current1 ? current.m_qname.substr(0, current.m_qname.length()-2) : current_other.m_qname.substr(0, current_other.m_qname.length()-2);
                    next_qname = has_next1 ? next.m_qname.substr(0, next.m_qname.length()-2) : next_other.m_qname.substr(0, next_other.m_qname.length()-2);
                    bool unique = current_qname != last_qname && current_qname != next_qname;
                    last_qname = current_qname;

                    AlignmentPair new_pair = AlignmentPair(
                        has_current1 ? std::optional<SamEntry>{current} : std::optional<SamEntry>{},
                        has_current2 ? std::optional<SamEntry>{current_other} : std::optional<SamEntry>{}
                    );

                    if (unique) {
                        unique_pairs.emplace_back(std::move(new_pair));
                    } else {
                        if (!pair_list.empty() && !SameRead(new_pair, pair_list.front())) {
                            pairs.emplace_back(std::move(pair_list));
                            pair_list.clear();
                        }
                        pair_list.emplace_back(new_pair);
                    }

                    std::swap(current, next);
                    std::swap(current_other, next_other);
                    has_current1 = has_next1;
                    has_current2 = has_next2;
                }


                AlignmentPair new_pair = AlignmentPair(
                        has_current1 ? std::optional<SamEntry>{current} : std::optional<SamEntry>{},
                        has_current2 ? std::optional<SamEntry>{current_other} : std::optional<SamEntry>{}
                );
                current_qname = has_current1 ? current.m_qname : current_other.m_qname;
                bool unique = current_qname != last_qname;
                if (unique) {
                    unique_pairs.emplace_back(std::move(new_pair));
                } else {
                    if (!pair_list.empty() && !SameRead(new_pair, pair_list.front())) {
                        pairs.emplace_back(std::move(pair_list));
                        pair_list.clear();
                    }
                    pair_list.emplace_back(new_pair);
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
                        gene.Load();
                        auto& ref = gene.Sequence();
//                        std::cout << "Extract: " << any.m_qname << " ---> " <<  any.m_rname << std::endl;

//                        PrintAlignment(any, ref);
                        ExtractSNPs(any, ref, snps, tid, gid);
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
                for (auto& pairlist : pair_lists) {
                    for (auto& pair : pairlist) {
                        auto &[ tid, gid ] = ExtractTaxidGeneid(pair.Any().m_rname);
                        auto &[ atid, agid ] = ExtractTaxidGeneid(pair.Any().m_qname);
//                        std::cerr << pair.Any().m_qname << " " << atid << " " << agid << std::endl;
                        std::cerr << "CAPTURE\t";
                        if (pair.IsPair()) {
                            auto &[ tid1, gid1 ] = ExtractTaxidGeneid(pair.First().m_rname);
                            auto &[ tid2, gid2 ] = ExtractTaxidGeneid(pair.Second().m_rname);
                            auto &[ ttid1, tgid1 ] = ExtractTaxidGeneid(pair.First().m_qname);
                            auto &[ ttid2, tgid2 ] = ExtractTaxidGeneid(pair.Second().m_qname);
                            CompressedCigarInfo(pair.First().m_cigar, info1);
                            CompressedCigarInfo(pair.Second().m_cigar, info2);
                            bool same_gene = gid1 == gid2;
                            int insert = std::max(pair.First().m_pos, pair.Second().m_pos) - std::min(pair.First().m_pos, pair.Second().m_pos);

                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag),  tid == ttid1, same_gene, insert, true, m_genome_loader.GetGenome(tid1).GetGene(gid1).GetLength()) << std::endl;
                        } else if (pair.HasFirst()) {
                            CompressedCigarInfo(pair.First().m_cigar, info1);
                            auto &[ ttid, tgid ] = ExtractTaxidGeneid(pair.First().m_qname);
                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), -1, -1, false, tid == ttid, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
                        } else if (pair.HasSecond()) {
                            CompressedCigarInfo(pair.Second().m_cigar, info2);
                            auto &[ ttid, tgid ] = ExtractTaxidGeneid(pair.Second().m_qname);
                            std::cerr << ErrorLine(pairlist_id, -1, -1, false, info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag), tid == ttid, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
                        }
                    }

                    pairlist_id++;
                }
            }

            void OutputErrorData(TruthSet const& set, std::vector<std::vector<AlignmentPair>> &pair_lists) {
                int pairlist_id = 0;
                CigarInfo info1;
                CigarInfo info2;
                for (auto& pairlist : pair_lists) {
                    for (auto& pair : pairlist) {
                        auto &[ tid, gid ] = ExtractTaxidGeneid(pair.Any().m_rname);
                        bool maybe_true = set.contains(tid);
                        std::cerr << "CAPTURE\t";
                        if (pair.IsPair()) {
                            auto &[ tid1, gid1 ] = ExtractTaxidGeneid(pair.First().m_rname);
                            auto &[ tid2, gid2 ] = ExtractTaxidGeneid(pair.Second().m_rname);
                            CompressedCigarInfo(pair.First().m_cigar, info1);
                            CompressedCigarInfo(pair.Second().m_cigar, info2);
                            bool same_gene = gid1 == gid2;
                            int insert = std::max(pair.First().m_pos, pair.Second().m_pos) - std::min(pair.First().m_pos, pair.Second().m_pos);

                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag),  maybe_true, same_gene, insert, true, m_genome_loader.GetGenome(tid1).GetGene(gid1).GetLength()) << std::endl;
                        } else if (pair.HasFirst()) {
                            CompressedCigarInfo(pair.First().m_cigar, info1);
                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), -1, -1, false, maybe_true, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
                        } else if (pair.HasSecond()) {
                            CompressedCigarInfo(pair.Second().m_cigar, info2);
                            std::cerr << ErrorLine(pairlist_id, -1, -1, false, info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag), maybe_true, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
                        }
                    }

                    pairlist_id++;
                }
            }

            MicrobialProfile ProfileWithTruth(TruthSet const& set) {
                MicrobialProfile profile(m_genome_loader);
                profile.SetName("Profile");

                size_t total_reads_profiled = 0;

                size_t uniques_total = 0;
                size_t uniques_fp = 0;
                size_t nonuniques_total = 0;
                size_t nonuniques_fp = 0;

                std::cout << "Uniques" << std::endl;
                for (int i = 0; i < ReadEvidenceQuality::SIZE; i++) {
                    auto qual = quality_enum_list[i];
                    auto& uniques = m_alignments.GetUniqueList(qual);
                    std::cout << qual << " " << uniques.size() << std::endl;
                    for (auto& ira : uniques) {
                        total_reads_profiled++;
                        if (!set.contains(ira.taxid)) {
                            uniques_fp++;
                        }
                        uniques_total++;
                        profile.AddRead(ira, true);
                    }

                }

                std::cout << "Non Unique Pairs" << std::endl;
                for (int i = 0; i < ReadEvidenceQuality::SIZE; i++) {
                    auto qual = quality_enum_list[i];
                    auto& non_uniques = m_alignments.GetNonUniquePairList(qual);
                    std::cout << qual << " " << non_uniques.size() << std::endl;
                    for (auto& ira_pair_list : non_uniques) {
                        bool fp = false;
                        if (ira_pair_list.front().first.has_value()) {
                            profile.AddRead(ira_pair_list.front().first.value(), false);
                            total_reads_profiled++;

                            // DEBUG
                            if (!set.contains(ira_pair_list.front().first.value().taxid)) {
                                fp = true;
                                nonuniques_fp++;
                            }
                            nonuniques_total++;
                            //DEBUG
                        }
                        if (ira_pair_list.front().second.has_value()) {
                            profile.AddRead(ira_pair_list.front().second.value(), false);
                            total_reads_profiled++;

                            // DEBUG
                            if (!set.contains(ira_pair_list.front().second.value().taxid)) {
                                fp = true;
                                nonuniques_fp++;
                            }
                            nonuniques_total++;
                            //DEBUG
                        }
//                        if (fp) {
//                            std::cout << "Non uniques: " << ira_pair_list.size() << std::endl;
//                        }
                        for (auto& ira_pair : ira_pair_list) {
//                            if (fp) {
//                                bool first_correct = ira_pair.first.has_value() && set.contains(ira_pair.first.value().taxid);
//                                bool second_correct = ira_pair.second.has_value() && set.contains(ira_pair.second.value().taxid);
//                                std::string tr = first_correct || second_correct ? "true" : "false";
//                                std::cout << (ira_pair.first.has_value() ? ira_pair.first->ToString() + " " + tr + "\n" : "");
//                                std::cout << (ira_pair.second.has_value() ? ira_pair.second->ToString() + " " + tr + "\n" : "");
//                                std::cout << "---" << std::endl;
//                            }
                            if (ira_pair.first.has_value() && ira_pair.second.has_value() && ira_pair.first->readid != ira_pair.second->readid) {
                                std::cout << "Erroneous" << std::endl;
                                exit(9);
                            }
                        }
//                        if (fp) Utils::Input();
                    }
                }

                std::cout << "Non Uniques" << std::endl;
                for (int i = 0; i < ReadEvidenceQuality::SIZE; i++) {
                    auto qual = quality_enum_list[i];
                    auto& non_uniques = m_alignments.GetNonUniqueList(qual);
                    std::cout << qual << " " << non_uniques.size() << std::endl;
                    for (auto& ira : non_uniques) {
                        profile.AddRead(ira.front(), false);
                        total_reads_profiled++;
                    }
                }

                std::cout << "Unique fp rate:     " << uniques_fp << "/" << uniques_total << "\t\t" << static_cast<double>(uniques_fp) / static_cast<double>(uniques_total) << std::endl;
                std::cout << "Non-Unique fp rate: " << nonuniques_fp << "/" << nonuniques_total << "\t\t" << static_cast<double>(nonuniques_fp) / static_cast<double>(nonuniques_total) << std::endl;

                std::cout << "Reads contributed to profiles: " << total_reads_profiled << std::endl;
                return profile;
            }

            MicrobialProfile Profile() {
                MicrobialProfile profile(m_genome_loader);
                profile.SetName("Profile");

                size_t total_reads_profiled = 0;

                std::cout << "Uniques" << std::endl;
                for (int i = 0; i < ReadEvidenceQuality::SIZE; i++) {
                    auto qual = quality_enum_list[i];
                    auto& uniques = m_alignments.GetUniqueList(qual);
                    std::cout << qual << " " << uniques.size() << std::endl;
                    for (auto& ira : uniques) {
                        total_reads_profiled++;
                        profile.AddRead(ira, true);
                    }
                }

                std::cout << "Non Unique Pairs" << std::endl;
                for (int i = 0; i < ReadEvidenceQuality::SIZE; i++) {
                    auto qual = quality_enum_list[i];
                    auto& non_uniques = m_alignments.GetNonUniquePairList(qual);
                    std::cout << qual << " " << non_uniques.size() << std::endl;
                    for (auto& ira_pair_list : non_uniques) {
                        if (ira_pair_list.front().first.has_value()) {
                            profile.AddRead(ira_pair_list.front().first.value(), false);
                            total_reads_profiled++;
                        }
                        if (ira_pair_list.front().second.has_value()) {
                            profile.AddRead(ira_pair_list.front().second.value(), false);
                            total_reads_profiled++;
                        }
                        std::cout << "Non uniques: " << ira_pair_list.size() << std::endl;
                        for (auto& ira_pair : ira_pair_list) {
                            std::cout << (ira_pair.first.has_value() ? ira_pair.first->ToString() : "") << std::endl;
                            std::cout << (ira_pair.second.has_value() ? ira_pair.second->ToString() : "") << std::endl;
                            std::cout << "---" << std::endl;
                            if (ira_pair.first.has_value() && ira_pair.second.has_value() && ira_pair.first->readid != ira_pair.second->readid) {
                                std::cout << "Erroneous" << std::endl;
                                exit(9);
                            }
                        }
                    }
                }

                std::cout << "Non Uniques" << std::endl;
                for (int i = 0; i < ReadEvidenceQuality::SIZE; i++) {
                    auto qual = quality_enum_list[i];
                    auto& non_uniques = m_alignments.GetNonUniqueList(qual);
                    std::cout << qual << " " << non_uniques.size() << std::endl;
                    for (auto& ira : non_uniques) {
                        profile.AddRead(ira.front(), false);
                        total_reads_profiled++;
                    }
                }

                std::cout << "Reads contributed to profiles: " << total_reads_profiled << std::endl;
                return profile;
            }


        };
        const Profiler::ReadEvidenceQuality Profiler::quality_enum_list[] = {
            Profiler::ReadEvidenceQuality::EXCELLENT,
            Profiler::ReadEvidenceQuality::GOOD,
            Profiler::ReadEvidenceQuality::MEDIUM,
            Profiler::ReadEvidenceQuality::BAD
        };

        const std::string Profiler::quality_map[] = {
            "Excellent",
            "Good",
            "Medium",
            "Bad"
        };
    }
}