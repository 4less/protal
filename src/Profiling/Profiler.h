//
// Created by fritsche on 07/10/22.
//

#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <sparse_map.h>
#include "LineSplitter.h"

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
            std::cout << "truth " << truth << std::endl;
            truths.insert(truth);
        }
        return truths;
    }

    namespace profiler {
        struct MicrobialProfile;
        struct InternalReadAlignment;
        struct Gene;
        struct Taxon;

        using TaxId = uint32_t;
        using GeneId = uint32_t;
        using GenePos = uint32_t;

        using InternalReadAlignmentList = std::vector<InternalReadAlignment>;
        using InternalNonUniqueReadAlignmentList = std::vector<InternalReadAlignment>;
        using InternalNonUniqueReadAlignmentListList = std::vector<InternalNonUniqueReadAlignmentList>;

        using GeneMap = tsl::sparse_map<uint32_t, Gene>;
        using TaxonMap = tsl::sparse_map<uint32_t, Taxon>;

        // Data structures

        class Gene {
        public:
            size_t m_mapped_reads = 0;

            void AddRead(GenePos pos) {
                m_mapped_reads++;
            }
        };

        class Taxon {
        private:
            std::string m_name;
            size_t m_total_hits = 0;
            double m_ani_sum = 0;
            GeneMap m_genes;

        public:
            void AddHit(GeneId geneid, GenePos genepos, double ani) {
                if (m_genes.contains(geneid)) {
                    m_genes[geneid] = Gene{};
                }
                m_genes[geneid].AddRead(genepos);
                m_total_hits++;
                m_ani_sum += ani;
            }

            size_t GeneNum() const {
                return m_genes.size();
            }

            double GetMeanANI() const {
                return m_ani_sum/m_total_hits;
            }

            size_t TotalHits() const {
                return m_total_hits;
            }

            void SetName(std::string name) {
                m_name = name;
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

        struct InternalReadAlignment {
            size_t readid;
            size_t pairid;
            TaxId taxid;
            GeneId geneid;
            GenePos genepos;
            int alignment_score;
            double alignment_ani;
            bool forward;

            std::string ToString() {
                std::string str;
                str += std::to_string(readid) + '\t';
                str += std::to_string(pairid) + '\t';
                str += std::to_string(taxid) + '\t';
                str += std::to_string(geneid) + '\t';
                str += std::to_string(genepos) + '\t';
                str += std::to_string(alignment_score) + '\t';
                str += std::to_string(alignment_ani) + '\t';
                str += std::to_string(forward);
                return str;
            }

            InternalReadAlignment(
                    size_t readid,
                    size_t pairid,
                    TaxId taxid,
                    GeneId geneid,
                    GenePos genepos,
                    bool forward,
                    int alignment_score,
                    double alignment_ani) :
                    readid(readid),
                    taxid(taxid),
                    geneid(geneid),
                    genepos(genepos),
                    forward(forward),
                    pairid(pairid),
                    alignment_score(alignment_score),
                    alignment_ani(alignment_ani) {};


            InternalReadAlignment(InternalReadAlignment const &other) :
                    readid(other.readid),
                    taxid(other.taxid),
                    geneid(other.geneid),
                    genepos(other.genepos),
                    forward(other.forward),
                    pairid(other.pairid),
                    alignment_score(other.alignment_score),
                    alignment_ani(other.alignment_ani) {};
        };

        class MicrobialProfile {
        public:
            void AddRead(InternalReadAlignment const& ira) {
                if (!m_taxa.contains(ira.taxid)) {
                    m_taxa.insert( { ira.taxid, Taxon() } );
                    m_taxa[ira.taxid].SetName(std::to_string(ira.taxid));
                }
                auto& taxon = m_taxa[ira.taxid];
                taxon.AddHit(ira.geneid, ira.genepos, ira.alignment_ani);
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
                std::cout << "Truth: ________________" << std::endl;
                for (auto& [key, taxon] : m_taxa) {
                    bool tp = set.contains(key);
                    std::cout << tp << "\t";
                    std::cout << key << "\t";
                    std::cout << taxon.GeneNum() << "\t";
                    std::cout << taxon.TotalHits() << "\t";
                    std::cout << taxon.GetMeanANI() << std::endl;
                }
            }

        private:
            std::string m_name;
            TaxonMap m_taxa;
        };

        class Profiler {
        public:

            enum ReadEvidenceQuality {
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
                size_t m_size = ReadEvidenceQuality::SIZE;


            public:
                AlignmentLists() {
                    for (int i = 0; i < m_size; i++)
                        m_unique_alignments.emplace_back(InternalReadAlignmentList{});
                    for (int i = 0; i < m_size; i++)
                        m_non_unique_alignments.emplace_back(InternalNonUniqueReadAlignmentListList{});
                };

                InternalReadAlignmentList &GetUniqueList(ReadEvidenceQuality quality) {
                    return m_unique_alignments[quality];
                }

                InternalNonUniqueReadAlignmentListList &GetNonUniqueList(ReadEvidenceQuality quality) {
                    return m_non_unique_alignments[quality];
                }
            };

        private:
            AlignmentLists m_alignments;
            InternalNonUniqueReadAlignmentList m_tmp;

        public:


            void
            ProcessUnique(InternalReadAlignment const &ira, ReadEvidenceQuality min_qual = ReadEvidenceQuality::GOOD) {
                auto qual = GetReadEvidenceQuality(ira.alignment_ani, ira.alignment_score);
                static int good_qual_counter = 0;
//                if (qual >= min_qual) {
//                    good_qual_counter++;
//                    std::cout << "good_unique_qual_counter: " << good_qual_counter << std::endl;
//                }
                m_alignments.GetUniqueList(qual).emplace_back(ira);
            }

            void ProcessNonUniqueWorker(std::vector<InternalReadAlignment> &non_uniques) {
                if (non_uniques.empty()) return;
                std::sort(non_uniques.begin(), non_uniques.end(),
                          [](InternalReadAlignment const &a, InternalReadAlignment const &b) {
                              return a.alignment_score > b.alignment_score;
                          });
                auto &best = non_uniques.front();
                auto qual = GetReadEvidenceQuality(best.alignment_ani, best.alignment_score);

//                std::cout << quality_map[qual];
//                std::cout << " These for a non-unique cluster with same read_name" << std::endl;
//                for (auto a : non_uniques) {
//                    std::cout << a.ToString() << std::endl;
//                }
//                std::cout << "-----" << std::endl;

                m_alignments.GetNonUniqueList(qual).emplace_back(std::move(non_uniques));
                non_uniques.clear();
            }

            void ProcessNonUnique(InternalReadAlignment const &ira) {
                static int good_qual_counter = 0;
                if (!m_tmp.empty() && m_tmp.front().readid != ira.readid) {
                    ProcessNonUniqueWorker(m_tmp);
                }
                m_tmp.emplace_back(ira);
            }

            void Process(InternalReadAlignment const &ira, bool unique) {
                if (unique) {
                    ProcessUnique(ira);
                } else {
                    ProcessNonUnique(ira);
                }
            }

            void RegisterRead(size_t read_id, TaxId taxid, GeneId geneid, GenePos genepos,
                              bool forward, size_t pairid, int alignment_score,
                              double alignment_ani, bool unique = true) {
                Process(InternalReadAlignment{
                        read_id, pairid, taxid, geneid, genepos, forward,
                        alignment_score, alignment_ani
                }, unique);
            }

            void SamFromTokens(std::vector<std::string>& tokens, SamEntry &sam) {
                sam.m_qname = tokens[0];
                sam.m_flag = std::stoul(tokens[1]);
                sam.m_rname = tokens[2];
                sam.m_pos = std::stoul(tokens[3]);
                sam.m_mapq = std::stoul(tokens[4]);
                sam.m_cigar = tokens[5];
                sam.m_rnext = tokens[6];
                sam.m_pnext = std::stoul(tokens[7]);
                sam.m_tlen = std::stol(tokens[8]);
                sam.m_seq = tokens[9];
                sam.m_qual = tokens[10];
            }

            MicrobialProfile FromSam(std::string &file_path) {
                std::ifstream file(file_path, std::ios::in);

                std::string delim = "\t";
                std::string line;
                std::vector<std::string> tokens;
                SamEntry current;
                SamEntry next;

                std::string last = "";


                if (!std::getline(file, line)) {
                    std::cerr << "Profile is empty" << std::endl;
                    return MicrobialProfile();
                }
                LineSplitter::Split(line, delim, tokens);
                SamFromTokens(tokens, current);

                size_t read_id = 0;
                while (std::getline(file, line)) {
                    LineSplitter::Split(line, delim, tokens);
                    SamFromTokens(tokens, next);
                    bool unique = current.m_qname != last && current.m_qname != next.m_qname;

                    LineSplitter::Split(current.m_rname, "_", tokens);

                    size_t taxid = std::stoul(tokens[0]);
                    size_t geneid = std::stoul(tokens[1]);
                    size_t pos = current.m_pos;
                    double ani = WFA2Wrapper::CigarANI(current.m_cigar);
                    int score = WFA2Wrapper::CigarScore(current.m_cigar);

                    read_id += current.m_qname != last;

                    RegisterRead(read_id, taxid, geneid, pos, true, 0, score, ani, unique);

                    last = current.m_qname;
                    std::swap(current, next);
                }
                ProcessNonUniqueWorker(m_tmp);

                bool unique = current.m_qname != last;
                LineSplitter::Split(current.m_rname, "_", tokens);
                size_t taxid = std::stoul(tokens[0]);
                size_t geneid = std::stoul(tokens[1]);
                size_t pos = current.m_pos;
                double ani = WFA2Wrapper::CigarANI(current.m_cigar);
                int score = WFA2Wrapper::CigarScore(current.m_cigar);
                RegisterRead(read_id, taxid, geneid, pos, true, 0, score, ani, unique);


                return Profile();
            }

            MicrobialProfile Profile() {
                MicrobialProfile profile;
                profile.SetName("Profile");

                std::cout << "Unique profile: " << std::endl;
                for (int i = 0; i < ReadEvidenceQuality::SIZE; i++) {
                    auto qual = quality_enum_list[i];
                    std::cout << "\n\nqual: " << quality_map[qual] << std::string(40, '_') << std::endl;
                    auto& uniques = m_alignments.GetUniqueList(qual);
                    for (auto& ira : uniques) {
                        profile.AddRead(ira);
                    }
                    std::cout << profile.ToString() << std::endl;
                    for (auto& [key, taxon] : profile.GetTaxa()) {
                        std::cout << taxon.ToString() << std::endl;
                    }
                    std::cout << "\n\nUNIQUE qual: " << quality_map[qual] << std::string(10, '^') << std::endl;
                    Utils::Input();
                }

                for (int i = 0; i < ReadEvidenceQuality::SIZE; i++) {
                    auto qual = quality_enum_list[i];
                    std::cout << "\n\nNONUNIQUE qual: " << quality_map[qual] << std::string(40, '_') << std::endl;
                    auto& uniques = m_alignments.GetNonUniqueList(qual);
                    for (auto& ira : uniques) {
                        profile.AddRead(ira.front());
                    }
                    std::cout << profile.ToString() << std::endl;
                    for (auto& [key, taxon] : profile.GetTaxa()) {
                        std::cout << taxon.ToString() << std::endl;
                    }
                    std::cout << "\n\nNONUNIQUE qual: " << quality_map[qual] << std::string(10, '^') << std::endl;

                    Utils::Input();
                }

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