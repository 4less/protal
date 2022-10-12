//
// Created by fritsche on 07/10/22.
//

#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <sparse_map.h>

namespace protal {
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
        GeneMap m_genes;

    public:
        void AddHit(GeneId geneid, GenePos genepos) {
            if (m_genes.contains(geneid)) {
                m_genes[geneid] = Gene{};
            }
            m_genes[geneid].AddRead(genepos);
        }
        size_t GeneNum() {
            return m_genes.size();
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
        void AddRead(InternalReadAlignment) {

        }
    private:
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


        void ProcessUnique(InternalReadAlignment const &ira, ReadEvidenceQuality min_qual = ReadEvidenceQuality::GOOD) {
            auto qual = GetReadEvidenceQuality(ira.alignment_ani, ira.alignment_score);

            if (qual >= min_qual) {

            }

            m_alignments.GetUniqueList(qual).emplace_back(ira);
        }

        void ProcessNonUnique(InternalReadAlignment const &ira) {
            if (!m_tmp.empty() && m_tmp.front().readid != ira.readid) {
                std::sort(m_tmp.begin(), m_tmp.end(),
                          [](InternalReadAlignment const &a, InternalReadAlignment const &b) {
                              return a.alignment_score < b.alignment_score;
                          });
                auto &best = m_tmp.front();
                auto qual = GetReadEvidenceQuality(best.alignment_ani, best.alignment_score);
                m_alignments.GetNonUniqueList(qual).emplace_back(m_tmp);
                m_tmp.clear();
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

    };
}