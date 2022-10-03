//
// Created by fritsche on 27/09/22.
//

#pragma once

#include <vector>
#include "BinaryClassifierEvaluator.h"

#include "SeedingStrategy.h"
#include "AlignmentOutputHandler.h"
#include "Constants.h"

namespace protal {
    using BCE = BinaryClassifierEvaluator;
    class CoreBenchmark {
        BCE m_seed_bce{"Seeds"};
        BCE m_anchor_bce{"Anchors"};
        BCE m_alignment_bce{"Alignments"};
        BCE m_best_alignment_bce{"Best_alignment"};
        BCE m_best_unique_alignment_bce{"Best_unique_alignment"};

        size_t m_anchor_fp_no_hit = 0;

    public:
        static std::pair<uint32_t, uint32_t> HeaderToTruth(std::string const& header) {
            if (header.empty()) {
                std::cerr << "Tried to extract taxonomic id and gene id from empty header" << std::endl;
            }
            uint32_t taxid = 0;
            uint32_t geneid = 0;

            size_t start = 0;
            size_t end = 0;
            while (end < header.length() && header[end] != '_') end++;

            if (end == header.length()) {
                std::cerr << "Tried to extract taxonomic id and gene id from header "  << header << std::endl;
                exit(9);
            }
            taxid = std::stoul(header.substr(start, end));
            end++;
            start = end;
            while (end < header.length() && isdigit(header[end])) {
                end++;
            }

            if (end - start == 0) {
                std::cerr << "Tried to extract taxonomic id and gene id from header "  << header << std::endl;
                exit(9);
            }
            geneid = std::stoul(header.substr(start, end));

            return { taxid, geneid };
        }

        void AddSeeds(SeedList const& seeds, size_t true_taxid, size_t true_geneid) {
            for (auto& seed : seeds) {
                if (seed.taxid == true_taxid && seed.geneid == true_geneid) {
                    m_seed_bce.tp++;
                    return;
                }
            }

            m_seed_bce.fn++;
        };

        void AddAnchors(AlignmentAnchorList const& anchors, size_t true_taxid, size_t true_geneid) {
            for (auto& anchor : anchors) {
                if (anchor.a.taxid == true_taxid && anchor.a.geneid == true_geneid) {
                    m_anchor_bce.tp++;
                    return;
                }
            }

            m_anchor_fp_no_hit += anchors.empty();

            m_anchor_bce.fn++;
        };

        void AddAlignmentResults(AlignmentResultList const& alignments, size_t true_taxid, size_t true_geneid) {
            if (alignments.empty()) {
                m_alignment_bce.fn++;
                m_best_alignment_bce.fn++;
                m_best_unique_alignment_bce.fn++;
                return;
            }

            auto& best_alignment = alignments[0];
            if (best_alignment.Taxid() == true_taxid && best_alignment.GeneId() == true_geneid) {
                m_best_alignment_bce.tp++;
            } else {
                m_best_alignment_bce.fp++;
            }

            if (alignments.size() == 1 || (alignments.size() > 1 && alignments[1].AlignmentScore() < best_alignment.AlignmentScore())) {

                if (best_alignment.Taxid() == true_taxid && best_alignment.GeneId() == true_geneid) {
                    m_best_unique_alignment_bce.tp++;
                } else {
                    m_best_unique_alignment_bce.fp++;
                }
            }

            for (auto& alignment : alignments) {
                if (alignment.Taxid() == true_taxid && alignment.GeneId() == true_geneid) {
                    m_alignment_bce.tp++;
                    break;
                }
            }
        };

        void AddBestAlignmentResult(AlignmentResult& alignment, size_t true_taxid, size_t true_geneid) {
            if (alignment.Taxid() == true_taxid && alignment.GeneId() == true_geneid) {
                m_best_alignment_bce.tp++;
            } else {
                m_best_alignment_bce.fp++;
            }
        }

        void AddNoBestAlignmentResult() {
            m_best_alignment_bce.fn++;
        }

        void Join(CoreBenchmark const& other) {
            m_seed_bce.Join(other.m_seed_bce);
            m_anchor_bce.Join(other.m_anchor_bce);
            m_alignment_bce.Join(other.m_alignment_bce);
            m_best_alignment_bce.Join(other.m_best_alignment_bce);
            m_best_unique_alignment_bce.Join(other.m_best_unique_alignment_bce);
            m_anchor_fp_no_hit += other.m_anchor_fp_no_hit;
        }

        void operator() (SeedList const& seeds, AlignmentAnchorList const& anchors, AlignmentResultList const& alignments, std::string const& header) {
            auto [tax_id_truth, gene_id_truth] = CoreBenchmark::HeaderToTruth(header);

            AddSeeds(seeds, tax_id_truth, gene_id_truth);
            AddAnchors(anchors, tax_id_truth, gene_id_truth);
            AddAlignmentResults(alignments, tax_id_truth, gene_id_truth);

        }

        void WriteRowStats() {
            m_seed_bce.WriteRowHeader();
            m_seed_bce.WriteRowStats();
            m_anchor_bce.WriteRowStats();
            m_alignment_bce.WriteRowStats();
            m_best_alignment_bce.WriteRowStats();
            m_best_unique_alignment_bce.WriteRowStats();
        }

        void Print() {
            std::cout << "Anchor stats: " << std::endl;
            m_anchor_bce.WriteStats();
            std::cout << "Alignment stats: " << std::endl;
            m_alignment_bce.WriteStats();
            std::cout << "False positive anchors without hit " << m_anchor_fp_no_hit << std::endl;
        }
    };
}

