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
        BCE m_anchor_bce;
        BCE m_alignment_bce;
        BCE m_best_alignment_bce;

        size_t m_anchor_fp_no_hit = 0;

    public:
        static std::pair<uint32_t, uint32_t> ExtractTruthFromHeader(std::string& header) {
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

        void AddAnchors(AlignmentAnchorList& anchors, size_t true_taxid, size_t true_geneid) {
            for (auto& anchor : anchors) {
                if (anchor.a.taxid == true_taxid && anchor.a.geneid == true_geneid) {
                    m_anchor_bce.tp++;
                    return;
                }
            }

            m_anchor_fp_no_hit += anchors.empty();

            m_anchor_bce.fn++;
        };

        void AddAlignmentResults(AlignmentAnchorList &anchors, AlignmentResultList& alignments, size_t true_taxid, size_t true_geneid) {
            for (auto& alignment : alignments) {
                if (alignment.Taxid() == true_taxid && alignment.GeneId() == true_geneid) {
                    m_alignment_bce.tp++;
                    return;
                }
            }

            if (!anchors.empty() && false) {
                std::cout << std::string(50, '-') << " False negative:" << std::endl;
                std::cout << " True taxon: " << true_taxid << " gene: " << true_geneid << std::endl;
                for (auto &alignment: alignments) {
                    std::cout << alignment.ToString() << std::endl;
                }
                std::cout << "\nAnchors:" << std::endl;
                for (auto &anchor: anchors) {
                    std::cout << anchor.a.ToString() << " " << anchor.b.ToString() << " " << anchor.hit_anchor_count
                              << std::endl;
                }
                std::cout << "\nAlignments:" << std::endl;
                for (auto &alignment: alignments) {
                    std::cout << alignment.ToString() << std::endl;
                }

                std::cout << "alignments " << alignments.size() << " / " << anchors.size() << " anchors" << std::endl;
                std::cout << std::string(50, '-') << " END" << std::endl;
            }
            m_alignment_bce.fn++;
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
            m_anchor_bce.Join(other.m_anchor_bce);
            m_alignment_bce.Join(other.m_alignment_bce);
            m_best_alignment_bce.Join(other.m_best_alignment_bce);
            m_anchor_fp_no_hit += other.m_anchor_fp_no_hit;
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

