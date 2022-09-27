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

        void AddAnchors(AlignmentAnchorList& anchors, size_t true_taxid, size_t true_geneid) {
            for (auto& anchor : anchors) {
                if (anchor.a.taxid == true_taxid && anchor.a.geneid == true_geneid) {
                    m_anchor_bce.tp++;
                    return;
                }
            }
            m_anchor_bce.fn;
        };

        void AddAlignmentResults(AlignmentResultList& alignments, size_t true_taxid, size_t true_geneid) {
            for (auto& alignment : alignments) {
                if (alignment.Taxid() == true_taxid && alignment.GeneId() == true_geneid) {
                    m_alignment_bce.tp++;
                    return;
                }
            }
            m_alignment_bce.fn;
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
        }
    };
}

