//
// Created by fritsche on 27/09/22.
//

#pragma once

#include <vector>
#include "BinaryClassifierEvaluator.h"

//#include "SeedingStrategy.h"
#include "AlignmentOutputHandler.h"
#include "Constants.h"
#include "ChainingStrategy.h"

namespace protal {
    using BCE = BinaryClassifierEvaluator;
    class CoreBenchmark {
        BCE m_seed_bce{"Seeds"};
        BCE m_anchor_bce{"Anchors"};
        BCE m_alignment_bce{"Alignments"};
        BCE m_best_alignment_bce{"Best_alignment"};
        BCE m_best_unique_alignment_bce{"Best_unique_alignment"};
        BCE m_alignment_pe_bce{"Alignments_PE"};
        BCE m_best_alignment_pe_bce{"Best_alignment_PE"};

        size_t m_anchor_fp_no_hit = 0;

        std::ofstream* ofs = nullptr;

    public:
        CoreBenchmark() {};

        CoreBenchmark(CoreBenchmark const& other) :
                ofs(other.ofs) {};

        void SetOutput(std::string const& output_file) {
            ofs = new ofstream(output_file, std::ios::app);
        }

        void DestroyOutput() {
            if (ofs) {
                if (ofs->is_open()) ofs->close();
                delete ofs;
            }
        }

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
            bool found = false;
            for (auto& seed : seeds) {
                if (seed.taxid == true_taxid && seed.geneid == true_geneid) {
                    m_seed_bce.tp++;
                    found = true;
                } else {
                    m_seed_bce.fp++;
                }
            }

            m_seed_bce.fn += !found;
        };

        void AddAnchors(AlignmentAnchorList const& anchors, size_t true_taxid, size_t true_geneid) {
            bool found = false;
            for (auto& anchor : anchors) {
                if (anchor.taxid == true_taxid && anchor.geneid == true_geneid) {
                    m_anchor_bce.tp++;
                    found = true;
                } else {
                    m_anchor_bce.fp++;
                }
            }

            m_anchor_fp_no_hit += anchors.empty();

            m_anchor_bce.fn += !found;
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
                m_best_alignment_bce.fn++;
            }

            if (alignments.size() == 1 || (alignments.size() > 1 && alignments[1].AlignmentScore() < best_alignment.AlignmentScore())) {

                if (best_alignment.Taxid() == true_taxid && best_alignment.GeneId() == true_geneid) {
                    m_best_unique_alignment_bce.tp++;
                } else {
                    m_best_unique_alignment_bce.fp++;
                    m_best_unique_alignment_bce.fn++;
                }
            }

            bool found = false;
            for (auto& alignment : alignments) {
                if (alignment.Taxid() == true_taxid && alignment.GeneId() == true_geneid) {
                    m_alignment_bce.tp++;
                    found = true;
                } else {
                    m_alignment_bce.fp++;
                }
            }
            m_alignment_bce.fn += !found;
        };

        void AddPairedAlignmentResults(PairedAlignmentResultList const& alignments, size_t true_taxid, size_t true_geneid) {
            if (alignments.empty()) {
//                m_alignment_pe_bce.fn++;
                m_best_alignment_pe_bce.fn++;
                return;
            }

            auto& best_alignment = alignments[0];
            auto best_taxid = best_alignment.first.IsSet() ? best_alignment.first.Taxid() : best_alignment.second.Taxid();
            auto best_geneid = best_alignment.first.IsSet() ? best_alignment.first.GeneId() : best_alignment.second.GeneId();
            if (best_taxid == true_taxid && best_geneid == true_geneid) {
                m_best_alignment_pe_bce.tp++;
            } else {
                m_best_alignment_pe_bce.fp++;
                m_best_alignment_pe_bce.fn++;
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
            m_best_alignment_pe_bce.Join(other.m_best_alignment_pe_bce);
            m_anchor_fp_no_hit += other.m_anchor_fp_no_hit;
        }

        void operator() (SeedList const& seeds, AlignmentAnchorList const& anchors, AlignmentResultList const& alignments, std::string const& header) {
            auto [tax_id_truth, gene_id_truth] = CoreBenchmark::HeaderToTruth(header);

            AddSeeds(seeds, tax_id_truth, gene_id_truth);
            AddAnchors(anchors, tax_id_truth, gene_id_truth);
            AddAlignmentResults(alignments, tax_id_truth, gene_id_truth);

        }

        void operator() (SeedList const& seeds, AlignmentAnchorList const& anchors, AlignmentResultList const& alignments, PairedAlignmentResultList const& alignment_pairs, std::string const& header) {
            auto [tax_id_truth, gene_id_truth] = CoreBenchmark::HeaderToTruth(header);

            AddSeeds(seeds, tax_id_truth, gene_id_truth);
            AddAnchors(anchors, tax_id_truth, gene_id_truth);
            AddAlignmentResults(alignments, tax_id_truth, gene_id_truth);
            AddPairedAlignmentResults(alignment_pairs, tax_id_truth, gene_id_truth);

        }

        void WriteRowStats(std::string prefix="") {
            m_seed_bce.WriteRowHeader(std::cout, prefix);
            m_seed_bce.WriteRowStats(std::cout, prefix);
            m_anchor_bce.WriteRowStats(std::cout, prefix);
            m_alignment_bce.WriteRowStats(std::cout, prefix);
            m_best_alignment_bce.WriteRowStats(std::cout, prefix);
            m_best_unique_alignment_bce.WriteRowStats(std::cout, prefix);
            m_best_alignment_pe_bce.WriteRowStats(std::cout, prefix);

            if (ofs) {
                if (ofs->tellp() == 0) {
                    m_seed_bce.WriteRowHeader(*ofs, prefix);
                }
                m_seed_bce.WriteRowStats(*ofs, prefix);
                m_anchor_bce.WriteRowStats(*ofs, prefix);
                m_alignment_bce.WriteRowStats(*ofs, prefix);
                m_best_alignment_bce.WriteRowStats(*ofs, prefix);
                m_best_unique_alignment_bce.WriteRowStats(*ofs, prefix);
            }
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

