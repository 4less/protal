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
        BCE m_best_anchor_bce{"Best_Anchor"};
        BCE m_alignment_bce{"Alignments"};
        BCE m_best_alignment_bce{"Best_alignment"};
        BCE m_best_unique_alignment_bce{"Best_unique_alignment"};
        BCE m_alignment_pe_bce{"Alignments_PE"};
        BCE m_best_alignment_pe_bce{"Best_alignment_PE"};

        size_t m_anchor_fp_no_hit = 0;

        std::ofstream* ofs = nullptr;

        size_t m_total_reads = 0;
        size_t m_total_pairs = 0;
        size_t m_correct_without_pair = 0;

    public:
        size_t equally_good_best_count = 0;

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

            if (!anchors.empty()) {
                auto& best = anchors[0];
                if (best.taxid == true_taxid && best.geneid == true_geneid) {
                    m_best_anchor_bce.tp++;
                } else {
                    m_best_anchor_bce.fp++;
                }
            }

            m_anchor_fp_no_hit += anchors.empty();
            m_best_anchor_bce.fn += anchors.empty();
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

        void ErrorOutput(SeedList const& seeds1, SeedList const& seeds2, AlignmentAnchorList const& anchors1, AlignmentAnchorList const& anchors2, AlignmentResultList const& alignment1, AlignmentResultList const& alignment2, PairedAlignmentResultList const& alignments, FastxRecord& record1, FastxRecord& record2) {
            auto [tax_id_truth, gene_id_truth] = CoreBenchmark::HeaderToTruth(record1.id);
            size_t true_taxid = tax_id_truth;
            size_t true_geneid = gene_id_truth;

            if (alignments.empty()) {
                std::cerr << "Seeds1:" << std::endl;
                for (auto& seed : seeds1) {
                    std::cerr << seed.ToString() << std::endl;
                }
                std::cerr << "Seeds2:" << std::endl;
                for (auto& seed : seeds2) {
                    std::cerr << seed.ToString() << std::endl;
                }
                std::cerr << "Anchors1:" << std::endl;
                for (auto& anchor : anchors1) {
                    std::cerr << anchor.ToString() << std::endl;
                }
                std::cerr << "Anchors2:" << std::endl;
                for (auto& anchor : anchors2) {
                    std::cerr << anchor.ToString() << std::endl;
                }
                return;
            }
            auto& best_alignment = alignments[0];
            auto best_taxid = best_alignment.first.IsSet() ? best_alignment.first.Taxid() : best_alignment.second.Taxid();
            auto best_geneid = best_alignment.first.IsSet() ? best_alignment.first.GeneId() : best_alignment.second.GeneId();

            bool equally_good_best = alignments.size() > 1 && Bitscore(alignments[0]) == Bitscore(alignments[1]);
            equally_good_best_count += equally_good_best;

            if (!equally_good_best && (best_taxid != true_taxid || best_geneid != true_geneid)) {
                std::cerr << "Truth: " << true_taxid << " Gene: " << true_geneid << std::endl;
                std::cerr << "Best : " << best_taxid << " Gene: " << best_geneid << std::endl;
                std::cerr << "Seeds1:" << std::endl;
                for (auto& seed : seeds1) {
                    std::cerr << seed.ToString() << std::endl;
                }
                std::cerr << "Seeds2:" << std::endl;
                for (auto& seed : seeds2) {
                    std::cerr << seed.ToString() << std::endl;
                }
                std::cerr << "Anchors1:" << std::endl;
                for (auto& anchor : anchors1) {
                    std::cerr << anchor.ToString() << std::endl;
                }
                std::cerr << "Anchors2:" << std::endl;
                for (auto& anchor : anchors2) {
                    std::cerr << anchor.ToString() << std::endl;
                }
                std::cerr << "Alignments1:" << std::endl;
                for (auto& alignment : alignment1) {
                    std::cerr << alignment.ToString() << std::endl;
                }
                std::cerr << "Alignments2:" << std::endl;
                for (auto& alignment : alignment2) {
                    std::cerr << alignment.ToString() << std::endl;
                }
                std::cerr << "Paired Alignments:" << std::endl;
                std::cerr << "Truth: " << true_taxid << " Gene: " << true_geneid << std::endl;
                std::cerr << "Best : " << best_taxid << " Gene: " << best_geneid << std::endl;
                for (auto& [ar1, ar2] : alignments) {
                    std::cerr << "BitScore: " << Bitscore({ ar1,  ar2}) << std::endl;
                    if (ar1.IsSet()) {
                        std::cerr << "1: " << ar1.ToString() << std::endl;
                    }
                    if (ar2.IsSet()) {
                        std::cerr << "2: " << ar2.ToString() << std::endl;
                    }
                    std::cerr << "-----" << std::endl;
                }
            }
        }

        void AddPairedAlignmentResults(PairedAlignmentResultList const& alignments, size_t true_taxid, size_t true_geneid) {
            if (alignments.empty()) {
                m_alignment_pe_bce.fn++;
                m_best_alignment_pe_bce.fn++;
                return;
            }

            bool equally_good_best = alignments.size() > 1 && Bitscore(alignments[0]) == Bitscore(alignments[1]);

            auto& best_alignment = alignments[0];
            auto best_taxid = best_alignment.first.IsSet() ? best_alignment.first.Taxid() : best_alignment.second.Taxid();
            auto best_geneid = best_alignment.first.IsSet() ? best_alignment.first.GeneId() : best_alignment.second.GeneId();
            if (equally_good_best) {
                m_best_alignment_pe_bce.tn++;
            } else if (best_taxid == true_taxid && best_geneid == true_geneid) {
                m_best_alignment_pe_bce.tp++;
            } else {
                m_best_alignment_pe_bce.fp++;
                m_best_alignment_pe_bce.fn++;
            }

            bool found = false;
            for (auto& [ar1, ar2] : alignments) {
                if ((ar1.IsSet() && ar1.Taxid() == true_taxid && ar1.GeneId() == true_geneid) &&
                    (ar1.IsSet() && ar1.Taxid() == true_taxid && ar1.GeneId() == true_geneid)) {
                    if (!found && !(ar1.IsSet() && ar2.IsSet())) {
                        m_correct_without_pair++;
                    }
                    found = true;
                    m_alignment_pe_bce.tp++;
                } else {
                    m_alignment_pe_bce.fp++;
                }
            }
            m_alignment_pe_bce.fn += !found;
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
            m_best_anchor_bce.Join(other.m_best_anchor_bce);
            m_alignment_bce.Join(other.m_alignment_bce);
            m_best_alignment_bce.Join(other.m_best_alignment_bce);
            m_best_unique_alignment_bce.Join(other.m_best_unique_alignment_bce);
            m_best_alignment_pe_bce.Join(other.m_best_alignment_pe_bce);
            m_alignment_pe_bce.Join(other.m_alignment_pe_bce);
            m_anchor_fp_no_hit += other.m_anchor_fp_no_hit;
            m_total_reads += other.m_total_reads;
            m_total_pairs += other.m_total_pairs;
            m_correct_without_pair += other.m_correct_without_pair;
            equally_good_best_count += other.equally_good_best_count;
        }

        void operator() (SeedList const& seeds, AlignmentAnchorList const& anchors, AlignmentResultList const& alignments, std::string const& header) {
            auto [tax_id_truth, gene_id_truth] = CoreBenchmark::HeaderToTruth(header);

            AddSeeds(seeds, tax_id_truth, gene_id_truth);
            AddAnchors(anchors, tax_id_truth, gene_id_truth);
            AddAlignmentResults(alignments, tax_id_truth, gene_id_truth);

            m_total_reads++;
        }

        void operator() (SeedList const& seeds, AlignmentAnchorList const& anchors, AlignmentResultList const& alignments, PairedAlignmentResultList const& alignment_pairs, std::string const& header) {
            auto [tax_id_truth, gene_id_truth] = CoreBenchmark::HeaderToTruth(header);

            AddSeeds(seeds, tax_id_truth, gene_id_truth);
            AddAnchors(anchors, tax_id_truth, gene_id_truth);
            AddAlignmentResults(alignments, tax_id_truth, gene_id_truth);
            AddPairedAlignmentResults(alignment_pairs, tax_id_truth, gene_id_truth);

            m_total_reads++;
            m_total_pairs++;
        }

        void WriteRowStatsToStream(std::ostream& os, std::string prefix="") {
            m_seed_bce.WriteRowStats(os, prefix);
            m_anchor_bce.WriteRowStats(os, prefix);
            m_best_anchor_bce.WriteRowStats(os, prefix);
            m_alignment_bce.WriteRowStats(os, prefix);
            m_best_alignment_bce.WriteRowStats(os, prefix);
            m_best_unique_alignment_bce.WriteRowStats(os, prefix);
            m_alignment_pe_bce.WriteRowStats(os, prefix);
            m_best_alignment_pe_bce.WriteRowStats(os, prefix);
        }

        void WriteRowStatsToFile(std::string file_path) {
            std::ofstream out(file_path, std::ios::out);
            WriteRowStatsToStream(out);
            out.close();
        }

        void WriteRowStats(std::string prefix="") {
            m_seed_bce.WriteRowHeader(std::cout, prefix);
            m_seed_bce.WriteRowStats(std::cout, prefix);
            m_anchor_bce.WriteRowStats(std::cout, prefix);
            m_best_anchor_bce.WriteRowStats(std::cout, prefix);
            m_alignment_bce.WriteRowStats(std::cout, prefix);
            m_best_alignment_bce.WriteRowStats(std::cout, prefix);
            m_best_unique_alignment_bce.WriteRowStats(std::cout, prefix);
            m_alignment_pe_bce.WriteRowStats(std::cout, prefix);
            m_best_alignment_pe_bce.WriteRowStats(std::cout, prefix);
            std::cout << "Total reads:                 " << m_total_reads << std::endl;
            std::cout << "Total pairs:                 " << m_total_pairs << std::endl;
            std::cout << "Total best PE is not paired: " << m_correct_without_pair << std::endl;
            std::cout << "Equally good best:           " << equally_good_best_count << std::endl;
            std::cout << "Equally good best rate:      " << static_cast<double>(equally_good_best_count)/m_total_pairs << std::endl;
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

