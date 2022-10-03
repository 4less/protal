//
// Created by fritsche on 05/09/22.
//

#pragma once

#include "WFA2Wrapper.h"
#include "AlignmentOutputHandler.h"
#include "GenomeLoader.h"
#include "SeedingStrategy.h"
#include "KmerUtils.h"
#include "FastAlignment.h"

namespace protal {
    class SimpleAlignmentHandler {
    private:
        WFA2Wrapper m_aligner;
        AlignmentResult m_alignment_result;
        GenomeLoader& m_genome_loader;
        size_t m_kmer_size = 15;


        static bool IsReverse(LookupResult const& first_anchor, LookupResult const& second_anchor) {
            return first_anchor.readpos > second_anchor.readpos;
        }

        inline void ReverseAnchorReadPos(LookupResult& anchor, size_t& read_length) {
            anchor.readpos = read_length - anchor.readpos - m_kmer_size;
        }

        inline void ReverseAnchorPairReadPos(LookupResult& anchor_a, LookupResult& anchor_b, size_t& read_length) {
            ReverseAnchorReadPos(anchor_a, read_length);
            ReverseAnchorReadPos(anchor_b, read_length);
        }

    public:
        size_t total_alignments = 0;
        size_t alignments_ani_lower93 = 0;
        size_t dummy = 0;
        Benchmark bm_alignment{ "Alignment" };


        SimpleAlignmentHandler(GenomeLoader& genome_loader, WFA2Wrapper& aligner, size_t kmer_size) :
                m_genome_loader(genome_loader),
                m_aligner(aligner),
                m_kmer_size(kmer_size) {

        }

        void operator() (AlignmentAnchorList& anchors, AlignmentResultList& results, std::string& sequence) {

            constexpr bool alignment_verbose = false;

            auto& fwd = sequence;
            auto rev = KmerUtils::ReverseComplement(fwd);

            char* read = nullptr;
            size_t read_len = sequence.length();
            bool reversed = false;

            std::sort(anchors.begin(), anchors.end(), [](AlignmentAnchor const& a, AlignmentAnchor const& b) {
                return a.hit_anchor_count > b.hit_anchor_count;
            });

            int take_top = 3;
//            int take_top = 10;
            int last_score = 0;

            for (auto& anchor : anchors) {

                if (IsReverse(anchor.a, anchor.b)) {
                    reversed = true;
                    ReverseAnchorPairReadPos(anchor.a, anchor.b, read_len);
                    read = const_cast<char *>(rev.c_str());
                } else {
                    read = const_cast<char *>(sequence.c_str());
                }

                int rpos_diff = static_cast<int>(anchor.b.readpos) - static_cast<int>(anchor.a.readpos);
                int gpos_diff = static_cast<int>(anchor.b.genepos) - static_cast<int>(anchor.a.genepos);

                int64_t anchor_distance_deviation = rpos_diff - gpos_diff;

                auto& genome = m_genome_loader.GetGenome(anchor.a.taxid);
                auto& gene = genome.GetGeneOMP(anchor.a.geneid);

                int abs_pos = anchor.a.genepos - anchor.a.readpos;
                dummy += abs_pos;


                size_t read_start = abs_pos < 0 ? -abs_pos : 0;
                size_t gene_start = std::max(abs_pos, 0);

                size_t overlap = std::min(gene.Sequence().length() - gene_start, read_len - read_start);

                total_alignments++;
                std::string query = std::string(read + read_start, overlap);
                std::string reference = gene.Sequence().substr(gene_start, overlap);
                dummy += query.length();



                double cigar_ani = 0;
                std::string cigar = "";
                int score = -1000;

                bool approximate_alignment = anchor_distance_deviation == 0;
//                bool approximate_alignment = true;

                bm_alignment.Start();
                if (approximate_alignment) {
                    auto [ascore, acigar] = FastAligner::FastAlign(query, reference);
                    cigar_ani = WFA2Wrapper::CigarANI(acigar);
                    cigar = acigar;
                    score = ascore;
//
//                    std::cout << "------------------------------------------" << std::endl;
//                    std::cout << anchor.a.ToString() << std::endl;
//                    std::cout << anchor.b.ToString() << std::endl;
//
//                    auto start = anchor.a.genepos;
//                    auto end = anchor.b.genepos;
//
//
//                    m_aligner.PrintAlignment();
//                    std::cout << "-----------------" << std::endl;
//                    std::cout << "Score: " << ascore << std::endl;
//                    std::cout << "Cigar: " << acigar << std::endl;
//                    std::cout << "-----------------" << std::endl;
//                    std::cout << cigar_ani << " vs " << cigar_ani2 << std::endl;
//
//                    if ((cigar_ani >= 0.93 || cigar_ani2 >= 0.93) && cigar != m_aligner.GetAligner().getAlignmentCigar()) {
//                        Utils::Input();
//                    }
                } else {
                    Benchmark bm_local{"alignment"};
                    m_aligner.Alignment(query, reference);
                    bm_local.Stop();

                    cigar_ani = WFA2Wrapper::CigarANI(m_aligner.GetAligner().getAlignmentCigar());
//                    std::cout << "anchor_distance_deviation: " << anchor_distance_deviation << std::endl;

                    cigar = m_aligner.GetAligner().getAlignmentCigar();
                    score = m_aligner.GetAligner().getAlignmentScore();

                    if constexpr(alignment_verbose) {
#pragma omp critical (alignment_debug)
                        {
                            std::cerr << "Alignment\t" << bm_local.GetDuration(Time::microseconds);
                            std::cerr << '\t' << query.length() << '\t' << reference.length() << '\t';
                            std::cerr << cigar_ani << '\t' << anchor.hit_anchor_count << '\n';
                        }
                    }
                }
                bm_alignment.Stop();


                m_alignment_result.Set(anchor.a.taxid, anchor.b.geneid, abs_pos, !reversed);
                m_alignment_result.Set(score, cigar);


                results.emplace_back(m_alignment_result);


                if (--take_top <= 0) {// && score < last_score) {
                    break;
                }
                last_score = score;
            }

            // Sort alignment results
            std::sort(results.begin(), results.end(), [](AlignmentResult const& a, AlignmentResult const& b) {
                return a.AlignmentScore() > b.AlignmentScore();
            });
        }
    };
}