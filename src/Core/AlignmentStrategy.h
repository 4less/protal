//
// Created by fritsche on 05/09/22.
//

#pragma once

#include "WFA2Wrapper.h"
#include "AlignmentOutputHandler.h"
#include "GenomeLoader.h"
#include "SeedingStrategy.h"
#include "KmerUtils.h"

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

            size_t take_top = 3;
            for (auto& anchor : anchors) {

                if (take_top-- == 0) break;

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

                bm_alignment.Start();
                Benchmark bm_local{"alignment"};
                m_aligner.Alignment(query, reference);
                bm_local.Stop();
                bm_alignment.Stop();

                auto cigar_ani = WFA2Wrapper::CigarANI(m_aligner.GetAligner().getAlignmentCigar());

                if constexpr(alignment_verbose) {
#pragma omp critical (alignment_debug)
                    {
                        std::cerr << "Alignment\t" << bm_local.GetDuration(Time::microseconds);
                        std::cerr << '\t' << query.length() << '\t' << reference.length() << '\t';
                        std::cerr << cigar_ani << '\t' << anchor.hit_anchor_count << '\n';

//                    double min_ani = static_cast<double>(anchor.hit_anchor_count)/reference.length();
//                    if (cigar_ani < min_ani) {
//                        std::cout << "Reversed: " << reversed << std::endl;
//                        std::cout << "Kmersize: " << m_kmer_size << std::endl;
//                        std::cout << "Gene start: " << gene_start << " - len: " << overlap << std::endl;
//                        std::cout << anchor.a.ToString() << " " << anchor.b.ToString() << std::endl;
//                        ReverseAnchorPairReadPos(anchor.a, anchor.b, read_len);
//                        std::cout << "Reversed: ";
//                        std::cout << anchor.a.ToString() << " " << anchor.b.ToString() << std::endl;
//                        std::cout << "min ani: " << min_ani << std::endl;
//                        m_aligner.PrintAlignment();
//                        Utils::Input();
//                    }
                    }
                }

                m_alignment_result.Set(anchor.a.taxid, anchor.b.geneid, abs_pos, !reversed);
                m_alignment_result.Set(m_aligner.GetAligner().getAlignmentScore(), m_aligner.GetAligner().getAlignmentCigar());


                if (cigar_ani > 0.93) {
                    results.emplace_back(m_alignment_result);
//
//                    std::cout << "CigarANI: " << cigar_ani << std::endl;
//                    m_aligner.PrintAlignment();
//                    Utils::Input();
                } else {
//                    int distance = (static_cast<int>(anchor.b.genepos) - anchor.a.genepos) - (static_cast<int>(anchor.b.readpos) - anchor.a.readpos);
//                    std::cout << "_______________________" << std::endl;
//                    std::cout << "Anchors: " << anchors.size() << '\t' << cigar_ani << '\t' << distance << "\trev? " << reversed << "\t\t" << anchor.a.ToString() << '\t' << anchor.b.ToString() << std::endl;
//                    std::cout << "___________Other anchors____________" << std::endl;
//                    for (auto& anchor : anchors) {
//                        std::cout << anchor.a.ToString() << " " << anchor.b.ToString() << std::endl;
//                    }
//                    std::cout << "____________Alignment___________" << std::endl;
//                    m_aligner.PrintAlignment();
                    alignments_ani_lower93++;
                }
            }
        }
    };
}