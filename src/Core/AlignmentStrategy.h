//
// Created by fritsche on 05/09/22.
//

#pragma once

#include "WFA2Wrapper.h"
#include "AlignmentOutputHandler.h"
#include "GenomeLoader.h"
#include "SeedingStrategy.h"
#include "ChainingStrategy.h"
#include "KmerUtils.h"
#include "FastAlignment.h"

namespace protal {
    class SimpleAlignmentHandler {
    private:
        WFA2Wrapper m_aligner;
        AlignmentResult m_alignment_result;
        GenomeLoader& m_genome_loader;
        size_t m_kmer_size = 15;

        size_t m_align_top = 3;
        double m_max_score_ani = 0;

        bool m_fastalign = false;


        static bool IsReverse(LookupResult const& first_anchor, LookupResult const& second_anchor) {
            return first_anchor.readpos > second_anchor.readpos;
        }

        inline void ReverseAnchorReadPos(LookupResult& anchor, size_t& read_length) {
            anchor.readpos = read_length - anchor.readpos - m_kmer_size;
        }

        inline void ReverseAnchorReadPos(ChainLink& seed, size_t& read_length) {
            seed.readpos = read_length - seed.readpos - seed.length;
        }

        inline void ReverseAnchorPairReadPos(LookupResult& anchor_a, LookupResult& anchor_b, size_t& read_length) {
            ReverseAnchorReadPos(anchor_a, read_length);
            ReverseAnchorReadPos(anchor_b, read_length);
        }

    public:
        size_t total_alignments = 0;
        size_t total_tail_alignments = 0;
        size_t total_tail_length = 0;
        Benchmark bm_alignment{ "Alignment" };
        Benchmark bm_seedext{ "Seed Extension" };
        size_t dummy = 0;


        SimpleAlignmentHandler(GenomeLoader& genome_loader, WFA2Wrapper& aligner, size_t kmer_size, size_t align_top, double max_score_ani, bool fastalign) :
                m_genome_loader(genome_loader),
                m_aligner(aligner),
                m_kmer_size(kmer_size),
                m_align_top(align_top),
                m_max_score_ani(max_score_ani),
                m_fastalign(fastalign) {};

        SimpleAlignmentHandler(SimpleAlignmentHandler const& other) :
                m_genome_loader(other.m_genome_loader),
                m_aligner(other.m_aligner),
                m_kmer_size(other.m_kmer_size),
                m_align_top(other.m_align_top),
                m_max_score_ani(other.m_max_score_ani),
                m_fastalign(other.m_fastalign) {};



        inline bool ReverseAnchor(AlignmentAnchor& anchor, size_t read_len) {
            bool reversed = false;
            if (IsReverse(anchor.a, anchor.b)) {
                reversed = true;
                ReverseAnchorPairReadPos(anchor.a, anchor.b, read_len);
            }
            return reversed;
        }

        inline void ReverseAnchor(CAlignmentAnchor& anchor, size_t read_len) {
            for (auto& chainlink : anchor.chain) {
                ReverseAnchorReadPos(chainlink, read_len);
            }
            std::reverse(anchor.chain.begin(), anchor.chain.end());
        }

        inline int SeedIndel(CAlignmentAnchor& anchor) {
            if (anchor.chain.size() == 1) return 0;
            int rpos_diff = static_cast<int>(anchor.Back().readpos) - static_cast<int>(anchor.Front().readpos);
            int gpos_diff = static_cast<int>(anchor.Back().genepos) - static_cast<int>(anchor.Front().genepos);
            return rpos_diff - gpos_diff;
        }



        size_t ExtendSeedLeft(Seed const& s, std::string const& query, std::string const& gene) {
            size_t extension = 0;
            size_t max_extension_len = std::min(s.readpos, s.genepos);
            for (int qpos = s.readpos, rpos = s.genepos;
                 extension < max_extension_len && query[qpos] == gene[rpos];
                 qpos--, rpos--) {
                extension++;
            }
            return extension;
        }

        size_t ExtendSeedRight(Seed const& s, size_t k, std::string const& query, std::string const& gene) {
            size_t extension = 0;
            size_t max_extension_len = std::min(
                    static_cast<uint32_t>(query.length() - s.readpos - k),
                    static_cast<uint32_t>(gene.length() - s.genepos - k));
            for (int qpos = s.readpos + k, rpos = s.genepos + k;
                 extension < max_extension_len && query[qpos] == gene[rpos];
                 qpos++, rpos++) {
                extension++;
            }
            return extension;
        }

        static std::pair<std::string_view, std::string_view>
        GetViewsForAlignment(Seed const& a, std::string const& query, std::string const& target) {

            auto min = std::min(a.readpos, a.genepos);
            auto qstart = a.readpos - min;
            auto tstart = a.genepos - min;
            auto overlap_length = std::min(query.length() - qstart, target.length() - tstart);

            if (tstart + overlap_length > target.length()) exit(9);
            if (qstart + overlap_length > query.length()) exit(10);

            return {
                std::string_view(query.c_str() + qstart, overlap_length),
                std::string_view(target.c_str() + tstart, overlap_length)
            };
        }

        static std::pair<size_t, size_t> ExtendSeed(Seed const& s, size_t k, std::string const& query, std::string const& gene) {
            size_t extension_left = 0;
            size_t extension_right = 0;

            size_t max_extension_len = std::min(s.readpos, s.genepos);
            for (int qpos = s.readpos - 1, rpos = s.genepos - 1;
                 extension_left < max_extension_len && query[qpos] == gene[rpos];
                 qpos--, rpos--) {
                extension_left++;
            }
            max_extension_len = std::min(
                    static_cast<uint32_t>(query.length() - s.readpos - k),
                    static_cast<uint32_t>(gene.length() - s.genepos - k));
            for (int qpos = s.readpos + k, rpos = s.genepos + k;
                 extension_right < max_extension_len && query[qpos] == gene[rpos];
                 qpos++, rpos++) {
                extension_right++;
            }
            return { extension_left, extension_right };
        }

        static std::pair<size_t, size_t> ExtendSeed(ChainLink& s, std::string const& query, std::string const& gene) {
            size_t extension_left = 0;
            size_t extension_right = 0;

            size_t max_extension_len = std::min(static_cast<uint32_t>(s.readpos), s.genepos);
            for (int qpos = s.readpos - 1, rpos = s.genepos - 1;
//                 qpos >= 0 && rpos >= 0 && // Maybe remove
                 extension_left < max_extension_len && query[qpos] == gene[rpos];
                 qpos--, rpos--) {
                extension_left++;
            }
            max_extension_len = std::min(
                    static_cast<uint32_t>(query.length() - s.readpos - s.length),
                    static_cast<uint32_t>(gene.length() - s.genepos - s.length));

            for (int qpos = s.readpos + s.length, rpos = s.genepos + s.length;
//                 qpos < query.length() && rpos < gene.length() &&
                 extension_right < max_extension_len && query[qpos] == gene[rpos];
                 qpos++, rpos++) {
                extension_right++;
            }
            assert(s.ReadStart() >= extension_left);
            assert(s.GeneStart() >= extension_left);
            assert(s.ReadEnd() + extension_right <= query.length());
            assert(s.GeneEnd() + extension_right <= gene.length());

            s.ExtendLength(extension_left, extension_right);
            return { extension_left, extension_right };
        }


        static int MaxScore(double min_ani, uint32_t overlap_length, int mismatch_penalty=4) {
            return (std::ceil((1 - min_ani) * static_cast<double>(overlap_length)) * mismatch_penalty) + 1;
        }



        void ExtendAllAnchors(AlignmentAnchorList const& anchors, std::string const& fwd, std::string const& rev) {
            for (auto anchor : anchors) {
                // copy anchor;
                bool reversed = !anchor.forward;
                auto query = reversed ? rev : fwd;

                // Get Resources
                auto& genome = m_genome_loader.GetGenome(anchor.taxid);
                auto& gene = genome.GetGeneOMP(anchor.geneid);


                for (auto i = 0; i < anchor.chain.size(); i++) {
                    auto& seed = anchor.chain[i];
                    auto [lefta, righta] = ExtendSeed(seed, query, gene.Sequence());
                    dummy += lefta + righta;
//                    if (i > 0 && seed.OverlapsWithLeft(anchor.chain[i-1])) {
//                        seed.Merge(anchor.chain[i-1].readpos, anchor.chain[i-1].length);
//                        anchor.chain.erase(anchor.chain.begin() + i);
//                        i--;
//                    }
                }
            }
        }

        inline std::tuple<bool, int, std::string> Align(CAlignmentAnchor& anchor, std::string const& query, std::string const& gene) {
            //                front_end_q                         back_start_q
            //                v                                   v
            //----------------SSSSSSSSS------------------SSSSSSSSS-------------------- Query
            //----------------SSSSSSSSS------------------SSSSSSSSS-------------------- Ref
            //                ^                                   ^
            //                front_end_r                         back_start_r

            constexpr bool debug = false;
            // Phase one = Extend left seed

//            std::cout << "---------------------------------------------" << std::endl;
//            std::cout << anchor.ToString() << std::endl;

            bool twoseeds = anchor.chain.size() > 1;

            assert(!anchor.chain.empty());
            auto [lefta, righta] = ExtendSeed(anchor.Front(), query, gene);

            if (twoseeds) {
                ExtendSeed(anchor.Back(), query, gene);
                bool overlap = false;
                for (auto i = 1; i < anchor.chain.size(); i++) {
                    if (anchor.chain[i].OverlapsWithLeft(anchor.chain[i-1])) {
//                        std::cout << "Merge " << anchor.chain[i-1].ToString() << " <<<< " << anchor.chain[i].ToString() << std::endl;
                        anchor.chain[i-1].Merge(anchor.chain[i].readpos, anchor.chain[i].length);
                        anchor.chain.erase(anchor.chain.begin() + i);
                        i--;
                        overlap = true;
                    }
                }

//                if (overlap) {
//                    std::cout << anchor.ToString() << std::endl;
//                    std::cout << anchor.Front().ReadStart() << " - " << anchor.Back().ReadEnd() << std::endl;
//                    m_aligner.Alignment(
//                            std::string_view(query.c_str() + anchor.Front().ReadStart(), anchor.Back().ReadEnd() - anchor.Front().ReadStart() ),
//                            std::string_view(gene.c_str() + anchor.Front().genepos, anchor.Back().ReadEnd() - anchor.Front().ReadStart())
//                            );
//                    m_aligner.PrintAlignment();
//
//                    Utils::Input();
//                }
            }



//            std::cout << "Extended: " << std::endl;
//            std::cout << anchor.ToString() << std::endl;

            twoseeds = anchor.chain.size() > 1;

            std::string mid_cigar = "";
            int mid_score = 0;
            int max_score = MaxScore(m_max_score_ani, anchor.UpdateLength(), 4);

            auto& first = anchor.Front();
            auto& second = anchor.Back();
            if (twoseeds) {
                auto gap = second.ReadStart() - first.ReadEnd();
                if (gap == 1) {
                    mid_cigar = 'X';
                    mid_score = -4;
                } else {
                    auto [score, cigar] = FastAligner::FastAlign(
                            std::string_view(query.c_str() + first.ReadEnd(), gap),
                            std::string_view(gene.c_str() + first.GeneEnd(), gap));
                    mid_score = score;
                    mid_cigar = cigar;
                }
                mid_cigar = std::string(first.length, 'M') + mid_cigar + std::string(second.length, 'M');
            } else {
                mid_cigar = std::string(first.length, 'M');
            }



            if constexpr(debug) {
                std::cout << "Total seed coverage " << anchor.total_length << std::endl;
                std::cout << "mid_score " << mid_score << std::endl;
            }

            max_score += mid_score;

            if (max_score < 0) {
                return { false, -1, "" };
            }


            // create stringviews
            auto left_length = std::min(first.ReadStart(), first.GeneStart());
            auto right_length = std::min(query.length() - second.ReadEnd(), gene.length() - second.GeneEnd());

            auto query_tail_left = std::string_view(query.c_str() + first.ReadStart() - left_length, left_length);
            auto target_tail_left = std::string_view(gene.c_str() + first.GeneStart() - left_length, left_length);
            auto query_tail_right = std::string_view(query.c_str() + second.ReadEnd(), right_length);
            auto target_tail_right = std::string_view(gene.c_str() + second.GeneEnd(), right_length);

            std::string lcig = "";
            std::string rcig = "";
            size_t lsco = 0;
            size_t rsco = 0;

            if (query_tail_left.length() + query_tail_right.length() + mid_cigar.length() != (second.ReadEnd() + right_length) - (first.ReadStart() - left_length)) {
                std::cout << "query length:     " << query.length() << std::endl;
                std::cout << "GeneLength:       " << gene.length() << std::endl;
                std::cout << "mid_len:          " << mid_cigar.length() << std::endl;
                std::cout << "Anchor:           " << anchor.Front().ToString() << " " << anchor.Back().ToString() << std::endl;
                exit(9);
            }

            bool success_a = false;
            bool success_b = false;
            if (query_tail_left.length() > 1) {
                auto [lscore, lsuccess, lcigar] = Optimal(query_tail_left, target_tail_left, max_score);
                success_a = lsuccess;
                total_tail_alignments++;
                total_tail_length += lcigar.length();
                lsco = lscore;
                lcig = lcigar;
            } else {
                success_a = true;
                if (query_tail_left.length() == 1)  {
                    lcig = 'X';
                    lsco = -4;
                }
            }

            // If tail alignment leads to read dropping under ani threshold return failed alignment
            if (!success_a) {
                return { false, 0, "" };
            }
            max_score += lsco;

            if (query_tail_right.length() > 1) {
                auto [rscore, rsuccess, rcigar] = Optimal(query_tail_right, target_tail_right, max_score);
                success_b = rsuccess;
                total_tail_alignments++;
                total_tail_length += rcigar.length();
                dummy += rscore;
                rcig = rcigar;
                rsco = rscore;
            } else {
                success_b = true;
                if (query_tail_right.length() == 1)  {
                    rcig = 'X';
                    rsco = -4;
                }
            }

            // If tail alignment leads to read dropping under ani threshold return failed alignment
            if (!success_b) {
                return { false, 0, "" };
            }

            std::string total_cigar = lcig + mid_cigar + rcig;


//            auto [rs, gs, len] = GetAlignmentPositions(anchor, gene.length(), query.length());
//
//            m_aligner.Alignment(
//                    std::string_view(query.c_str() + rs, len),
//                    std::string_view(gene.c_str() + gs, len));
//            auto ref_cig = m_aligner.GetAlignmentCigar();
//            if (ref_cig != total_cigar && !((ref_cig.find('D') < ref_cig.length()) && ((ref_cig.find('I') < ref_cig.length())))) {
//                std::cout << anchor.ToString() << std::endl;
//                ExtendSeed(second, query, gene);
//                std::cout << anchor.ToString() << std::endl;
//                std::cout << total_cigar << std::endl;
//                std::cout << lcig << " " << mid_cigar << " " << rcig << std::endl;
//                m_aligner.PrintAlignment();
//                Utils::Input();
//            }

            return { true, mid_score + lsco + rsco, total_cigar };
        }

        static std::tuple<size_t, size_t, size_t> GetAlignmentPositions(CAlignmentAnchor& anchor, size_t gene_length, size_t read_length) {
            int abs_pos = anchor.Front().genepos - anchor.Front().readpos;
            size_t read_start = abs_pos < 0 ? -abs_pos : 0;
            size_t gene_start = std::max(abs_pos, 0);
            size_t overlap = std::min(gene_length - gene_start, read_length - read_start);
            return { read_start, gene_start, overlap };
        }

        inline std::pair<int, std::string> Approximate(std::string &query, std::string &gene) {
//            auto [ascore, acigar] = FastAligner::FastAlign(query, gene);
//            cigar_ani = WFA2Wrapper::CigarANI(acigar);
//            cigar = acigar;
//            score = ascore;


            return { 0, "" };
        }

        inline std::tuple<int, bool, std::string> Optimal(std::string_view query, std::string_view gene, int max_score) {
            m_aligner.Alignment(query, gene, max_score);
            return { m_aligner.GetAligner().getAlignmentScore(), m_aligner.Success(), m_aligner.GetAligner().getAlignmentCigar() };
        }

        void operator() (AlignmentAnchorList& anchors, AlignmentResultList& results, std::string& sequence) {

            constexpr bool alignment_verbose = false;

            auto& fwd = sequence;
            auto rev = KmerUtils::ReverseComplement(fwd);

            bool reversed = false;

            size_t read_len = sequence.length();


            std::sort(anchors.begin(), anchors.end(), [](CAlignmentAnchor const& a, CAlignmentAnchor const& b) {
                return a.total_length > b.total_length;
            });

            int take_top = m_align_top;



            int last_score = 0;


            bm_seedext.Start();
            ExtendAllAnchors(anchors, fwd, rev);
            bm_seedext.Stop();

            std::sort(anchors.begin(), anchors.end(), [](CAlignmentAnchor& a, CAlignmentAnchor& b) {
               return a.UpdateLength() > b.UpdateLength();
            });

            for (auto& anchor : anchors) {
//                if (anchor.chain.size() == 1) continue;

                auto& read = anchor.forward ? fwd : rev;

                // Is indel between anchor seeds?
                int anchor_indels = SeedIndel(anchor);

                // Get Resources
                auto& genome = m_genome_loader.GetGenome(anchor.taxid);
                auto& gene = genome.GetGeneOMP(anchor.geneid);

                // Absolute read positioning with respect to gene
                int abs_pos = anchor.Front().genepos - anchor.Front().readpos;


                size_t read_start = abs_pos < 0 ? -abs_pos : 0;
                size_t gene_start = std::max(abs_pos, 0);
                size_t overlap = std::min(gene.Sequence().length() - gene_start, read_len - read_start);

                std::string query = std::string(read.c_str() + read_start, overlap);
                std::string reference = gene.Sequence().substr(gene_start, overlap);


                double cigar_ani = 0;
                std::string cigar = "";
                int score;

                // anchor_indels == 0 means no indels likely between seeds of anchor
                bool approximate_alignment = m_fastalign && anchor_indels == 0;

                bm_alignment.Start();
                if (approximate_alignment) {
                    auto [successt, scoret, cigart] = Align(anchor, read, gene.Sequence());
                    if (!successt)  {
                        bm_alignment.Stop();
                        goto alignment_end;
                    }

                    cigar = cigart;
                    score = scoret;
                    cigar_ani = WFA2Wrapper::CigarANI(cigart);
                } else {
                    Benchmark bm_local{"alignment"};
                    m_aligner.Alignment(query, reference, MaxScore(m_max_score_ani, query.length()));
                    bm_local.Stop();

                    cigar_ani = WFA2Wrapper::CigarANI(m_aligner.GetAligner().getAlignmentCigar());

                    cigar = m_aligner.GetAligner().getAlignmentCigar();
                    score = m_aligner.GetAligner().getAlignmentScore();

//                    m_aligner.PrintAlignment();

                    if constexpr(alignment_verbose) {
#pragma omp critical (alignment_debug)
                        {
                            std::cerr << "Alignment\t" << bm_local.GetDuration(Time::microseconds);
                            std::cerr << '\t' << query.length() << '\t' << reference.length() << '\t';
                            std::cerr << cigar_ani << '\t' << anchor.total_length << '\n';
                        }
                    }
                }
                bm_alignment.Stop();

                if (cigar_ani >= m_max_score_ani) {
                    m_alignment_result.Set(anchor.taxid, anchor.geneid, abs_pos, !reversed);
                    m_alignment_result.Set(score, cigar);
                    results.emplace_back(m_alignment_result);
                    total_alignments++;
                }
                alignment_end:
                if (--take_top <= 0) {
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