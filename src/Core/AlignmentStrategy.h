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

        size_t m_align_top = 3;
        double m_max_score_ani = 0;

        bool m_fastalign = false;


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
        size_t total_tail_alignments = 0;
        size_t total_tail_length = 0;
        Benchmark bm_alignment{ "Alignment" };
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

        inline int SeedIndel(AlignmentAnchor& anchor) {
            int rpos_diff = static_cast<int>(anchor.b.readpos) - static_cast<int>(anchor.a.readpos);
            int gpos_diff = static_cast<int>(anchor.b.genepos) - static_cast<int>(anchor.a.genepos);
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

        static int MaxScore(double min_ani, uint32_t overlap_length, int mismatch_penalty=4) {
            return (std::ceil((1 - min_ani) * static_cast<double>(overlap_length)) * mismatch_penalty) + 1;
        }


        inline std::tuple<bool, int, std::string> Align(AlignmentAnchor const& anchor, std::string const& query, std::string const& gene) {
            //                front_end_q                         back_start_q
            //                v                                   v
            //----------------SSSSSSSSS------------------SSSSSSSSS-------------------- Query
            //----------------SSSSSSSSS------------------SSSSSSSSS-------------------- Ref
            //                ^                                   ^
            //                front_end_r                         back_start_r

            constexpr bool debug = false;
            // Phase one = Extend left seed

            auto [lefta, righta] = ExtendSeed(anchor.a, m_kmer_size, query, gene);
            auto [leftb, rightb] = ExtendSeed(anchor.b, m_kmer_size, query, gene);

            size_t qstartm = anchor.a.readpos + m_kmer_size + righta;
            size_t qendm = anchor.b.readpos - leftb;
            size_t rstartm = anchor.a.genepos + m_kmer_size + righta;
            size_t rendm = anchor.b.genepos - leftb;


            size_t seed_a_size = lefta + m_kmer_size + righta;
            size_t seed_b_size = leftb + m_kmer_size + rightb;

            size_t qendl = anchor.a.readpos - lefta;
            size_t rendl = anchor.a.genepos - lefta;
            size_t qstartl = qendl <= rendl ? 0 : qendl - rendl;
            size_t rstartl = qendl >= rendl ? 0 : rendl - qendl;
            size_t ql_len = qendl - qstartl;
            size_t rl_len = rendl - rstartl;

            size_t qstartr = anchor.b.readpos + m_kmer_size + rightb;
            size_t rstartr = anchor.b.genepos + m_kmer_size + rightb;
            size_t qtail = query.length() - qstartr;
            size_t rtail = gene.length() - rstartr;

            size_t total_seed_cov = anchor.a.readpos + m_kmer_size + righta > anchor.b.readpos - leftb ?
                                    qstartr - qendl :
                                    seed_a_size+seed_b_size;

            size_t qendr = qtail <= rtail ? query.length() : qstartr + rtail;
            size_t rendr = qtail >= rtail ? gene.length() : rstartr + qtail;
            size_t qr_len = qendr - qstartr;
            size_t rr_len = rendr - rstartr;

            size_t overlap = qendr - qstartl;

            int max_score = MaxScore(m_max_score_ani, overlap, 4);

            std::string mid_cigar = "";
            int mid_score = 0;
            if (qendm > qstartm) {
                if constexpr(debug) {
                    std::cout << std::string_view(query.c_str() + qstartm, qendm - qstartm) << std::endl;
                    std::cout << std::string_view(gene.c_str() + rstartm, rendm - rstartm) << std::endl;
                }
                auto [score, cigar] = FastAligner::FastAlign(
                        std::string_view(query.c_str() + qstartm, qendm - qstartm),
                        std::string_view(gene.c_str() + rstartm, rendm - rstartm));
                mid_score = score;
                mid_cigar = std::string(seed_a_size, 'M') + cigar + std::string(seed_b_size, 'M');
            } else {
                mid_cigar = std::string(total_seed_cov, 'M');
            }

            if constexpr(debug) {
                std::cout << "Total seed coverage " << total_seed_cov << std::endl;
                std::cout << "mid_score " << mid_score << std::endl;
            }

            max_score += mid_score;

            if (max_score < 0) {
                return { false, -1, "" };
            }

            // create stringviews
            auto q_ltail = std::string_view(query.c_str() + qstartl, ql_len);
            auto r_ltail = std::string_view(gene.c_str() + rstartl, rl_len);
            auto q_rtail = std::string_view(query.c_str() + qstartr, qr_len);
            auto r_rtail = std::string_view(gene.c_str() + rstartr, rr_len);

            std::string lcig = "";
            std::string rcig = "";
            size_t lsco = 0;
            size_t rsco = 0;



            if (q_ltail.length() + q_rtail.length() + mid_cigar.length() != qendr - qstartl) {
                std::cout << "query length:     " << query.length() << std::endl;
                std::cout << "GeneLength:       " << gene.length() << std::endl;
                std::cout << "qendr - qstartl:  " << qendr - qstartl << std::endl;
                std::cout << "left_len:         " << ql_len << std::endl;
                std::cout << "right_len:        " << qr_len << std::endl;
                std::cout << "mid_len:          " << mid_cigar.length() << std::endl;
                std::cout << "Total seed cov:   " << total_seed_cov << std::endl;
                std::cout << "Anchor:           " << anchor.a.ToString() << " " << anchor.b.ToString() << std::endl;

                std::cout << q_ltail << " " << std::string_view(query.c_str() + qendl, qstartr -qendl) << " " << q_rtail << std::endl;
                std::cout << std::string(q_ltail.length() + 1, ' ') << mid_cigar << std::endl;
                exit(9);
            }

            bool success_a = false;
            bool success_b = false;
            if (qendl - qstartl > 1) {
                Benchmark bm_local("");
                auto [lscore, lsuccess, lcigar] = Optimal(q_ltail, r_ltail, max_score);
//                auto [lscore, lsuccess, lcigar] = Optimal(q_ltail, r_ltail, max_ani, q_ltail.length() + q_rtail.length(), mid_score);
                success_a = lsuccess;

                bm_local.Stop();
                total_tail_alignments++;
                total_tail_length += lcigar.length();
                lsco = lscore;
                lcig = lcigar;
            } else {
                success_a = true;
                if (qendl - qstartl)  {
                    lcig = 'X';
                    lsco = -4;
                }
            }



            // If tail alignment leads to read dropping under ani threshold return failed alignment
            if (!success_a) {
                return { false, 0, "" };
            }

            max_score += lsco;


            if (qendr - qstartr > 1) {
                Benchmark bm_local("");
                auto [rscore, rsuccess, rcigar] = Optimal(q_rtail, r_rtail, max_score);
//                auto [rscore, rsuccess, rcigar] = Optimal(q_rtail, r_rtail, max_ani, mid_cigar.length(), mid_score);
                success_b = rsuccess;

                bm_local.Stop();
                total_tail_alignments++;
                total_tail_length += rcigar.length();
                dummy += rscore;
                rcig = rcigar;
                rsco = rscore;
            } else {
                success_b = true;
                if (qendr - qstartr) {
                    rcig = 'X';
                    lsco = -4;
                }
            }

            // If tail alignment leads to read dropping under ani threshold return failed alignment
            if (!success_b) {
                return { false, 0, "" };
            }

            std::string total_cigar = lcig + mid_cigar + rcig;


            dummy += qendr;
            dummy += rendr;

            if constexpr(debug) {
                if (true) {
                    std::cout << "\n--------------------------------------------" << std::endl;
                    std::cout << q_ltail << " " << q_rtail << " " << qendl - qstartl << std::endl;
                    std::cout << r_ltail << " " << r_rtail << " " << qendr - qstartr << std::endl;
                    std::cout << "suca " << success_a << " " << success_b << " sub" << std::endl;
                    auto [rs, gs, len] = GetAlignmentPositions(anchor, gene.length(), query.length());
                    auto query_view = std::string_view(query.c_str() + rs, len);
                    auto gene_view = std::string_view(gene.c_str() + gs, len);
                    m_aligner.Alignment(
                            query_view,
                            gene_view, max_score);
                    double ani = WFA2Wrapper::CigarANI(m_aligner.GetAlignmentCigar());

                    std::cout << anchor.a.ToString() << std::endl;
                    std::cout << "<- " << lefta << " " << righta << " ->" << std::endl;
                    std::cout << anchor.b.ToString() << std::endl;
                    std::cout << "<- " << leftb << " " << rightb << " ->" << std::endl;
                    std::cout << "QUERY: " << q_ltail << " " << q_rtail << std::endl;
                    std::cout << "Cigar: " << lcig << " " << rcig << std::endl;
                    std::cout << "Ref:   " << r_ltail << " " << r_rtail << std::endl;

                    if (m_aligner.Success()) {
                        std::cout << "Ani: " << ani << std::endl;
                        std::cout << total_cigar << " " << WFA2Wrapper::CigarANI(total_cigar) << std::endl;
                        m_aligner.PrintAlignment();

                    } else {
                        std::cout << "no success" << std::endl;
                        std::cout << std::string(query_view.data(), query_view.length()) << std::endl;
                        std::cout << query_view << std::endl;
                        std::cout << std::string(gene_view.data(), gene_view.length()) << std::endl;
                        std::cout << gene_view << std::endl;
                        exit(8);
                    }
                    if (total_cigar != m_aligner.GetAlignmentCigar()) {
                        Utils::Input();
                    }
                }
            }

            return { true, mid_score + lsco + rsco, total_cigar };
        }

        static std::tuple<size_t, size_t, size_t> GetAlignmentPositions(AlignmentAnchor const& anchor, size_t gene_length, size_t read_length) {
            int abs_pos = anchor.a.genepos - anchor.a.readpos;
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

//            char* read = nullptr;
            bool reversed = false;

            size_t read_len = sequence.length();


            std::sort(anchors.begin(), anchors.end(), [](AlignmentAnchor const& a, AlignmentAnchor const& b) {
                return a.hit_anchor_count > b.hit_anchor_count;
            });

            int take_top = m_align_top;



            int last_score = 0;

            for (auto& anchor : anchors) {
                reversed = ReverseAnchor(anchor, read_len);
                auto& read = reversed ?
                       rev :
                       fwd;

                // Is indel between anchor seeds?
                int anchor_indels = SeedIndel(anchor);

                // Get Resources
                auto& genome = m_genome_loader.GetGenome(anchor.a.taxid);
                auto& gene = genome.GetGeneOMP(anchor.a.geneid);

                // Absolute read positioning with respect to gene
                int abs_pos = anchor.a.genepos - anchor.a.readpos;


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

                if (cigar_ani >= m_max_score_ani) {
                    m_alignment_result.Set(anchor.a.taxid, anchor.b.geneid, abs_pos, !reversed);
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