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
        size_t m_max_score_ani = 0;


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


        SimpleAlignmentHandler(GenomeLoader& genome_loader, WFA2Wrapper& aligner, size_t kmer_size, size_t align_top, size_t max_score_ani) :
                m_genome_loader(genome_loader),
                m_aligner(aligner),
                m_kmer_size(kmer_size),
                m_align_top(align_top),
                m_max_score_ani(max_score_ani) {}


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

        size_t MaxScore(size_t total_length, double target_ani) {
            return target_ani * total_length;
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


        inline std::pair<int, std::string> Align(AlignmentAnchor const& anchor, std::string const& query, std::string const& gene) {
            //                front_end_q                         back_start_q
            //                v                                   v
            //----------------SSSSSSSSS------------------SSSSSSSSS-------------------- Query
            //----------------SSSSSSSSS------------------SSSSSSSSS-------------------- Ref
            //                ^                                   ^
            //                front_end_r                         back_start_r

            constexpr bool debug = true;
            // Phase one = Extend left seed


            auto [lefta, righta] = ExtendSeed(anchor.a, m_kmer_size, query, gene);
            auto [leftb, rightb] = ExtendSeed(anchor.b, m_kmer_size, query, gene);



            size_t qendl = anchor.a.readpos - lefta;
            size_t rendl = anchor.a.genepos - lefta;
//            std::cout << qendl << " " << rendl << std::endl;
            size_t qstartl = qendl <= rendl ? 0 : qendl - rendl;
            size_t rstartl = qendl >= rendl ? 0 : rendl - qendl;

//            std::cout << "align query: " << qstartl << " - " << qendl << std::endl;
//            std::cout << "align ref:   " << rstartl << " - " << rendl << std::endl;

            size_t qstartr = anchor.b.readpos + m_kmer_size + rightb;
            size_t rstartr = anchor.b.genepos + m_kmer_size + rightb;
            size_t qtail = query.length() - qstartr;
            size_t rtail = gene.length() - rstartr;

            size_t qendr = qtail <= rtail ? query.length() : qstartr + rtail;
            size_t rendr = qtail >= rtail ? gene.length() : rstartr + qtail;

            auto q_ltail = std::string_view(query.c_str() + qstartl, qendl - qstartl);
            auto r_ltail = std::string_view(gene.c_str() + rstartl, rendl - rstartl);
            auto q_rtail = std::string_view(query.c_str() + qstartr, qendr - qstartr);
            auto r_rtail = std::string_view(gene.c_str() + rstartr, rendr - rstartr);
            std::string lcig = "";
            std::string rcig = "";
            double max_ani = 0.9;

            bool success_a = false;
            bool success_b = false;
            if (qendl - qstartl > 1) {
                Benchmark bm_local("");
                auto [lscore, lsuccess, lcigar] = Optimal(q_ltail, r_ltail, max_ani);
                success_a = lsuccess;
                bm_local.Stop();
                total_tail_alignments++;
                total_tail_length += lcigar.length();
                dummy += lscore;
                lcig = lcigar;
                auto qlen = qendl - qstartl;
                auto rlen = rendl - rstartl;

//#pragma omp critical (alignment_debug)
//                    {
//                        std::cerr << "Alignment\t" << bm_local.GetDuration(Time::microseconds);
//                        std::cerr << '\t' << qlen << '\t' << rlen << '\t';
//                        std::cerr << WFA2Wrapper::CigarANI(lcigar) << '\t' << anchor.hit_anchor_count << '\n';
//                    }
            } else {
                success_a = true;
            }


            if (qendr - qstartr > 1) {
                Benchmark bm_local("");
                auto [rscore, rsuccess, rcigar] = Optimal(q_rtail, r_rtail, max_ani);
                success_b = rsuccess;
                bm_local.Stop();
                total_tail_alignments++;
                total_tail_length += rcigar.length();
                dummy += rscore;
                rcig = rcigar;

                auto qlen = qendl - qstartl;
                auto rlen = rendl - rstartl;
//
//#pragma omp critical (alignment_debug)
//                    {
//                        std::cerr << "Alignment\t" << bm_local.GetDuration(Time::microseconds);
//                        std::cerr << '\t' << qlen << '\t' << rlen << '\t';
//                        std::cerr << WFA2Wrapper::CigarANI(rcigar) << '\t' << anchor.hit_anchor_count << '\n';
//                    }

            } else {
                success_b = true;
            }

            bool success = success_a && success_b;

            dummy += qendr;
            dummy += rendr;

            if constexpr(debug) {
                if (!success) {
                    std::cout << "\n\n--------------------------------------------" << std::endl;
                    auto [rs, gs, len] = GetAlignmentPositions(anchor, gene.length(), query.length());
                    auto query_view = std::string_view(query.c_str() + rs, len);
                    auto gene_view = std::string_view(gene.c_str() + gs, len);
//                    auto query_view = std::string(query.c_str() + rs, len);
//                    auto gene_view = std::string(gene.c_str() + gs, len);
                    m_aligner.Alignment1(
                            query_view,
                            gene_view);
                    double ani = WFA2Wrapper::CigarANI(m_aligner.GetAlignmentCigar());

                    std::cout << anchor.a.ToString() << std::endl;
                    std::cout << "<- " << lefta << " " << righta << " ->" << std::endl;
                    std::cout << anchor.b.ToString() << std::endl;
                    std::cout << "<- " << leftb << " " << rightb << " ->" << std::endl;
                    std::cout << "QUERY: " << q_ltail << " " << q_rtail << std::endl;
                    std::cout << "Cigar: " << lcig << " " << rcig << std::endl;
                    std::cout << "Ref:   " << r_ltail << " " << r_rtail << std::endl;
                    if (ani > 0.9) {
                        if (m_aligner.Success())
                            m_aligner.PrintAlignment();
                        else {
                            std::cout << "no success" << std::endl;
                            std::cout << std::string(query_view.data(), query_view.length()) << std::endl;
                            std::cout << query_view << std::endl;
                            std::cout << std::string(gene_view.data(), gene_view.length()) << std::endl;
                            std::cout << gene_view << std::endl;
                        }
                        Utils::Input();
                    }
                }
            }
            return { 0, "" };
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

        inline std::tuple<int, bool, std::string> Optimal(std::string_view query, std::string_view gene, double min_ani=0.9) {
            m_aligner.Alignment(query, gene, min_ani);
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


                total_alignments++;

                double cigar_ani = 0;
                std::string cigar = "";
                int score = -1000;

                // anchor_indels == 0 means no indels likely between seeds of anchor
                bool approximate_alignment = anchor_indels == 0;
//                bool approximate_alignment = false;

//                auto q = std::string(read, read_len);
//                auto [lefta, righta] = ExtendSeed(anchor.a, 15, q, gene.Sequence());
//                auto [leftb, rightb] = ExtendSeed(anchor.b, 15, q, gene.Sequence());
//                dummy += (lefta + m_kmer_size + righta) + (leftb + m_kmer_size + rightb);

                bm_alignment.Start();
                Align(anchor, read, gene.Sequence());
                if (approximate_alignment) {
                    auto [ascore, acigar] = FastAligner::FastAlign(query, reference);
                    cigar_ani = WFA2Wrapper::CigarANI(acigar);
                    cigar = acigar;
                    score = ascore;
                } else {
                    Benchmark bm_local{"alignment"};
                    m_aligner.Alignment(query, reference);
                    bm_local.Stop();

                    cigar_ani = WFA2Wrapper::CigarANI(m_aligner.GetAligner().getAlignmentCigar());

                    cigar = m_aligner.GetAligner().getAlignmentCigar();
                    score = m_aligner.GetAligner().getAlignmentScore();


//                    std::cout << std::string(69, '-') << std::endl;
//                    std::cout << anchor.a.ToString() << " " << anchor.b.ToString() << std::endl;
//                    std::cout << std::string(29, '-') << "Extend first" << std::endl;
//                    std::cout << lefta  << " " << m_kmer_size << " " << righta <<  " -> " << (lefta + m_kmer_size + righta)  << std::endl;
//                    std::cout << std::string(29, '-') << "Extend second" << std::endl;
//                    std::cout << leftb  << " " << m_kmer_size << " " << rightb <<  " -> " << (leftb + m_kmer_size + rightb) << std::endl;
//                    std::cout << std::string(29, '-') << "Alignment" << std::endl;
//                    m_aligner.PrintAlignment();
//                    std::cout << anchor.hit_anchor_count << " -> " << (lefta + m_kmer_size + righta) + (leftb + m_kmer_size + rightb) << std::endl;

//                    Align(anchor, read, gene.Sequence());
//                    Utils::Input();

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