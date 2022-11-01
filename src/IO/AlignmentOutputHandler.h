//
// Created by fritsche on 01/09/22.
//

#pragma once

#include "WFA2Wrapper.h"
#include <cstdint>
#include <string>
#include <vector>
#include "Constants.h"
#include "BufferedOutput.h"
#include "VarkitInterface.h"
#include "KmerLookup.h"
#include "FastxReader.h"
#include "SamHandler.h"
#include "SNP.h"
#include <htslib/sam.h>

namespace protal {
    class AlignmentResult {
        static constexpr int DEFAULT_ALIGNMENT_SCORE = INT32_MIN;
        static constexpr int DEFAULT_VALUE = UINT32_MAX;

        /* Only stores crucial alignment information, no information about header sequence etc
         * This is to save space and it is not necessary as this information can be passed down by another function.
         */
        int m_alignment_score = DEFAULT_ALIGNMENT_SCORE;

        uint32_t m_taxid = DEFAULT_VALUE;
        uint32_t m_geneid = DEFAULT_VALUE;
        uint32_t m_genepos = DEFAULT_VALUE;

        bool m_forward = false;
        std::string m_cigar;
        std::string m_compressed_cigar;



    public:
        AlignmentResult(int alignment_score, std::string& cigar, uint32_t taxid, uint32_t geneid, uint32_t genepos, bool forward) :
                m_alignment_score(alignment_score),
                m_cigar(cigar),
                m_taxid(taxid),
                m_geneid(geneid),
                m_genepos(genepos),
                m_forward(forward) {};

        AlignmentResult(int alignment_score, std::string&& cigar, uint32_t taxid, uint32_t geneid, uint32_t genepos, bool forward) :
                m_alignment_score(alignment_score),
                m_cigar(std::move(cigar)),
                m_taxid(taxid),
                m_geneid(geneid),
                m_genepos(genepos),
                m_forward(forward) {};

        AlignmentResult() {};

        std::string& Cigar() {
            return m_cigar;
        }

        std::string Cigar() const {
            return m_cigar;
        }

        uint32_t Taxid() const {
            return m_taxid;
        }

        uint32_t GeneId() const {
            return m_geneid;
        }

        uint32_t GenePos() const {
            return m_genepos;
        }

        bool Forward() const {
            return m_forward;
        }

        int AlignmentScore() const {
            return m_alignment_score;
        }


        void Set(int score, std::string& cigar) {
            m_cigar = cigar;
            m_alignment_score = score;
        }

        void Set(int score, std::string&& cigar) {
            m_cigar = cigar;
            m_alignment_score = score;
        }

        void Set(size_t taxid, size_t geneid, size_t abs_pos, bool forward) {
            m_taxid = taxid;
            m_geneid = geneid;
            m_genepos = abs_pos;
            m_forward = forward;
        }


        void Reset() {
            m_alignment_score = DEFAULT_ALIGNMENT_SCORE;
            m_cigar.clear();
        }

        bool IsSet() const {
            return m_alignment_score != DEFAULT_ALIGNMENT_SCORE;
        }

        std::string ToString() {
            std::string str = "";
            str += std::to_string(AlignmentScore()) + '\t';
            str += Cigar() + '\t';
            str += std::to_string(Taxid()) + '\t';
            str += std::to_string(GeneId()) + '\t';
            str += std::to_string(GenePos()) + '\t';
            str += std::to_string(Forward());
            return str;
        }
    };

    static int ScorePairedAlignment(AlignmentResult const& a, AlignmentResult const& b) {
        int score = 0;
        if (a.IsSet() && b.IsSet()) {
            return (a.AlignmentScore() + b.AlignmentScore()) / 2.2;
        } else if (a.IsSet()) {
            return a.AlignmentScore() - 2;
        } else if (b.IsSet()) {
            return b.AlignmentScore() - 2;
        }
        return score;
    }


    static int ScorePairedAlignment(PairedAlignment const& pair) {
        return ScorePairedAlignment(pair.first, pair.second);
    }

    static void ArtoSAM(SamEntry &sam, AlignmentResult const& ar, AlignmentInfo &info, FastxRecord &record) {
        sam.m_qname = record.id;
        sam.m_flag = 0;
        sam.m_rname = std::to_string(ar.Taxid()) + "_" + std::to_string(ar.GeneId());
        sam.m_pos = ar.GenePos() + info.gene_start_offset;
        sam.m_mapq = ar.AlignmentScore();
        sam.m_cigar = info.compressed_cigar;
        sam.m_rnext = "*";
        sam.m_pnext = 0;
        sam.m_tlen = ar.Cigar().length();
        sam.m_seq = ar.Forward() ? record.sequence : KmerUtils::ReverseComplement(record.sequence);
        sam.m_qual = record.quality;
    }

    using SNPList = std::vector<SNP>;
    /*
     * Protal Output Handler
     */
    class ProtalOutputHandler {
    private:
        std::ostream& m_sam_os;
        BufferedStringOutput m_sam_output;
        SamEntry m_sam;

        AlignmentInfo m_info;

        double m_min_cigar_ani = 0;
    public:
        size_t alignments = 0;

        ProtalOutputHandler(std::ostream& sam_os, size_t varkit_buffer_capacity, size_t sam_buffer_capacity, double min_cigar_ani=0.0f) :
                m_sam_os(sam_os),
                m_sam_output(sam_buffer_capacity),
                m_min_cigar_ani(min_cigar_ani) {}

        ProtalOutputHandler(ProtalOutputHandler const& other) :
                m_sam_os(other.m_sam_os),
                m_sam_output(other.m_sam_output.Capacity()),
                m_min_cigar_ani(other.m_min_cigar_ani) {}

        ~ProtalOutputHandler() {
#pragma omp critical(sam_output)
            m_sam_output.Write(m_sam_os);
        }


        void operator () (AlignmentResultList& alignment_results, FastxRecord& record, size_t read_id=0, bool first_pair=true) {
            if (alignment_results.empty()) return;
//            bool multiple_best = alignment_results.size() > 1 && alignment_results[0].AlignmentScore() == alignment_results[1].AlignmentScore();
//
//            if (multiple_best) {
//                return;
//            }

            auto& best = alignment_results.front();

            if (WFA2Wrapper::CigarANI(best.Cigar()) < m_min_cigar_ani) {
                return;
            }

            /* ##############################################################
             * Output alignments.
             */
            bool first = true;
            for (auto& ar :  alignment_results) {
                WFA2Wrapper::GetAlignmentInfo(m_info, ar.Cigar());
                ArtoSAM(m_sam, ar, m_info, record);
                Flag::SetPairedEnd(m_sam.m_flag, false, false, first_pair, !first_pair);
                Flag::SetRead2Unmapped(m_sam.m_flag, true);
                Flag::SetRead1ReverseComplement(m_sam.m_flag, ar.Forward());
                Flag::SetNotPrimaryAlignment(m_sam.m_flag, !first);

                alignments++;
                if (!m_sam_output.Write(m_sam.ToString())) {
#pragma omp critical(sam_output)
                    m_sam_output.Write(m_sam_os);
                }
                first = false;
            }
        }
    };


    /*
     * Protal Output Handler
     */
    class ProtalPairedOutputHandler {
    private:
        std::ostream& m_sam_os;
        BufferedStringOutput m_sam_output;
        SamEntry m_sam1;
        SamEntry m_sam2;

        AlignmentInfo m_info;

        double m_min_cigar_ani = 0;
    public:
        size_t alignments = 0;

        ProtalPairedOutputHandler(std::ostream& sam_os, size_t varkit_buffer_capacity, size_t sam_buffer_capacity, double min_cigar_ani=0.0f) :
                m_sam_os(sam_os),
                m_sam_output(sam_buffer_capacity),
                m_min_cigar_ani(min_cigar_ani) {}

        ProtalPairedOutputHandler(ProtalPairedOutputHandler const& other) :
                m_sam_os(other.m_sam_os),
                m_sam_output(other.m_sam_output.Capacity()),
                m_min_cigar_ani(other.m_min_cigar_ani) {}

        ~ProtalPairedOutputHandler() {
#pragma omp critical(sam_output)
            m_sam_output.Write(m_sam_os);
        }


        void operator () (PairedAlignmentResultList& alignment_results, FastxRecord& record1, FastxRecord& record2, size_t read_id=0, bool first_pair=true) {
            if (alignment_results.empty()) return;
//            bool multiple_best = alignment_results.size() > 1 && alignment_results[0].AlignmentScore() == alignment_results[1].AlignmentScore();
//
//            if (multiple_best) {
//                return;
//            }

            auto& best = alignment_results.front();

//            size_t anum = 0;
//            std::cout << record1.id << std::endl;
//            for (auto& [ar1, ar2] :  alignment_results) {
//                std::cout << anum++ << " SCORE: " << ScorePairedAlignment(ar1, ar2) << '\t' << ar1.AlignmentScore() << "{"<< ar1.Taxid() << "," << ar1.GeneId() << "}" << " ";
//                std::cout << ar2.AlignmentScore() << "{"<< ar2.Taxid() << "," << ar2.GeneId() << "}" << std::endl;
//            }
//            Utils::Input();

            if (WFA2Wrapper::CigarANI(best.first.Cigar()) < m_min_cigar_ani) {
                return;
            }



            /* ##############################################################
             * Output alignments.
             */
            bool first = true;
            for (auto& [ar1, ar2] :  alignment_results) {
                bool both = ar1.IsSet() && ar2.IsSet();

                if (ar1.IsSet()) {

                    WFA2Wrapper::GetAlignmentInfo(m_info, ar1.Cigar());
                    ArtoSAM(m_sam1, ar1, m_info, record1);
                    Flag::SetPairedEnd(m_sam1.m_flag, true, both, true);
                    Flag::SetRead2Unmapped(m_sam1.m_flag, !ar2.IsSet());
                    Flag::SetRead1ReverseComplement(m_sam1.m_flag, ar1.Forward());
                    Flag::SetNotPrimaryAlignment(m_sam1.m_flag, !first);
                }
                if (ar2.IsSet()) {
                    auto len = std::count_if(ar2.Cigar().begin(), ar2.Cigar().end(), [](char c) {
                        return c != 'I';
                    });

                    if (len != record2.sequence.length()) {
                        std::cout << len << std::endl;
                        std::cout << record2.sequence << std::endl;
                        std::cout << ar2.Cigar() << std::endl;
                        exit(11);
                    }
                    WFA2Wrapper::GetAlignmentInfo(m_info, ar2.Cigar());
                    ArtoSAM(m_sam2, ar2, m_info, record2);
                    Flag::SetPairedEnd(m_sam2.m_flag, true, both, false, true);
                    Flag::SetRead1Unmapped(m_sam2.m_flag, !ar1.IsSet());
                    Flag::SetRead2ReverseComplement(m_sam2.m_flag, ar2.Forward());
                    Flag::SetNotPrimaryAlignment(m_sam2.m_flag, !first);
                }
                if (both) {
                    m_sam1.m_rnext = m_sam2.m_rname;
                    m_sam2.m_rnext = m_sam1.m_rname;
                    m_sam1.m_pnext = m_sam2.m_pos;
                    m_sam2.m_pnext = m_sam1.m_pos;
                }


                first = false;
                alignments++;
                if (ar1.IsSet() && !m_sam_output.Write(m_sam1.ToString())) {
#pragma omp critical(sam_output)
                    m_sam_output.Write(m_sam_os);
                }
                if (ar2.IsSet() && !m_sam_output.Write(m_sam2.ToString())) {
#pragma omp critical(sam_output)
                    m_sam_output.Write(m_sam_os);
                }
            }
        }
    };
}