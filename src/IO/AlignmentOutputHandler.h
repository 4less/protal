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

    /*
     * Varkit Output Handler
     *
     */
    class VarkitOutputHandler {
    private:
        std::ostream& m_os;
        std::ostream& m_sam_os;

        BufferedStringOutput m_sam_output;
        BufferedOutput<ClassificationLine> m_output;
        ClassificationLine m_line;
        SamEntry m_sam;
        AlignmentInfo m_info;

    public:
        size_t alignments = 0;

        VarkitOutputHandler(std::ostream& os, std::ostream& sam_os, size_t varkit_buffer_capacity, size_t sam_buffer_capacity) :
                m_os(os),
                m_sam_os(sam_os),
                m_output(varkit_buffer_capacity),
                m_sam_output(sam_buffer_capacity) {}

        VarkitOutputHandler(VarkitOutputHandler const& other) :
                m_os(other.m_os),
                m_sam_os(other.m_sam_os),
                m_output(other.m_output.Capacity()),
                m_sam_output(other.m_sam_output.Capacity()) {}

        ~VarkitOutputHandler() {
#pragma omp critical(varkit_output)
            m_output.Write(m_os);
#pragma omp critical(sam_output)
            m_sam_output.Write(m_sam_os);
        }


        void operator () (AlignmentResultList& alignment_results, FastxRecord& record) {
            if (alignment_results.empty()) return;

            int best_score = alignment_results[0].AlignmentScore();
            size_t best_score_idx = 0;
            size_t best_score_occurences = 0;

            for (auto i = 1; i < alignment_results.size(); i++) {
                auto& result = alignment_results[i];
                if (result.AlignmentScore() == best_score) {
                    best_score_occurences++;
                }
                if (result.AlignmentScore() < best_score) {
                    best_score_occurences = 1;
                    best_score = result.AlignmentScore();
                    best_score_idx = i;
                }
            }

            if (best_score_occurences > 1) {
                return;
            }

            m_line.Reset();

            auto& best = alignment_results[best_score_idx];

            WFA2Wrapper::GetAlignmentInfo(m_info, best.Cigar());

            ToClassificationLine(m_line, best.Taxid(), best.GeneId(), best.GenePos() + m_info.alignment_start, m_info.alignment_length, 1, m_info.alignment_ani);

            m_sam.m_qname = record.id;
            m_sam.m_flag = 0;
            m_sam.m_rname = std::to_string(best.Taxid()) + "_" + std::to_string(best.GeneId());
            m_sam.m_pos = best.GenePos() + m_info.alignment_start;
            m_sam.m_mapq = best.AlignmentScore();
            m_sam.m_cigar = best.Cigar();
            m_sam.m_rnext = "*";
            m_sam.m_pnext = 0;
            m_sam.m_tlen = best.Cigar().length();
            m_sam.m_seq = record.sequence;
            m_sam.m_qual = record.quality;

//            std::cout << m_sam.ToString() << std::endl;

            alignments++;
            if (!m_output.Write(m_line)) {
#pragma omp critical(varkit_output)
                m_output.Write(m_os);
            }
            if (!m_sam_output.Write(m_sam.ToString())) {
#pragma omp critical(sam_output)
                m_sam_output.Write(m_sam_os);
            }
        }
    };
}