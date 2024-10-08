//
// Created by fritsche on 04/11/22.
//

#pragma once

#include "Constants.h"
#include "SamHandler.h"
//#include "WFA2Wrapper.h"
//#include "WFA2Wrapper2.h"
#include <cstdio>
#include <cstddef>
#include "AlignmentUtils.h"

namespace protal {
    static const std::pair<TaxId, GeneId> ExtractTaxidGeneid(std::string &ref) {
        int delim_pos = -1;
        while (ref[++delim_pos] != '_');
        int next_delim_pos = delim_pos;
        while (next_delim_pos < ref.size() && ref[++next_delim_pos] != '_');
        return { stoul(ref.substr(0, delim_pos)), stoul(ref.substr(delim_pos+1, next_delim_pos - delim_pos - 1)) };
    }

    struct InternalReadAlignment {
        size_t readid;
        size_t pairid;              //128
        TaxId taxid;
        GeneId geneid;
        GenePos genepos;
        int alignment_score;        //256
        double alignment_ani;
        uint32_t alignment_length;
        bool forward;
        bool should_be_paired = false;
        bool is_paired = false;
        bool benchmark = false;
        bool correct = false;


        std::string ToString() const {
            std::string str;
            str += std::to_string(readid) + '\t';
            str += std::to_string(pairid) + '\t';
            str += std::to_string(taxid) + '\t';
            str += std::to_string(geneid) + '\t';
            str += std::to_string(genepos) + '\t';
            str += std::to_string(alignment_score) + '\t';
            str += std::to_string(alignment_ani) + '\t';
            str += std::to_string(alignment_length) + '\t';
            str += std::to_string(forward);
            return str;
        }

        InternalReadAlignment(
                size_t readid,
                size_t pairid,
                TaxId taxid,
                GeneId geneid,
                GenePos genepos,
                bool forward,
                int alignment_score,
                double alignment_ani,
                uint32_t alignment_length) :
                readid(readid),
                taxid(taxid),
                geneid(geneid),
                genepos(genepos),
                forward(forward),
                pairid(pairid),
                alignment_score(alignment_score),
                alignment_ani(alignment_ani),
                alignment_length(alignment_length) {};

        InternalReadAlignment(size_t readid, SamEntry &sam, bool benchmark = false) :
                readid(readid),
                taxid(0),
                geneid(0),
                genepos(sam.m_pos),
                forward(false),
                benchmark(benchmark),
                pairid(Flag::IsRead2(sam.m_flag))
        {
            auto &[ tid, gid ] = ExtractTaxidGeneid(sam.m_rname);
            taxid = tid;
            geneid = gid;

            CigarInfo info;
            CompressedCigarInfo(sam.m_cigar, info);

            if (info.clipped_alignment_length == 0 || info.unclipped_alignment_length == 0) {
                std::cout << info.ToString() << std::endl;
                exit(1);
            }
            alignment_score = info.Score();
            alignment_ani = info.Ani();
            alignment_length = info.clipped_alignment_length;

            if (benchmark) {
                correct = sam.m_qname.substr(0, sam.m_rname.length()) == sam.m_rname;
            }
        };


        InternalReadAlignment(InternalReadAlignment const &other) :
                readid(other.readid),
                taxid(other.taxid),
                geneid(other.geneid),
                genepos(other.genepos),
                forward(other.forward),
                pairid(other.pairid),
                alignment_score(other.alignment_score),
                alignment_ani(other.alignment_ani),
                alignment_length(other.alignment_length) {};
    };

    enum PairedStatus {
        PAIRED,
        FIRST_ONLY,
        SECOND_ONLY,
        NONE
    };

    static PairedStatus PairedStatus(OptIRAPair const& pair) {
        if (pair.first.has_value() && pair.second.has_value())
            return PairedStatus::PAIRED;
        if (pair.first.has_value())
            return PairedStatus::FIRST_ONLY;
        if (pair.second.has_value())
            return PairedStatus::SECOND_ONLY;
        return PairedStatus::NONE;
    }

    static bool IsAmbiguous(InternalNonUniqueReadOptIRAPairList const& list) {
        if (list.size() <= 1) return false;

        auto& pair1 = list[0];
        auto& pair2 = list[1];
        auto status1 = PairedStatus(pair1);
        auto status2 = PairedStatus(pair2);

        bool both_paired = status1 == PAIRED && status2 == PAIRED;
        bool both_single =
                status1 == FIRST_ONLY && status2 == FIRST_ONLY ||
                status1 == FIRST_ONLY && status2 == SECOND_ONLY ||
                status1 == SECOND_ONLY && status2 == FIRST_ONLY ||
                status1 == SECOND_ONLY && status2 == SECOND_ONLY;
        bool one_paired_one_single =
                status1 == PAIRED && status2 == FIRST_ONLY ||
                status1 == PAIRED && status2 == SECOND_ONLY ||
                status1 == FIRST_ONLY && status2 == PAIRED ||
                status1 == SECOND_ONLY && status2 == PAIRED;

        if (both_paired) {
            return  pair1.first.value().alignment_score == pair2.first.value().alignment_score &&
                    pair1.second.value().alignment_score == pair2.second.value().alignment_score;
        }
        if (both_single) {
            auto& first = pair1.first.has_value() ? pair1.first.value() : pair1.second.value();
            auto& second = pair2.first.has_value() ? pair2.first.value() : pair2.second.value();
            return first.alignment_score == second.alignment_score;
        }

        return false;
    }
}