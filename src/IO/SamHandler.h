//
// Created by fritsche on 26/08/22.
//

#ifndef PROTAL_SAMHANDLER_H
#define PROTAL_SAMHANDLER_H


#include <string>

namespace protal {
    using QNAME_t = std::string;
    using FLAG_t = uint16_t;
    using RNAME_t = std::string;
    using POS_t = uint32_t;
    using MAPQ_t = uint8_t;
    using CIGAR_t = std::string;
    using TLEN_t = int64_t;
    using SEQ_t = std::string;
    using QUAL_t = std::string;

    class Flag {
    public:
        static void SetPairedEnd(FLAG_t &flag, bool ispaired, bool bothalign, bool is_read1=false, bool is_read2=false) {
            flag |= ispaired;
            flag |= (is_read1 << 6);
            flag |= (is_read2 << 7);
        }

        static void SetPairBothAlign(FLAG_t &flag, bool both_align) {
            flag |= (both_align << 1);
        }

        static void SetRead1Unmapped(FLAG_t &flag, bool unmapped) {
            flag |= (unmapped << 2);
        }

        static void SetRead2Unmapped(FLAG_t &flag, bool unmapped) {
            flag |= (unmapped << 3);
        }

        static void SetRead1ReverseComplement(FLAG_t &flag, bool rc) {
            flag |= (rc << 4);
        }

        static void SetRead2ReverseComplement(FLAG_t &flag, bool rc) {
            flag |= (rc << 5);
        }

        static void SetNotPrimaryAlignment(FLAG_t &flag, bool not_primary) {
            flag |= (not_primary << 8);
        }

        static void SetAlignmentFailsQuality(FLAG_t &flag, bool fails) {
            flag |= (fails << 9);
        }

        static void SetDuplicate(FLAG_t &flag, bool duplicate) {
            flag |= (duplicate << 10);
        }

        static void SetSupplementaryAlignment(FLAG_t &flag, bool is_supplementary) {
            flag |= (is_supplementary << 11);
        }


        static bool IsPaired(FLAG_t flag) {
            return flag & 1;
        }

        static bool IsPairBothALign(FLAG_t flag) {
            return flag & (1 << 1);
        }

        static bool IsRead1Unmapped(FLAG_t flag) {
            return flag & (1 << 2);
        }

        static bool IsRead2Unmapped(FLAG_t flag) {
            return flag & (1 << 3);
        }

        static bool IsRead1ReverseComplement(FLAG_t flag) {
            return flag & (1 << 4);
        }

        static bool IsRead2ReverseComplement(FLAG_t flag) {
            return flag & (1 << 5);
        }

        static bool IsRead1(FLAG_t flag) {
            return flag & (1 << 6);
        }

        static bool IsRead2(FLAG_t flag) {
            return flag & (1 << 7);
        }

        static bool IsNotPrimaryAlignment(FLAG_t flag) {
            return flag & (1 << 8);
        }

        static bool IsAlignmentFailsQuality(FLAG_t flag) {
            return flag & (1 << 9);
        }

        static bool IsDuplicate(FLAG_t &flag) {
            return flag & (1 << 10);
        }

        static bool IsSupplementaryAlignment(FLAG_t &flag) {
            return flag & (1 << 11);
        }
    };

    struct SamEntry {
        QNAME_t m_qname;
        FLAG_t m_flag;
        RNAME_t m_rname;
        POS_t m_pos;
        MAPQ_t m_mapq;
        CIGAR_t m_cigar;
        RNAME_t m_rnext;
        POS_t m_pnext;
        TLEN_t m_tlen;
        SEQ_t m_seq;
        QUAL_t m_qual;

        std::string ToString() {
            return  m_qname + '\t'
                    + std::to_string(m_flag) + '\t'
                    + m_rname + '\t'
                    + std::to_string(m_pos) + '\t'
                    + std::to_string(m_mapq) + '\t'
                    + m_cigar + '\t'
                    + m_rnext + '\t'
                    + std::to_string(m_pnext) + '\t'
                    + std::to_string(m_tlen) + '\t'
                    + m_seq + '\t'
                    + m_qual;
        }
    };

    class SamHandler {

    };
}


#endif //PROTAL_SAMHANDLER_H
