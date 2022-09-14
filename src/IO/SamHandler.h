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
