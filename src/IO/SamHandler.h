//
// Created by fritsche on 26/08/22.
//

#ifndef PROTAL_SAMHANDLER_H
#define PROTAL_SAMHANDLER_H


#include <string>

namespace protal {
    using QNAME_t = std::string;
    using FLAG_t = int;
    using RNAME_t = std::string;
    using POS_t = size_t;
    using MAPQ_t = uint8_t;
    using CIGAR_t = std::string;
    using ISIZE_t = size_t;
    using SEQ_t = std::string;
    using QUAL_t = std::string;

    struct SamEntry {
        QNAME_t m_qname;
        FLAG_t m_flag;
        RNAME_t m_rname;
        POS_t m_pos;
        MAPQ_t m_mapq;
        CIGAR_t m_cigar;
        RNAME_t m_mrnm;
        POS_t m_mpos;
        ISIZE_t m_isize;
        SEQ_t m_seq;
        QUAL_t m_qual;

    };

    class SamHandler {

    };
}


#endif //PROTAL_SAMHANDLER_H
