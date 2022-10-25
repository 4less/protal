//
// Created by fritsche on 24/10/22.
//

#pragma once

#include <cstdint>
#include <cstddef>

namespace protal {
    class SNP {
    public:
        using SNPFlag = char;
        using Orientation = bool;
        using Base = char;
        using FirstRead = bool;
        using ReadId = size_t;
        using SSize = uint16_t;
        using Qual = uint8_t;
        using TaxId = uint32_t;
        using GeneId = uint32_t;
        using SNPPos = uint32_t;

        SNPFlag type = 0;                     //8bit
        Orientation orientation = false;      //8bit
        Base reference = 0;                   //8bit
        Base variant = 0;                     //8bit //32
        Qual quality = 0;                     //8bit
        FirstRead first_read = false;         //8bit
        SSize structural_size = 01;           //16bit //64
        ReadId readid = 0;                    //64bit //128
        TaxId taxid = 0;                      //32bit
        GeneId geneid = 0;                    //32bit //192
        SNPPos snp_pos = 0;                   //32bit

//        SNP(SNPFlag type, Orientation orientation, Base reference, Base qual, FirstRead isfirst,
//            SSize ssize, ReadId readid, TaxId taxid, GeneId geneid, SNPPos snppos) :
//                m_type(type),
//                m_orientation(orientation),
//                m_reference(reference),
//                m_quality(qual),
//                m_first_read(isfirst),
//                m_structural_size(ssize),
//                m_readid(readid),
//                m_taxid(taxid),
//                m_geneid(geneid),
//                m_snp_pos(snppos) {}

        size_t inline Key() const {
            return (static_cast<size_t>(geneid) << 32ll) | static_cast<size_t>(snp_pos);
        }

        static inline size_t Key(GeneId geneid, SNPPos snppos) {
            return (static_cast<size_t>(geneid) << 32ll) | static_cast<size_t>(snppos);
        }

        SNP(){}
    };
}