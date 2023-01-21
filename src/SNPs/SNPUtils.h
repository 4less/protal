//
// Created by fritsche on 11/11/22.
//

#pragma once
#include "SNP.h"
#include "SamHandler.h"
#include "AlignmentUtils.h"

namespace protal {
    using SNPList = std::vector<SNP>;

    class SNPNode {

    };

    class SNPNetwork {
        
    };

    static void ExtractSNPs(SamEntry const& sam, std::string const& reference, SNPList &snps, int taxid, int geneid, int read_id = 0) {
        snps.clear();
        int qpos = 0;
        int rpos = sam.m_pos;

        int count = 0;
        char op = ' ';


        std::string query = "";
        std::string ref = "";

        std::string align = "";
        std::string cigar = "";
        int cpos = 0;
        bool faulty = false;
        SNP snp;
        snp.taxid = taxid;
        snp.geneid = geneid;
        snp.readid = read_id;
        snp.orientation = sam.IsReversed();
        bool is_variant = false;

        while (NextCompressedCigar(cpos, sam.m_cigar, count, op)) {
            is_variant = false;
            if (op == 'M') {
                for (auto i = 0; i < count; i++) {
                    if (sam.m_seq[qpos + i] != reference[rpos+i]) {
                        faulty = true;
                    }
                }
            }

            if (op == 'I') {
                is_variant = true;
                snp.type = 2;
                snp.SetStructural(sam.m_seq.substr(qpos, count));
                snp.structural_size = count;
//                std::cout << "Insertion of " << count << " (" << sam.m_seq.substr(qpos, count) << ") at " << qpos << ", " << rpos << std::endl;
            } else if (op == 'D') {
                is_variant = true;
                snp.type = 4;
                snp.SetStructural(sam.m_seq.substr(qpos, count));
                snp.structural_size = count;
//                std::cout << "Deletion of  " << count << " (" << reference.substr(rpos, count) << ") at " << qpos << ", " << rpos << std::endl;
            } else if (op == 'X') {
                is_variant = true;
                snp.type = 1;
                snp.reference = reference[rpos];
                snp.snp_pos = rpos;
                snp.variant = sam.m_seq[qpos];
                snp.quality = sam.m_qual[qpos];
//                std::cout << sam.m_seq[qpos] << " -> " << reference[rpos] << " (" << qpos << ", " << rpos << ")" << std::endl;

            }
            if (is_variant) {
                snps.emplace_back(std::move(snp));
            }

            qpos += (op != 'I') * count;
            rpos += (!(op == 'D' || op == 'S')) * count;
        }
//        for (auto snp : snps) {
//            std::cout << snp.ToString() << std::endl;
//        }

        if (faulty) {
            exit(3);
        }
    }
}