//
// Created by fritsche on 11/11/22.
//

#pragma once
#include "SNP.h"
#include "SamHandler.h"
#include "AlignmentUtils.h"
#include "Variant.h"
#include "robin_map.h"

namespace protal {
    using SNPList = std::vector<SNP>;

    using VariantBin = std::vector<Variant>;
    using Variants = tsl::robin_map<uint32_t, VariantBin>;

    static bool ExtractSNPs(SamEntry const& sam, std::string const& reference, SNPList &snps, int taxid, int geneid, int read_id = 0) {
//        snps.clear();
        int qpos = 0;
        int rpos = sam.m_pos - 1;

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
                    if (sam.m_seq[qpos + i] != 'N' && reference[rpos+i] != 'N' && sam.m_seq[qpos + i] != reference[rpos+i]) {
                        std::cerr << sam.m_seq[qpos + i] << "-" << reference[rpos+i] << std::endl;
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

            qpos += (op != 'D') * count;
            rpos += (!(op == 'I' || op == 'S')) * count;

//            if (faulty) exit(123);
        }

//        if (faulty) {
//            std::cerr << "FAULTY ---------------------------------" << std::endl;
//            PrintAlignment(sam, reference, std::cerr);
//            for (auto snp : snps) {
//                std::cerr << snp.ToString() << std::endl;
//            }
//            std::cerr << sam.m_seq << std::endl;
//            std::cerr << "Faulty sam: \n" << sam.ToString() << std::endl;
//            return false;
//        }
        return !faulty;
    }


    static bool IsAlignmentValid(AlignmentInfo const& info, std::string const& query, std::string const& reference) {
//        snps.clear();
        int qpos = 0;
        int rpos = info.gene_alignment_start;

        int count = 0;
        char op = ' ';


        int cpos = 0;
        bool faulty = false;

        while (NextCompressedCigar(cpos, info.compressed_cigar, count, op)) {
            if (op == 'M') {
                for (auto i = 0; i < count; i++) {
                    if (qpos+i < 0 || qpos+i >= query.size()) {
                        std::cerr << "Access out of bounds in IsAlignmentValid: " << qpos+i << " read len: " << query.size() << std::endl;
                    }
                    if (rpos+i < 0 || rpos+i >= reference.size()) {
                        std::cerr << "Access out of bounds in IsAlignmentValid: " << rpos+i << " read len: " << reference.size() << std::endl;
                    }

                    if (query[qpos + i] != 'N' && reference[rpos+i] != 'N' && query[qpos + i] != reference[rpos+i]) {
                        std::cerr << query[qpos + i] << "-" << reference[rpos+i] << '\t' << qpos+i << "-" << rpos+i << '\t' << query.size() << "-" << reference.size() << std::endl;
                        faulty = true;
                    }
                }
            }

            qpos += (op != 'D') * count;
            rpos += (!(op == 'I' || op == 'S')) * count;

        }

        return !faulty;
    }


}