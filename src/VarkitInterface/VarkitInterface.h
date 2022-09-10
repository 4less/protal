//
// Created by fritsche on 01/09/22.
//

#pragma once

#include <string>
#include <cstring>

namespace protal {


    struct ClassificationLine {
        static constexpr int GENEID_UNSET = 0;
        static constexpr uint32_t TAXID_UNSET = 0;
        static constexpr int GENEPOS_UNSET = INT32_MIN;
        uint32_t record_id = 0;
        uint32_t mutation_count = 0;

        uint32_t taxid = TAXID_UNSET;
        int geneid = GENEID_UNSET;
        int genepos = GENEPOS_UNSET;

        uint32_t total_hits_best = 0;
        uint32_t total_hits = 0;
        uint32_t leaf_hits_best = 0;
        uint32_t read_length = 0;
        uint32_t overlap_with_gene = 0;

        double rank_confidence = 0;
        double candidate_confidence = 0;
        double classification_confidence = 0;

        char classified = 'U';
        bool forward = false;
        uint8_t read_num = 0;

        double predicted_ani = -1;

        static const size_t buffer_size = 20;
        char buffer[buffer_size];

        void SetHeader(std::string &header) {
            size_t len = buffer_size - 1;
            if (header.length() < buffer_size - 1)
                len = header.length();
            std::memcpy(buffer, header.c_str(), len);
            buffer[len] = '\0';
        }

        void Reset() {
            record_id = 0;
            mutation_count = 0;
            taxid = TAXID_UNSET;
            geneid = GENEID_UNSET;
            genepos = GENEPOS_UNSET;
            total_hits_best = 0;
            total_hits = 0;
            leaf_hits_best = 0;
            read_length = 0;
            overlap_with_gene = 0;
            rank_confidence = 0.0;
            candidate_confidence = 0.0;
            classification_confidence = 0.0;
            predicted_ani = 0;
            classified = 'U';
            read_num = 0;
        }

        bool HasGene() const {
            return geneid != GENEID_UNSET;
        }
        bool HasPos() const {
            return genepos != GENEPOS_UNSET;
        }

        std::string ToString() {
            return (std::to_string(record_id) + " mutations: "
                    + std::to_string(mutation_count) + " taxid: "
                    + std::to_string(taxid) + " geneid: "
                    + std::to_string(geneid) + " genepos: "
                    + std::to_string(genepos) + " totalhitsbest: "
                    + std::to_string(total_hits_best) + " leafhitsbest: "
                    + std::to_string(leaf_hits_best) + " totalhits: "
                    + std::to_string(total_hits) + " rl: "
                    + std::to_string(read_length) + " ov: "
                    + std::to_string(overlap_with_gene) + " header: "
                    + std::string(buffer));



        }

    };

    /**
     * Mandatory fields used by varkits Profiler:
     * - taxid
     * - geneid
     * - genepos
     * - read_length
     * - total_hits_best (for read filter. Just set to high)
     * - (overlap with gene gets computed)
     */
    static void ToClassificationLine(ClassificationLine& line, size_t taxid, size_t geneid, size_t genepos, size_t read_length, size_t record_id, double ani, size_t total_hits_best=15) {
        line.Reset();
        line.taxid = taxid;
        line.geneid = geneid;
        line.genepos = genepos;
        line.read_length = read_length;
        line.record_id = record_id;
        line.predicted_ani = ani;
        line.total_hits_best = total_hits_best;
    }
}