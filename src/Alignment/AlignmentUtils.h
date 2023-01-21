//
// Created by fritsche on 10/11/22.
//

#pragma once

#include <iostream>
#include "assert.h"
#include "SamHandler.h"
namespace protal {

    static bool NextCompressedCigar(int& pos, const std::string& ccigar, int& count, char& operation) {
        if (pos == ccigar.length()) return false;
        if (!std::isdigit(ccigar[pos])) {
            std::cerr << "Malformed compressed cigar " << ccigar << std::endl;
            exit(90);
        }

        int cpos = pos;

        while (cpos < ccigar.length() && std::isdigit(ccigar[cpos])) {
            cpos++;
        }


        count = std::stoi(ccigar.substr(pos, cpos - pos));
        operation = ccigar[cpos];
        pos = cpos+1;

        return true;
    }

    struct CigarInfo {
        int32_t clipped_alignment_length = 0;
        int32_t unclipped_alignment_length = 0;
        int16_t insertions_open = 0;
        int16_t insertions = 0;
        int16_t deletions_open = 0;
        int16_t deletions = 0;
        int16_t softclipped = 0;
        int16_t hardclipped = 0;
        int16_t mismatches = 0;
        int16_t matches = 0;

        void Reset() {
            clipped_alignment_length = 0;
            unclipped_alignment_length = 0;
            insertions_open = 0;
            insertions = 0;
            deletions_open = 0;
            deletions = 0;
            softclipped = 0;
            hardclipped = 0;
            mismatches = 0;
            matches = 0;
        }

        double Ani() const {
            if (clipped_alignment_length == 0) {
                std::cout << ToString() << std::endl;
            }
            return static_cast<double>(matches + insertions + deletions) / static_cast<double>(clipped_alignment_length);
        }

        int Score(int mismatch_penalty = 4, int gap_open_penalty = 6, int gap_extend_penalty = 2) const {
            return -(mismatch_penalty * mismatches + gap_open_penalty * (insertions_open + deletions_open) + gap_extend_penalty * (insertions + deletions));
        }

        std::string ToString() const {
            std::string str;
            str += "clipped_alignment_length: " + std::to_string(clipped_alignment_length) + '\n';
            str += "unclipped_alignment_length: " + std::to_string(unclipped_alignment_length) + '\n';
            str += "insertions_open: " + std::to_string(insertions_open) + '\n';
            str += "insertions: " + std::to_string(insertions) + '\n';
            str += "deletions_open: " + std::to_string(deletions_open) + '\n';
            str += "deletions: " + std::to_string(deletions) + '\n';
            str += "softclipped: " + std::to_string(softclipped) + '\n';
            str += "hardclipped: " + std::to_string(hardclipped) + '\n';
            str += "mismatches: " + std::to_string(softclipped) + '\n';
            str += "matches: " + std::to_string(hardclipped);
            return str;
        }
    };



    struct AlignmentInfo {
        uint32_t read_start_offset = 0;
        uint32_t gene_start_offset = 0;
        uint32_t alignment_length = 0;
        float alignment_ani = 0.0f;

        uint16_t mismatches = 0;
        uint16_t deletions = 0;
        uint16_t insertions = 0;
        uint16_t iblocks = 0;
        uint16_t dblocks = 0;
        int16_t alignment_score = 0;
        int16_t alignment_qscore = 0;
        int16_t mismatch_quality_sum = 0;

        std::string cigar;
        std::string compressed_cigar;

        std::string ToString() const {
            std::string str = "Read Offset: " + std::to_string(read_start_offset);
            str += ", Gene Start: " + std::to_string(gene_start_offset);
            str += ", AL: " + std::to_string(alignment_length);
            str += ", AS: " + std::to_string(alignment_score);
            str += ", QS: " + std::to_string(alignment_qscore);
            str += ", MMQS: " + std::to_string(mismatch_quality_sum);
            str += ", ANI: " + std::to_string(alignment_ani);
            str += ", CCIG: " + compressed_cigar;
            return str;
        }
    };

    struct AlignmentPair {
        std::optional<SamEntry> first;
        std::optional<SamEntry> second;

        AlignmentPair(std::optional<SamEntry>& first, std::optional<SamEntry>& second) :
                first(first),
                second(second) {
        }

        AlignmentPair(std::optional<SamEntry>&& first, std::optional<SamEntry>&& second) :
                first(first),
                second(second) {
        }

        AlignmentPair() {};

        bool HasFirst() {
            return first.has_value();
        }

        bool HasSecond() {
            return second.has_value();
        }

        bool IsPair() {
            return HasFirst() && HasSecond();
        }

        bool IsSingle() {
            return HasFirst() ^ HasSecond();
        }

        SamEntry& First() {
            assert(first.has_value());
            return first.value();
        }

        SamEntry& Second() {
            assert(second.has_value());
            return second.value();
        }

        SamEntry& Single() {
            return HasFirst() ? First() : Second();
        }

        SamEntry& Any() {
            return Single();
        }

        std::pair<SamEntry&,SamEntry&> Pair() {
            return { First(), Second() };
        }
    };

    static double CigarANI(std::string cigar) {
        if (cigar.empty()) return 0;
        size_t matches = std::count_if(cigar.begin(), cigar.end(), [](char const &c) {
            return c == 'M';
        });
        size_t mismatches = std::count_if(cigar.begin(), cigar.end(), [](char const &c) {
            return c == 'X';
        });
        size_t indel = std::count_if(cigar.begin(), cigar.end(), [](char const &c) {
            return c == 'I' || c == 'D';
        });

        return (static_cast<double>(matches) / (matches + mismatches + indel));
    }

    static void CompressedCigarInfo(std::string cigar, CigarInfo& info) {
        info.Reset();
        if (cigar.empty()) return;

        int count;
        char c;
        int digit_start = -1;

        int total_count = 0;
        for (auto i = 0; i < cigar.length(); i++) {
            if (!std::isdigit(cigar[i])) {
                count = stoi(cigar.substr(digit_start, i-digit_start));
                c = cigar[i];
                digit_start = -1;
                if (c == 'S') info.softclipped += count;
                else if (c == 'H') info.hardclipped += count;
                else if (c == 'X') info.mismatches += count;
                else if (c == 'M') info.matches += count;
                else if (c == 'I') {
                    info.insertions += count;
                    info.insertions_open++;
                }
                else if (c == 'D') {
                    info.deletions += count;
                    info.deletions_open++;
                }
                total_count += count;
            } else if (digit_start == -1) {
                digit_start = i;
            }
        }

        info.clipped_alignment_length = total_count - info.softclipped;
        info.unclipped_alignment_length = total_count;

//            std::cout << total_count << " " << info.Score() << " " << info.Ani() << std::endl;
    }

    static int CompressedCigarScore(std::string cigar, int mismatch_pen = 4, int gapopen_pen = 6, int gapext_pen = 2, int match_score = 0) {
        if (cigar.empty()) return 0;

        int score = 0;
        int count;
        char c;
        int digit_start = -1;


        for (auto i = 0; i < cigar.length(); i++) {
            if (!std::isdigit(cigar[i])) {
                count = stoi(cigar.substr(digit_start, i-digit_start));
                c = cigar[i];
                digit_start = -1;
                if (c == 'X') score -= mismatch_pen * count;
                if (c == 'M') score += match_score;
                if (c == 'I' || c == 'D') score -= (gapopen_pen + (count-1) * gapext_pen);
            } else if (digit_start == -1) {
                digit_start = i;
            }
        }
        return score;
    }

    static int CigarScore(std::string cigar, int mismatch_pen = 4, int gapopen_pen = 6, int gapext_pen = 2, int match_score=0) {
        int score = 0;
        char last = 'M';
        for (char c : cigar) {
            if (c == 'M') score += match_score;
            if (c == 'X') score -= 4;
            if (c == 'I' && last == 'I') score -= gapext_pen;
            if (c == 'I' && last != 'I') score -= gapopen_pen;
            if (c == 'D' && last == 'D') score -= gapext_pen;
            if (c == 'D' && last != 'D') score -= gapopen_pen;
            last = c;
        }
        return score;
    }

    static void GetExtendedAlignmentInfo(AlignmentInfo& info, std::string const& cigar, std::string& sequence, std::string& quality, int dove_left_read=0, int dove_left_reference=0, int dove_right_read=0, int dove_right_reference=0, int mismatch_pen = 4, int mismatch_pen_min = 2, int mismatch_pen_max = 6, int gapopen_pen = 6, int gapext_pen = 2) {
        size_t alignment_start = 0;
        size_t alignment_length = 0;
        size_t indel_count = 0;
        size_t match_count = 0;
        size_t mismatch_count = 0;

        info.compressed_cigar = "";
        info.deletions = 0;
        info.insertions = 0;

        int cigar_start_read = 0;
        while (cigar[cigar_start_read] == 'D' && cigar_start_read < dove_left_read) {
            cigar_start_read++;
            assert(cigar_start_read < cigar.length());
        }

        int cigar_start_gene = 0;
        while (cigar[cigar_start_gene] == 'I' && cigar_start_gene < dove_left_reference) {
            cigar_start_gene++;
            assert(cigar_start_gene < cigar.length());
        }

        int cigar_end_read = cigar.length() - 1;
        assert(cigar_end_read < cigar.length() && cigar_end_read >= 0);
        while (cigar[cigar_end_read] == 'D' && cigar_end_read >= (cigar.length() - dove_right_read)) {
            cigar_end_read--;
            assert(cigar_end_read < cigar.length() && cigar_end_read >= 0);
        }

        int cigar_end_gene = cigar.length() - 1;
//            std::cout << "Old cigar end: " << cigar_end_gene << std::endl;
        assert(cigar_end_gene < cigar.length() && cigar_end_gene >= 0);
        while (cigar[cigar_end_gene] == 'I' && cigar_end_gene >= (cigar.length() - dove_right_reference)) {
            cigar_end_gene--;
            assert(cigar_end_gene < cigar.length() && cigar_end_gene >= 0);
        }
//            std::cout << "New cigar end: " << cigar_end_gene << std::endl;

        cigar_end_read++;
        cigar_end_gene++;

        int cigar_start = std::max(cigar_start_read, cigar_start_gene);
        int cigar_end = std::min(cigar_end_read, cigar_end_gene);
//            std::cout << cigar_start_read << ", "  << cigar_start_gene << ", " << cigar_end_read << ", " << cigar_end_gene << std::endl;
//            std::cout << "Cigar start: " << cigar_start << ",  Cigar end: " << cigar_end << std::endl;
        info.cigar = cigar.substr(cigar_start, cigar_end - cigar_start);

        // Get compressed cigar
        char last_instruction = ' ';
        size_t instruction_counter = 0;

        int mismatch_pen_diff = mismatch_pen_max - mismatch_pen_min;
        int ascore = 0;
        int qascore = 0;

        int mismatch_quality_sum = 0;

        info.compressed_cigar = "";
        for (auto i = 0; i < info.cigar.length(); i++) {
            char c = info.cigar[i];
            // Compressed cigar
            if (info.cigar[i] != last_instruction) {
                if (instruction_counter > 0) {
                    info.compressed_cigar += std::to_string(instruction_counter) + last_instruction;

                    if (last_instruction == 'X') {
                        ascore -= mismatch_pen * instruction_counter;
                        qascore -= mismatch_pen * instruction_counter;
                    }
                    if (last_instruction == 'I' || last_instruction == 'D') {
                        info.iblocks += (last_instruction == 'I');
                        info.dblocks += (last_instruction == 'D');
                        ascore -= (gapopen_pen + (instruction_counter-1) * gapext_pen);
                        qascore -= (gapopen_pen + (instruction_counter-1) * gapext_pen);
                    }
                }
                last_instruction = info.cigar[i];
                instruction_counter = 1;
            } else {
                instruction_counter++;
            }
            // compressed cigar end

            bool insertion = c == 'I';
            bool deletion = c == 'D';
            info.insertions += insertion;
            info.deletions += deletion;

            indel_count += deletion || insertion;
            match_count += c == 'M';
            mismatch_count  += c == 'X';

            if (c == 'X') {
                int pos = cigar_start_read + i;
                mismatch_quality_sum += quality[pos];
//                std::cout << "X: " << sequence[pos] << " " << quality[pos] << " " << std::to_string(quality[pos]) << std::endl;
            }

            alignment_length += deletion;
            alignment_length -= (insertion ||  c == 'H' || c == 'S');
            alignment_length++;
        }

        // Get remainder of operations
        if (instruction_counter > 0) {
            info.compressed_cigar += std::to_string(instruction_counter) + last_instruction;
            if (last_instruction == 'X') {
                ascore -= mismatch_pen * instruction_counter;
                qascore -= mismatch_pen * instruction_counter;
            }
            if (last_instruction == 'I' || last_instruction == 'D') {
                ascore -= (gapopen_pen + (instruction_counter-1) * gapext_pen);
                qascore -= (gapopen_pen + (instruction_counter-1) * gapext_pen);
            }
        }

        info.alignment_score = ascore;
        info.alignment_qscore = qascore;
        info.mismatch_quality_sum = mismatch_quality_sum;

        info.read_start_offset = cigar_start_read;
        info.gene_start_offset = cigar_start_gene;

        info.alignment_length = alignment_length;
        info.alignment_ani = (static_cast<double>(match_count) / (match_count + mismatch_count + indel_count));
    }


    static void GetAlignmentInfo(AlignmentInfo& info, std::string const& cigar, int dove_left_read=0, int dove_left_reference=0, int dove_right_read=0, int dove_right_reference=0, int mismatch_pen = 4, int mismatch_pen_min = 2, int mismatch_pen_max = 6, int gapopen_pen = 6, int gapext_pen = 2) {
        size_t alignment_start = 0;
        size_t alignment_length = 0;
        size_t indel_count = 0;
        size_t match_count = 0;
        size_t mismatch_count = 0;

        info.compressed_cigar = "";
        info.deletions = 0;
        info.insertions = 0;

        int cigar_start_read = 0;
        while (cigar[cigar_start_read] == 'D' && cigar_start_read < dove_left_read) {
            cigar_start_read++;
            assert(cigar_start_read < cigar.length());
        }

        int cigar_start_gene = 0;
        while (cigar[cigar_start_gene] == 'I' && cigar_start_gene < dove_left_reference) {
            cigar_start_gene++;
            assert(cigar_start_gene < cigar.length());
        }

        int cigar_end_read = cigar.length() - 1;
        assert(cigar_end_read < cigar.length() && cigar_end_read >= 0);
        while (cigar[cigar_end_read] == 'D' && cigar_end_read >= (cigar.length() - dove_right_read)) {
            cigar_end_read--;
            assert(cigar_end_read < cigar.length() && cigar_end_read >= 0);
        }

        int cigar_end_gene = cigar.length() - 1;
//            std::cout << "Old cigar end: " << cigar_end_gene << std::endl;
        assert(cigar_end_gene < cigar.length() && cigar_end_gene >= 0);
        while (cigar[cigar_end_gene] == 'I' && cigar_end_gene >= (cigar.length() - dove_right_reference)) {
            cigar_end_gene--;
            assert(cigar_end_gene < cigar.length() && cigar_end_gene >= 0);
        }
//            std::cout << "New cigar end: " << cigar_end_gene << std::endl;

        cigar_end_read++;
        cigar_end_gene++;

        int cigar_start = std::max(cigar_start_read, cigar_start_gene);
        int cigar_end = std::min(cigar_end_read, cigar_end_gene);
//            std::cout << cigar_start_read << ", "  << cigar_start_gene << ", " << cigar_end_read << ", " << cigar_end_gene << std::endl;
//            std::cout << "Cigar start: " << cigar_start << ",  Cigar end: " << cigar_end << std::endl;
        info.cigar = cigar.substr(cigar_start, cigar_end - cigar_start);

        // Get compressed cigar
        char last_instruction = ' ';
        size_t instruction_counter = 0;

        int mismatch_pen_diff = mismatch_pen_max - mismatch_pen_min;
        int ascore = 0;
        int qascore = 0;

        info.compressed_cigar = "";
        for (auto i = 0; i < info.cigar.length(); i++) {
            char c = info.cigar[i];
            // Compressed cigar
            if (info.cigar[i] != last_instruction) {
                if (instruction_counter > 0) {
                    info.compressed_cigar += std::to_string(instruction_counter) + last_instruction;

                    if (last_instruction == 'X') {
                        ascore -= mismatch_pen * instruction_counter;
                        qascore -= mismatch_pen * instruction_counter;
                    }
                    if (last_instruction == 'I' || last_instruction == 'D') {
                        ascore -= (gapopen_pen + (instruction_counter-1) * gapext_pen);
                        qascore -= (gapopen_pen + (instruction_counter-1) * gapext_pen);
                    }
                }
                last_instruction = info.cigar[i];
                instruction_counter = 1;
            } else {
                instruction_counter++;
            }
            // compressed cigar end

            bool insertion = c == 'I';
            bool deletion = c == 'D';
            info.insertions += insertion;
            info.deletions += deletion;

            indel_count += deletion || insertion;
            match_count += c == 'M';
            mismatch_count  += c == 'X';

            alignment_length += deletion;
            alignment_length -= (insertion ||  c == 'H' || c == 'S');
            alignment_length++;
        }

        // Get remainder of operations
        if (instruction_counter > 0) {
            info.compressed_cigar += std::to_string(instruction_counter) + last_instruction;
            if (last_instruction == 'X') {
                ascore -= mismatch_pen * instruction_counter;
                qascore -= mismatch_pen * instruction_counter;
            }
            if (last_instruction == 'I' || last_instruction == 'D') {
                ascore -= (gapopen_pen + (instruction_counter-1) * gapext_pen);
                qascore -= (gapopen_pen + (instruction_counter-1) * gapext_pen);
            }
        }

        info.alignment_score = ascore;
        info.alignment_qscore = qascore;

        info.read_start_offset = cigar_start_read;
        info.gene_start_offset = cigar_start_gene;

        info.alignment_length = alignment_length;
        info.alignment_ani = (static_cast<double>(match_count) / (match_count + mismatch_count + indel_count));
    }

    static AlignmentInfo GetAlignmentInfo(std::string const& cigar) {
        AlignmentInfo info;
        GetAlignmentInfo(info, cigar);
        return info;
    }

    static void PrintAlignment(SamEntry const& sam, std::string const& reference, std::ostream &out = std::cout) {
        int qpos = 0;
        int rpos = sam.m_pos;

        int count = 0;
        char op = ' ';
//        std::cout << sam.ToString() << std::endl;
//        std::cout << "Reflength: " << reference.length() << std::endl;
//        std::cout << sam.m_seq << std::endl;
//        std::cout << reference.substr(rpos, 50) << std::endl;


        std::string query = "";
        std::string ref = "";

        std::string align = "";
        std::string cigar = "";
        int cpos = 0;
        bool hasid = false;
        bool faulty = false;
        while (NextCompressedCigar(cpos, sam.m_cigar, count, op)) {
//            std::cout << op << " (" << count << ")" << std::endl;
            align += op == 'M' ? std::string(count, '|') : std::string(count, ' ');
            cigar += std::string(count, op);

            if (op == 'M') {
                for (auto i = 0; i < count; i++) {
                    if (sam.m_seq[qpos + i] != reference[rpos+i]) {
                        faulty = true;
                    }
                }
            }

            if (op == 'I') {
                query += std::string(count, '-');
                ref += reference.substr(rpos, count);
                rpos += count;
                hasid = true;
            } else if (op == 'D' || op == 'S') {
                query += sam.m_seq.substr(qpos, count);
                ref += std::string(count, '-');
                qpos += count;

                hasid = true;
            } else {
                query += sam.m_seq.substr(qpos, count);
                ref += reference.substr(rpos, count);

                qpos += count;
                rpos += count;
            }
        }
        out << "-------------------------" << std::endl;
        out << "Align: " << sam.m_qname << " to " << sam.m_rname << std::endl;
        out << "Reversed? " << sam.IsReversed() << std::endl;
        out << cigar << std::endl;
        out << query << std::endl;
        out << align << std::endl;
        out << ref << std::endl;
        out << "-------------------------" << std::endl;

        if (faulty) {
            exit(3);
            Utils::Input();
        }
    }

}