//
// Created by fritsche on 10/11/22.
//

#pragma once

//#include "AlignmentOutputHandler.h"
#include <iostream>
#include "assert.h"
#include "SamHandler.h"
#include "Utils.h"
#include "Constants.h"
#include <tuple>

namespace protal {
    struct AlignmentInfo {
        uint16_t matches = 0;
        uint16_t mismatches = 0;
        uint16_t deletions = 0;
        uint16_t insertions = 0;
        uint16_t insertion_blocks = 0;
        uint16_t deletion_blocks = 0;
        uint16_t softclips = 0;
        uint16_t hardclips = 0;

        int alignment_score = 0;
        int gene_alignment_start = 0;
        uint32_t alignment_length = 0;
        float alignment_ani = 0.0f;

        std::string cigar;
        std::string compressed_cigar;

        // Deprecated
        uint32_t read_start_offset = 0;
        uint16_t cigar_start = 0;
        uint16_t cigar_end = 0;
        uint16_t cigar_start_read = 0;
        uint16_t cigar_end_read = 0;
        uint16_t cigar_start_gene = 0;
        uint16_t cigar_end_gene = 0;

        void Reset() {
            read_start_offset = 0;
            gene_alignment_start = 0;
            alignment_length = 0;
            alignment_ani = 0.0f;

            matches = 0;
            mismatches = 0;
            deletions = 0;
            insertions = 0;
            insertion_blocks = 0;
            deletion_blocks = 0;
            softclips = 0;
            hardclips = 0;
            alignment_score = 0;

            cigar_start = 0;
            cigar_end = 0;

            cigar.clear();
            compressed_cigar.clear();
        }

        bool Valid(size_t read_length) {
            auto len = std::count_if(cigar.begin(), cigar.end(), [](char c) {
                return c != 'D';
            });
            return len != read_length;
        }

        void GetInstructionCountsAndCompress() {
            ResetCigarStats();
            compressed_cigar.clear();

            bool swap_indel = false;

            char last_instruction = ' ';
            size_t instruction_counter = 0;

            for (auto i = 0; i < cigar.length(); i++) {
                auto c = cigar[i];
                // Compressed cigar
                if (c != last_instruction) {
                    if (instruction_counter > 0) {
                        insertion_blocks += last_instruction == 'I';
                        deletion_blocks += last_instruction == 'D';
                        if (swap_indel) {
                            if (last_instruction == 'D') last_instruction = 'I';
                            else if (last_instruction == 'I') last_instruction = 'D';
                        }
                        compressed_cigar += std::to_string(instruction_counter) + last_instruction;
                    }
                    last_instruction = c;
                    instruction_counter = 1;
                } else {
                    instruction_counter++;
                }
                // compressed cigar end
                insertions += c == 'I';
                deletions += c == 'D';
                matches += c == 'M';
                mismatches += c == 'X';
                softclips += c == 'S';
                hardclips += c == 'H';
            }

            // Get remainder of operations
            if (instruction_counter > 0) {
                insertion_blocks += last_instruction == 'I';
                deletion_blocks += last_instruction == 'D';
                if (swap_indel) {
                    if (last_instruction == 'D') last_instruction = 'I';
                    else if (last_instruction == 'I') last_instruction = 'D';
                }
                compressed_cigar += std::to_string(instruction_counter) + last_instruction;
            }
            alignment_length = cigar.length() - softclips;
        }


        void ResetCigarStats() {
            matches = 0;
            mismatches = 0;
            deletions = 0;
            insertions = 0;
            insertion_blocks = 0;
            deletion_blocks = 0;
            softclips = 0;
            hardclips = 0;
        }

        std::string ToString() const {
            std::string str = "Read Offset: " + std::to_string(read_start_offset);
            str += ", Gene Start: " + std::to_string(gene_alignment_start);
            str += ", AL: " + std::to_string(alignment_length);
            str += ", AS: " + std::to_string(Score());
            str += ", ANI: " + std::to_string(alignment_ani);
            str += ", matches: " + std::to_string(matches);
            str += ", mismatches: " + std::to_string(mismatches);
            str += ", deletions: " + std::to_string(deletions);
            str += ", insertions: " + std::to_string(insertions);
            str += ", insertion_blocks: " + std::to_string(insertion_blocks);
            str += ", deletion_blocks: " + std::to_string(deletion_blocks);
            str += ", softclips: " + std::to_string(softclips);
            str += ", hardclips: " + std::to_string(hardclips);
            str += ", cigar_start_read: " + std::to_string(cigar_start_read);
            str += ", cigar_end_read: " + std::to_string(cigar_end_read);
            str += ", cigar_start_ref: " + std::to_string(cigar_start_gene);
            str += ", cigar_end_ref: " + std::to_string(cigar_end_gene);
            return str;
        }
//
//        std::string GetClippedCigar() const {
//            return cigar.substr(cigar_start, cigar_end - cigar_start);
//        }
//
//        size_t GetClippedCigarLength() const {
//            return cigar_end - cigar_start;
//        }

        int GetEffectiveAlignmentLength() const {
            return cigar.length() - softclips;
        }

        int Score(int match_score = 0, int mismatch_penalty = 4, int gap_open_penalty = 6, int gap_extend_penalty = 2) const {
            return (matches * match_score) - (mismatch_penalty * mismatches + gap_open_penalty * (insertion_blocks + deletion_blocks) + gap_extend_penalty * (insertions + deletions));
        }

        void UpdateScore(int match_score = 0, int mismatch_penalty = 4, int gap_open_penalty = 6, int gap_extend_penalty = 2) {
            alignment_score = Score(match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty);
        }

        double GetProxyANI(int match_score = 0, int mismatch_penalty = 4, int gap_open_penalty = 6, int gap_extend_penalty = 2) const {
            auto score = match_score == 0 ?
                    1.0 + (static_cast<double>(Score(match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty)) / (mismatch_penalty * GetEffectiveAlignmentLength())) :
                         static_cast<double>(Score(match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty)) / (GetEffectiveAlignmentLength() * match_score);
            return static_cast<double>(score);
        }
    };

    class AlignmentResult {
        static constexpr int DEFAULT_VALUE = UINT32_MAX;

        uint32_t m_taxid = DEFAULT_VALUE;
        uint32_t m_geneid = DEFAULT_VALUE;
        int32_t m_genepos = INT32_MAX;

        bool m_forward = false;

        AlignmentInfo m_info;

    public:
        AlignmentResult(int alignment_score, uint32_t taxid, uint32_t geneid, int32_t genepos, bool forward) :
                m_taxid(taxid),
                m_geneid(geneid),
                m_genepos(genepos),
                m_forward(forward) {};

        AlignmentResult() {};

        std::string& Cigar() {
            return m_info.cigar;
        }

        std::string Cigar() const {
            return m_info.cigar;
        }

        uint32_t Taxid() const {
            return m_taxid;
        }

        uint32_t GeneId() const {
            return m_geneid;
        }

        int32_t GenePos() const {
            return m_genepos;
        }

        bool Forward() const {
            return m_forward;
        }

        int AlignmentScore() const {
            return m_info.Score();
        }

        AlignmentInfo & GetAlignmentInfo() {
            return m_info;
        }

        const AlignmentInfo& GetAlignmentInfo() const {
            return m_info;
        }

        void Set(size_t taxid, size_t geneid, int32_t abs_pos, bool forward) {
            m_taxid = taxid;
            m_geneid = geneid;
            m_genepos = abs_pos;
            m_forward = forward;
        }

        bool IsSet() const {
            return m_taxid != DEFAULT_VALUE;
        }

        void Reset() {
            m_taxid = DEFAULT_VALUE;
        }

        std::string ToString() const {
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

        std::string compressed_cigar;

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


    static void ExtractCigarInfo(std::string cigar, CigarInfo& info) {
        info.Reset();
        if (cigar.empty()) return;

        int instruction_counter = 0;
        char last_instruction = ' ';

        info.compressed_cigar = "";

        for (auto i = 0; i < cigar.length(); i++) {
            char c = cigar[i];
            // Compressed cigar
            if (cigar[i] != last_instruction) {
                if (instruction_counter > 0) {
                    info.compressed_cigar += std::to_string(instruction_counter) + last_instruction;
                }
                last_instruction = cigar[i];
                instruction_counter = 1;
            } else {
                instruction_counter++;
            }
            // compressed cigar end

            bool insertion = c == 'I';
            bool deletion = c == 'D';
            bool match = c == 'M';
            bool mismatch = c == 'X';
            bool hardclipped = c == 'H';
            bool softclipped = c == 'S';

            info.insertions += insertion;
            info.deletions += deletion;
            info.matches += match;
            info.hardclipped += hardclipped;
            info.softclipped += softclipped;
        }
        auto alignment_length = cigar.length() + info.deletions - info.insertions;
        info.unclipped_alignment_length = alignment_length;
        info.clipped_alignment_length = alignment_length - info.softclipped - info.hardclipped;
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


    static std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> GetLeadingTrailingIndels(std::string const& cigar) {
        uint32_t del_left = 0;
        uint32_t del_right = 0;
        uint32_t ins_left = 0;
        uint32_t ins_right = 0;

        while (cigar[del_left] == 'D') del_left++;
        while (cigar[cigar.length() - 1 - del_right] == 'D') del_right++;
        while (cigar[ins_left] == 'I') ins_left++;
        while (cigar[cigar.length() - 1 - ins_right] == 'I') ins_right++;
        return { del_left, del_right, ins_left, ins_right };
    }

    /**
     * Turn leading or trailing bases into softclip characters. With WFA this is necessary,
     * as WFA never outputs S or H for Softclipped/Hardclipped respectively.
     * @param cigar (uncompressed)
     * @param left (leading deletions)
     * @param right (trailing deletions)
     */
    static void Softclip(std::string& cigar, int left, int right) {
        for (auto i = 0; i < left; i++) cigar[i] = 'S';
        for (auto i = 0; i < right; i++) cigar[cigar.length() - 1 - i] = 'S';
    }

    /**
     * Hardclipping is only to clip away leading and trailing
     * Insertions as they do not belong to the core read. These insertions and deletions occur
     * when more reference is provided for an ends-free alignment to allow for insertions and deletions
     * within the alignment.
     * @param cigar (uncompressed)
     * @param left (leading insertions)
     * @param right (trailing insertions)
     * @return New clipped cigar string
     */
    static std::string Hardclip(std::string const& cigar, int left, int right) {
        return cigar.substr(left, cigar.length() - left - right);
    }

    static void PostProcessAlignment(std::string const& cigar, AlignmentInfo& info, size_t read_length, size_t total_reference_length, int reference_start_pos, int reference_end_pos, int relative_read_position) {
        auto [del_left, del_right, ins_left, ins_right] = GetLeadingTrailingIndels(cigar);

//        size_t new_cigar_start = ins_left;
//        size_t new_cigar_end = cigar.length() - ins_right;

//        info.cigar = Hardclip(cigar, ins_left, ins_right);
//        Softclip(info.cigar, del_left, del_right);

        info.gene_alignment_start = reference_start_pos + ins_left;

        size_t new_cigar_start = del_left;
        size_t new_cigar_end = cigar.length() - del_right;

        info.cigar = Hardclip(cigar, del_left, del_right);
        Softclip(info.cigar, ins_left, ins_right);

        info.gene_alignment_start = reference_start_pos + del_left;
//        std::cout << "GENE ALIGNMENT START " << info.gene_alignment_start << " from: " << reference_start_pos << " " << del_left << std::endl;
        info.GetInstructionCountsAndCompress();
        info.alignment_ani = (static_cast<double>(info.matches) / (info.matches + info.mismatches + info.insertions + info.deletions));
    }

    static AlignmentInfo GetAlignmentInfo(std::string const& cigar){
        AlignmentInfo info;
        info.cigar = cigar;
        info.GetInstructionCountsAndCompress();
        return info;
    }

    static AlignmentInfo GetAlignmentInfo(std::string&& cigar){
        AlignmentInfo info;
        info.cigar = cigar;
        info.GetInstructionCountsAndCompress();
        return info;
    }

    static void PrintAlignment(std::string const& cigar, std::string& query, std::string const& reference, std::ostream &out = std::cout) {
        int qpos = 0;
        int rpos = 0;
//        std::cout << sam.ToString() << std::endl;
//        std::cout << "Reflength: " << reference.length() << std::endl;
//        std::cout << sam.m_seq << std::endl;
//        std::cout << reference.substr(rpos, 50) << std::endl;



        std::string query_str = "";
        std::string ref = "";

        std::string align = "";
        std::string cigar_str = "";
        int cpos = 0;
        bool has_indel = false;
        bool faulty = false;
        auto i = 0;
        for (auto c : cigar) {
            align += c == 'M' ? '|' : ' ';
            cigar_str += c;

            if (c == 'M') {
                if (query[qpos + i] != reference[rpos+i]) {
                    faulty = true;
                }
            }

            if (c == 'I') {
                query_str += '-';
                ref += reference[rpos];
                rpos++;
                has_indel = true;
            } else if (c == 'D' || c == 'S') {
                query_str += query[qpos];
                ref += '-';
                qpos++;

                has_indel = true;
            } else {
                query_str += query[qpos];
                ref += reference[rpos];

                qpos++;
                rpos++;
            }
            i++;
        }
        if (faulty) {
            exit(3);
            Utils::Input();
        }
    }

    static bool PrintAlignment(SamEntry const& sam, std::string const& reference, std::ostream &out = std::cout) {
        int qpos = 0;
        int rpos = sam.m_pos - 1;

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

            if (op == 'D') {
                query += std::string(count, '-');
                ref += reference.substr(rpos, count);
                rpos += count;
                hasid = true;
            } else if (op == 'I' || op == 'S') {
                query += sam.m_seq.substr(qpos, count);
                ref += std::string(count, '-');
                qpos += count;
                // Maybe it can fix the problems
                //rpos += count;

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
            std::cout << "Faulty within print" << std::endl;
            return false;
        }
        return true;
    }

    static int Bitscore(PairedAlignment const& a) {
        auto score1 = a.first.IsSet() ? a.first.GetAlignmentInfo().Score(2, 3, 1, 2) : 0;
        auto score2 = a.second.IsSet() ? a.second.GetAlignmentInfo().Score(2, 3, 1, 2) : 0;
        return score1 + score2;
    }

    int MAPQv1(int s1, int s2) {
        //Assumes alignments come sorted and best is on top
        return 40.0 * (1.0-static_cast<double>(s2)/s1) * log10(static_cast<double>(s1));
    }

    int MAPQv2(int s1, int s2) {
        //Assumes alignments come sorted and best is on top
        if (s1 ==  s2) return 0;
        return 1 + (40.0 * (1.0-static_cast<double>(s2)/s1) * log10(static_cast<double>(s1)));
    }

    int MAPQv1(PairedAlignmentResultList& paired_alignment_results) {
        auto best_score = Bitscore(paired_alignment_results[0]);
        auto second_best_score = paired_alignment_results.size() > 1 ? Bitscore(paired_alignment_results[1]) : 0;

        //Assumes alignments come sorted and best is on top
        return MAPQv2(best_score, second_best_score);
    }

    std::tuple<int, int, int> MAPQv1Debug(PairedAlignmentResultList& paired_alignment_results) {
        auto best_score = Bitscore(paired_alignment_results[0]);
        auto second_best_score = paired_alignment_results.size() > 1 ? Bitscore(paired_alignment_results[1]) : 0;

        auto mapq = MAPQv2(best_score, second_best_score);
        //Assumes alignments come sorted and best is on top
        return {mapq, best_score, second_best_score};
    }
}