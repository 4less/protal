//
// Created by joachim on 08/06/2020.
//

#ifndef VARKIT_KMERUTILS_H
#define VARKIT_KMERUTILS_H


#include <cstdint>
#include <string>
#include <FastxReader.h>
#include <vector>
#include <iostream>

using namespace std;

namespace KmerUtils {
    static std::pair<size_t, size_t> ExtractHeaderInformation(std::string const& header) {
        int underscore = header.find('_');
        if (underscore < 0) {
            std::cerr << "Must contain underscore. " << header << std::endl;
            exit(9);
        }
        int start = header[0] == '>' ? 1 : 0;
        size_t taxid = std::stoul(header.substr(start, underscore - start));
        size_t geneid = std::stoul(header.substr(underscore+1, header.size() - underscore - 1));
        return { taxid, geneid };
    }

    static const std::pair<uint32_t, uint32_t> ExtractTaxIdGeneId(std::string &ref) {
        int delim_pos = -1;
        while (ref[++delim_pos] != '_');
        int next_delim_pos = delim_pos;
        while (next_delim_pos < ref.size() && ref[++next_delim_pos] != '_');
        return { stoul(ref.substr(0, delim_pos)), stoul(ref.substr(delim_pos+1, next_delim_pos - delim_pos - 1)) };
    }

    static void SortVarkitFasta(std::vector<FastxRecord>& fasta) {
        std::sort(fasta.begin(), fasta.end(), [](FastxRecord const& a, FastxRecord const& b) constexpr {
            auto [a_taxid, a_geneid] = ExtractHeaderInformation(a.header);
            auto [b_taxid, b_geneid] = ExtractHeaderInformation(b.header);
            if (a_taxid < b_taxid) return true;
            else if (b_taxid < a_taxid) return false;
            else {
                if (a_geneid < b_geneid) return true;
                else return false;
            }
            return true;
        });
    }

    static inline uint8_t getBitFromBase(char base) {
        switch (base) {
            case 'A':
                return (0);
            case 'C':
                return (1);
            case 'G':
                return (2);
            case 'T':
                return (3);
        }
        return (-1);
    }
    
    static inline uint8_t getBitFromBase(const char* base) {
        switch((char)(*base)) {
            case 'A':
                return (0);
            case 'C':
                return (1);
            case 'G':
                return (2);
            case 'T':
                return (3);
        }
        return(-1);
    }
    
    static inline uint8_t getBitFromBaseC(char base) {
        switch(base) {
            case 'A':
                return (3);
            case 'C':
                return (2);
            case 'G':
                return (1);
            case 'T':
                return (0);
        }
        return(-1);
    }
    
    static inline uint8_t getBitFromBaseC(const char* base) {
        switch((char)(*base)) {
            case 'A':
                return (3);
            case 'C':
                return (2);
            case 'G':
                return (1);
            case 'T':
                return (0);
        }
        return(-1);
    }
    
    static inline char getBaseFromBit(uint8_t bits) {
        switch(bits) {
            case 0:
                return ('A');
            case 1:
                return ('C');
            case 2:
                return ('G');
            case 3:
                return ('T');
        }
        return('#');
    }

    static inline uint64_t BaseToInt(char const& base) {
        switch(base) {
            case 'A':
                return (0);
            case 'C':
                return (1);
            case 'G':
                return (2);
            case 'T':
                return (3);
        }
        return(4);
    }

    static inline uint64_t BaseToIntC(char const& base) {
        switch(base) {
            case 'A':
                return (3);
            case 'C':
                return (2);
            case 'G':
                return (1);
            case 'T':
                return (0);
        }
        return(4);
    }

    static inline uint64_t BaseToInt(char const &base, int replace_n) {
        switch(base) {
            case 'A':
                return (0);
            case 'C':
                return (1);
            case 'G':
                return (2);
            case 'T':
                return (3);
        }
        return(replace_n);
    }

    static inline uint64_t BaseToIntC(char const& base, int replace_n) {
        switch(base) {
            case 'A':
                return (3);
            case 'C':
                return (2);
            case 'G':
                return (1);
            case 'T':
                return (0);
        }
        return(replace_n);
    }

    static inline char IntToBase(uint8_t bits) {
        switch(bits) {
            case 0:
                return ('A');
            case 1:
                return ('C');
            case 2:
                return ('G');
            case 3:
                return ('T');
        }
        return('#');
    }

    static std::string ToString(size_t kmer, size_t bits) {
        std::string result (bits/2, '_');
//        for (int i = bits/2 - 1; i >= 0; i--) {
        for (int i = 0; i < bits/2; i++)
            result[i] = KmerUtils::getBaseFromBit((kmer >> (2llu * (bits/2 - i - 1))) & 3llu);
        return result;
    }

    static bool IsBitSet(size_t n, size_t pos) {
        return (n >> (64 - pos - 1)) & 1;
    }

    static std::string ToBitString(size_t n) {
        std::string result = "";
        for (size_t i = 0; i < 64; i++) {
            result += to_string(IsBitSet(n, i));
        }
        return result;
    }

    static inline char Complement(char c) {
        switch(c) {
            case 'A': return ('T');
            case 'C': return ('G');
            case 'G': return ('C');
            case 'T': return ('A');
        }
        return ('N');
    }
    
    static std::string bytesToDNA(uint8_t *kmer, int len) {
        std::string s;
        for (int i = 0; i < len; i++) {
            s += getBaseFromBit( (kmer[i/4] >> (2*(3-(i%4)))) & 0x03);
        }
        return s;
    }
    
    std::string ReverseComplement(std::string forward) {
        std::string reverse = "";
        const char * seq = forward.c_str();
        for (int i = forward.length()-1; i >= 0; i--) {
            reverse += Complement(forward[i]);
        }
        return reverse;
    }
//    static std::string ReverseComplement(std::string &forward) {
//        std::string reverse = "";
//        const char * seq = forward.c_str();
//        for (int i = forward.length()-1; i >= 0; i--) {
//            reverse += Complement(forward[i]);
//        }
//        return reverse;
//    }

    static std::string ExpandShape(std::string kmer, const bool* shape, size_t shape_size) {
        std::string result = "";

        size_t kmer_pos = 0;
        for (int i = 0; i < shape_size; i++) {
            if (!shape[i]) {
                result += kmer[kmer_pos];
                kmer_pos++;
            } else {
                result += '_';
            }
        }
        return result;
    }
    
    static string dnaToBitString(string kmer) {
        const char* cseq = kmer.c_str();
        std::string int_string = "";
        for (int i = 0; i < kmer.length(); i++) {
            int_string += to_string(getBitFromBase(cseq[i]));
        }
        return int_string;
    }
    
    static void dnaToBytes(string kmer, uint8_t * key) {
        const char* cstr = kmer.c_str();
        
        for (int i = 0; i < kmer.length(); i++) {
            if (i % 4 == 0) key[i/4] = 0;
            key[i/4] = key[i/4] | KmerUtils::getBitFromBase(cstr[i]) << (2*(3-(i%4)));
        }
    }
}


#endif
