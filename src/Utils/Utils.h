//
// Created by joachim on 20/07/2020.
//

#pragma once

#include <vector>
#include <string>
#include <iomanip>
#include <sys/stat.h>
#include <unordered_map>
#include <cstring>
#include <robin_set.h>
#include <iostream>
#include <emmintrin.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>

namespace Utils {
    class Histogram {
        std::vector<uint32_t> frequencies{0, 0};

    public:
        void Join(Histogram const& other) {
            if (frequencies.size() < other.frequencies.size()) {
                frequencies.resize(other.frequencies.size(), 0);
            }
            for (auto i = 0; i < other.frequencies.size(); i++) {
                frequencies[i] += other.frequencies[i];
            }
        }
        void AddObservation(uint32_t observation) {
            if (observation >= frequencies.size()) {
                frequencies.resize(observation+1, 0);
            }
            frequencies[observation]++;
        }

        void ToTSV(std::string& path) {
            std::ofstream out(path, std::ios::out);
            ToStream(out);
            out.close();
        }
        void ToTSV(std::string&& path) {
            std::ofstream out(path.c_str(), std::ios::out);
            ToStream(out);
            out.close();
        }

        void ToStream(std::ostream& os) {
            for (auto i = 0; i < frequencies.size(); i++) {
                os << i << '\t' << frequencies[i] << '\n';
            }
        }
    };

    static void shift(uint8_t * array, size_t length, size_t by) {
        for (int i = 0; i < length-1; i++) {
            array[i] <<= by;
            array[i] |= array[i + 1] >> (8-by);
        }
        array[length-1] <<= by;
    }

    static size_t nChoosek( size_t n, size_t k )
    {
        if (k > n) return 0;
        if (k * 2 > n) k = n-k;
        if (k == 0) return 1;

        int result = n;
        for( size_t i = 2; i <= k; ++i ) {
            result /= i;
            result *= (n-i+1);
        }
        return result;
    }

    static std::string Input() {
        std::string input;
        std::cout << "Input: ";
        std::cin >> input;
        return input;
    }

    struct pair_hash
    {
        template <class T1, class T2>
        std::size_t operator() (const std::pair<T1, T2> &pair) const {
            return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
        }
    };

    template<typename T>
    static void PrintVector(std::vector<T> vector, std::ostream& os) {
        if (vector.size() == 0) {
            os << "[]" << std::endl;
            return;
        }
        os << "[" << vector[0].ToString();

        for (auto i = 1; i < vector.size(); i++) {
            os << ", " << vector[i].ToString();
        }
        os << "]" << std::endl;
    }


    static std::string PrintU8(uint8_t* e, size_t size) {
        std::string res = "[";
        for (auto i = 0; i < size; i++) {
            res += std::to_string(i);
            if (size < size - 1)
                res += ", ";
        }
        res += "]";
        return res;
    }

    static std::string GetSharedPrefix(std::string a, std::string b) {
        for (int i = 0; i < std::min(a.size(), b.size()); i++)
            if (a[i] != b[i]) return a.substr(0, i);
        return "";
    }

    static bool HasSuffix(std::string s, std::string suffix) {
        if (suffix.length() > s.length()) return false;
        for (int i = 0; i < suffix.length(); i++) {
            auto str_idx = s.length() - 1 - i;
            auto suffix_idx = suffix.length() - 1 - i;
            if (s[str_idx] != suffix[suffix_idx]) return false;
        }
        return true;
    }

    static std::string RemoveSuffix(std::string s, std::string suffix) {
        if (suffix.length() > s.length()) return s;
        for (int i = 0; i < suffix.length(); i++) {
            auto str_idx = s.length() - 1 - i;
            auto suffix_idx = suffix.length() - 1 - i;
            if (s[str_idx] != suffix[suffix_idx]) return s;
        }
        return s.substr(0, s.length() - suffix.length());
    }

    static std::string GetDirectory(std::string s) {
        return s.substr(0, s.find_last_of('/'));
    }

    static std::string GetBasename(std::string s) {
        return s.substr(s.find_last_of('/')+1, s.length());
    }

    static std::string PadString(std::string s, size_t pad_to, char pad_with) {
        if (s.length() > pad_to) return s;
        return s + std::string(pad_to - s.length(), pad_with);
    }

    static std::string PrintU8(tsl::robin_set<uint8_t> &s) {
        std::string res = "[";
        size_t idx = 0;
        for (auto e : s) {
            res += std::to_string((uint32_t) e);
            if (idx < s.size() - 1)
                res += ", ";
            idx++;
        }
        res += "]";
        return res;
    }

    template <typename Range, typename Value = typename Range::value_type>
    std::string Join(Range const& elements, const char *const delimiter) {
        std::ostringstream os;
        auto b = std::begin(elements), e = std::end(elements);

        if (b != e) {
            std::copy(b, std::prev(e), std::ostream_iterator<Value>(os, delimiter));
            b = std::prev(e);
        }
        if (b != e) {
            os << *b;
        }

        return os.str();
    }


    //            typename std::enable_if<std::is_function<Func>::value>
    template <typename Range, typename Value = typename Range::value_type, typename Func>
    std::string Join(Range const& elements, const char *const delimiter, Func func) {
        std::ostringstream os;
        auto b = std::begin(elements), e = std::end(elements);

        if (b != e) {
            std::for_each(b, std::prev(e), [&os, &func, &delimiter](Value const& value) {
                os << func(value) << delimiter;
            });
            b = std::prev(e);
        }
        if (b != e) {
            os << func(*b);
        }

        return os.str();
    }

    static std::vector<std::string> split(const std::string& s, std::string delimiter)
    {
        std::vector<std::string> tokens;
        std::string token;
        
        std::string line = s;
        int lastpos = 0;
        int pos = 0;
        
        while ((pos = line.find(delimiter, lastpos)) != std::string::npos) {
            token = line.substr(lastpos, pos-lastpos);
            tokens.push_back(token);
            lastpos = pos+delimiter.length();
        }
        token = line.substr(lastpos, line.length()-lastpos);
        tokens.push_back(token);
        return std::move(tokens);
    }


    static std::string strip_path(std::string& str) {
        auto idx = str.find_last_of("/");
        return str.substr(idx, str.length());
    }

    template <typename Sequence, typename Pred>
    static Sequence& trim_end(Sequence& seq, Pred pred) {
        auto last = std::find_if_not(seq.rbegin(),
                                     seq.rend(),
                                     pred);
        seq.erase(last.base(), seq.end());
        return seq;
    }

    static std::string& trim_end(std::string& str, const char tr) {
        return trim_end(str, [&tr](const char c){ return c == tr; });
    }

    static std::string split(std::string str, char c) {
        return str.substr(0, str.find_first_of(c));
    }

    static std::string join(std::vector<std::string> list, std::string sep="") {
        std::string result = "";
        if (list.empty()) return result;

        result = list[0];
        for (int i = 1; i < list.size(); i++) {
            result += sep + list[i];
        }
        return result;
    }
    
    static void split(std::vector<std::string> &tokens, const std::string& s, std::string delimiter)
    {
        tokens.clear();

        std::string token;
        
        std::string line = s;
        int lastpos = 0;
        int pos = 0;
        
        while ((pos = line.find(delimiter, lastpos)) != std::string::npos) {
            token = line.substr(lastpos, pos-lastpos);
            tokens.push_back(token);
            lastpos = pos+delimiter.length();
        }
        token = line.substr(lastpos, line.length()-lastpos);
        tokens.push_back(token);
    }
    
    template <class T>
    static bool vectorContains(std::vector<T> v, T element) {
        return (v.find(element) == v.end());
    }

    static std::string FormatSeconds(size_t s) {
        std::string result = "";

        uint hours = s/3600;
        s -= hours * 3600;
        uint mins = s/60;
        s -= mins * 60;
        uint secs = s;

        if (hours) result += std::to_string(hours) + "h ";
        if (mins) result += std::to_string(mins) + "m ";
        if (secs) result += std::to_string(secs) + "s";
        return result;
    }


    static std::string FormatMilliseconds(size_t ms) {
        std::string result = "";

        uint hours = ms/3600000;
        ms -= hours * 3600000;
        uint mins = ms/60000;
        ms -= mins * 60000;
        uint secs = ms/1000;
        ms -= secs * 1000;

        if (hours) result += std::to_string(hours) + "h ";
        if (mins) result += std::to_string(mins) + "m ";
        if (secs) result += std::to_string(secs) + "s ";
        if (ms) result += std::to_string(ms) + "ms";
        return result;
    }
    
    static inline bool exists (const std::string& name) {
        if (FILE *file = fopen(name.c_str(), "r")) {
            fclose(file);
            return true;
        } else {
            return false;
        }
    }
    
    static std::string to_string_precision(double l, int precision) {
        std::ostringstream strs;
        strs << std::fixed << std::setprecision(precision)  << l;
        return strs.str();
    }
    
    template <class T, class U>
    static bool mapContains(std::unordered_map<T, U> v, T element) {
        return (v.find(element) == v.end());
    }
    
    static uint64_t swapEndianess(uint64_t &value, uint32_t bytes) {
        uint8_t swap[8];
        memset(swap, 0 ,8);
        for (int i = 7; i > (7-bytes); i--) {
            swap[7 - i] = ((uint8_t *) &value)[i];
        }
        return *((uint64_t*)swap);
    }
    
    static void swapEndianess(uint8_t* value, int32_t bytes) {
        uint8_t cache;
        for (int i = 7; i > (7-bytes) && (i >= 4); i--) {
            cache = value[i];
            value[i] = value[7-i];
            value[7-i] = cache;
        }
    }
    
    static std::string stripExtension(std::string file) {
        int index = file.rfind('.');
        int slash_i = file.rfind('/');
        if (index > 0 && slash_i > 0 && index > slash_i) {
            return file.substr(slash_i+1, index-slash_i-1);
        }
        return file;
    }

    static long GetFileSize(std::string filename)
    {
        struct stat stat_buf;
        int rc = stat(filename.c_str(), &stat_buf);
        return rc == 0 ? stat_buf.st_size : -1;
    }



//// reference implementation
//    int FastHammondRef(const char *s, const char *t, int length)
//    {
//        int result = 0;
//        int i;
//
//        for (i = 0; i < length; ++i)
//        {
//            if (s[i] == t[i])
//                result++;
//        }
//        return result;
//    }
//
//// optimised implementation
//    int FastHammond(const char *s, const char *t, int length)
//    {
//        int result = 0;
//        int i;
//
//        __m128i vsum = _mm_set1_epi32(0);
//        for (i = 0; i < length - 15; i += 16)
//        {
//            __m128i vs, vt, v, vh, vl, vtemp;
//
//            vs = _mm_loadu_si128((__m128i *)&s[i]); // load 16 chars from input
//            vt = _mm_loadu_si128((__m128i *)&t[i]);
//            v = _mm_cmpeq_epi8(vs, vt);             // compare
//            vh = _mm_unpackhi_epi8(v, v);           // unpack compare result into 2 x 8 x 16 bit vectors
//            vl = _mm_unpacklo_epi8(v, v);
//            vtemp = _mm_madd_epi16(vh, vh);         // accumulate 16 bit vectors into 4 x 32 bit partial sums
//            vsum = _mm_add_epi32(vsum, vtemp);
//            vtemp = _mm_madd_epi16(vl, vl);
//            vsum = _mm_add_epi32(vsum, vtemp);
//        }
//
//        // get sum of 4 x 32 bit partial sums
//        vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 8));
//        vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 4));
//        result = _mm_cvtsi128_si32(vsum);
//
//        // handle any residual bytes ( < 16)
//        if (i < length)
//        {
//            result += FastHammondRef(&s[i], &t[i], length - i);
//        }
//
//        return result;
//    }
}