//
// Created by fritsche on 18/10/22.
//

#pragma once

#include <vector>
#include <string>
#include <stdexcept>

namespace protal {
    class LineSplitter {
        std::vector<std::string> m_tokens;
        std::string m_delimiter = "\t";

    public:
        static bool IsDelim(std::string &str, std::string& delim, int pos) {
            for (int i = 0; i < delim.length() && pos+i < str.length(); i++) {
                if (str[pos+i] != delim[i]) return false;
            }
            return true;
        }

        static void Split(std::string& line, std::string& delimiter, std::vector<std::string> &tokens) {
            tokens.clear();
            if (line.empty()) return;

            size_t delim_size = delimiter.length();

            size_t start = 0;

            for (int i = 0; i <= line.length() - delim_size; i++) {
                if (IsDelim(line, delimiter, i)) {
                    tokens.emplace_back(std::string(line.c_str() + start, i  - start));
                    i += delim_size;
                    start = i;
                }
            }
            tokens.emplace_back(std::string(line.c_str() + start, line.length() - start));
        }

        std::vector<std::string>& Tokens() {
            return m_tokens;
        }

        static void Split(std::string& line, std::string&& delimiter, std::vector<std::string> &tokens) {
            tokens.clear();

            size_t delim_size = delimiter.length();

            size_t start = 0;

            if (line.length() < delimiter.length()) return;

            for (int i = 0; i <= line.length() - delim_size; i++) {
                if (IsDelim(line, delimiter, i)) {
                    tokens.emplace_back(std::string(line.c_str() + start, i  - start));
                    i += delim_size;
                    start = i;
                }
            }
            tokens.emplace_back(std::string(line.c_str() + start, line.length() - start));
        }


        void Split(std::string line) {
            Split(line, m_delimiter, m_tokens);
        }

        std::string operator[] (size_t index) {
            if (index >= m_tokens.size()) {
                throw std::out_of_range("Line splitter: trying to access index " + std::to_string(index) + " of split line. Number of tokens is " + std::to_string(m_tokens.size()));
            }
            return m_tokens[index];
        }
    };
}