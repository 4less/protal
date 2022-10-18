//
// Created by fritsche on 18/10/22.
//

#pragma once

#include <vector>
#include <string>

namespace protal {
    class LineSplitter {
        std::vector<std::string> m_tokens;
        std::string m_delimiter = "\t";

    public:
        static bool IsDelim(std::string &str, std::string& delim, int pos) {
            for (int i = 0; i < delim.length(); i++) {
                if (str[pos+i] != delim[i]) return false;
            }
            return true;
        }

        static void Split(std::string& line, std::string& delimiter, std::vector<std::string> &tokens) {
            tokens.clear();

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

        void Split(std::string line) {
            Split(line, m_delimiter, m_tokens);
        }

    };
}