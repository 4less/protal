//
// Created by fritsche on 05/09/23.
//

#pragma once

#include <vector>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include "Utils.h"

namespace protal {
    using ClosestNeighborMatrix = std::vector<std::vector<uint32_t>>;
    class TreeFilter {
    private:
        ClosestNeighborMatrix m_close_neighbors;

        static void Load(std::string& filepath) {
            ClosestNeighborMatrix matrix;
            std::vector<uint32_t> row(3, 0);
            try {
                std::ifstream is(filepath);
                std::string line;
                std::vector<std::string> tokens;
                if (is.is_open()) {
                    // File is successfully opened
                    // Perform operations on the file
                    while (std::getline(is, line)) {
                        Utils::split(tokens, line, "\t");
                        row.clear();
                        if (line.size() != 3) {
                            throw std::runtime_error(filepath + " has more or less than 3 columns: \nCurrent line:\n" + line);
                        }
                        auto row_index = std::stoul(tokens[0]);
                        for (auto i = 1; i < tokens.size(); i++) {
                            if (row_index >= matrix.size()) matrix.resize(row_index+1);
                            row.emplace_back(std::stoul(tokens[i]));
                        }
                        matrix[row_index] = row;
                    }

                    is.close();
                } else {
                    // File failed to open
                    throw std::runtime_error("Failed to open file: " + filepath);
                }
            } catch (const std::exception& e) {
                // Exception occurred during file opening
                std::cerr << "Error: " << e.what() << std::endl;
            }
        }

    public:
        TreeFilter(std::string file) {
            Load(file);
        }
    };
}