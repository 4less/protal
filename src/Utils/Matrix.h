//
// Created by fritsche on 27/07/23.
//

#ifndef PROTAL_MATRIX_H
#define PROTAL_MATRIX_H

#include <cstdio>
#include <vector>
#include <robin_map.h>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <math.h>

using Key = std::string;
using Value = size_t;

//template<typename K, typename V>
//using Map = std::unordered_map<K, V>;

template<typename K, typename V>
using Map = tsl::robin_map<K, V>;

template<typename T>
class Matrix {
    template<typename K, typename V>
    using Map = std::unordered_map<K, V>;
    using NamesDict = Map<std::string, size_t>;
    using NamesList = std::vector<std::string>;
    using RowType = std::vector<T>;
    using MatrixType = std::vector<RowType>;

    NamesDict m_row_names_to_index;
    NamesDict m_col_names_to_index;
    NamesList m_row_names;
    NamesList m_col_names;
    MatrixType m_matrix;

    bool any_set = false;

public:
    Matrix(size_t row_size, size_t col_size, T default_value) :
            m_matrix(std::vector<RowType>(row_size, RowType(col_size, default_value))) {};
    Matrix() {};

    void SetNames(NamesList& from, NamesList& to_list, NamesDict& to_dict) {
        to_list = from;
        to_dict.clear();
        for (auto i = 0; i < to_list.size(); i++) {
            to_dict.insert({to_list[i], i});
        }
    }

    void AddName(std::string new_name) {
        m_row_names_to_index.insert( { new_name, m_row_names.size() });
        m_col_names_to_index.insert( { new_name, m_col_names.size() });
        m_row_names.emplace_back(new_name);
        m_col_names.emplace_back(new_name);

        m_matrix.resize(m_matrix.size() + 1, RowType(m_matrix.size(), NAN));
        std::for_each(m_matrix.begin(), m_matrix.end(), [](RowType& row) {
            row.resize(row.size() + 1, NAN);
        });
    }

    void SetColNames(NamesList& col_names) {
        SetNames(col_names, m_col_names, m_col_names_to_index);
    }
    void SetRowNames(NamesList& row_names) {
        SetNames(row_names, m_row_names, m_row_names_to_index);
    }
    void SetColNames(NamesList&& col_names) {
        SetNames(col_names, m_col_names, m_col_names_to_index);
    }
    void SetRowNames(NamesList&& row_names) {
        SetNames(row_names, m_row_names, m_row_names_to_index);
    }

    std::string GetName(NamesList const& names, size_t index) const {
        if (names.empty()) {
            return std::to_string(index);
        } else {
            return names[index];
        }
    }

    bool HasName(std::string& name) {
        return m_row_names_to_index.contains(name);
    }

    bool HasIdx(uint32_t index) {
        return index < m_row_names.size();
    }

    bool AnySet() const {
        return any_set;
    }

    void SetValue(size_t row_idx, size_t col_idx, T value) {
        if (row_idx >= m_matrix.size()) {
            std::cerr << "row_idx out of bounds." << row_idx << std::endl;
            exit(9);
        }
        if (m_matrix.empty() || col_idx >= m_matrix.front().size()) {
            std::cerr << "col_idx out of bounds." << col_idx << std::endl;
            exit(9);
        }
        m_matrix[row_idx][col_idx] = value;
        any_set = true;
    }

    void SetValue(std::string& row_name, std::string& col_name, T value) {
        auto row_idx_find = m_row_names_to_index.find(row_name);
        auto col_idx_find = m_col_names_to_index.find(col_name);

        if (row_idx_find == m_row_names_to_index.end()) {
            std::cerr << "row has no name " << row_name << std::endl;
            exit(9);
        }
        if (col_idx_find == m_col_names_to_index.end()) {
            std::cerr << "col has no name " << col_name << std::endl;
            exit(9);
        }
        SetValue(row_idx_find->second, col_idx_find->second, value);
    }

    void SetValue(size_t row_idx, size_t col_idx, T value, bool lower_triangle) {
        if (row_idx >= m_matrix.size()) {
            std::cerr << "row_idx out of bounds." << row_idx << std::endl;
            exit(9);
        }
        if (m_matrix.empty() || col_idx >= m_matrix.front().size()) {
            std::cerr << "col_idx out of bounds." << col_idx << std::endl;
            exit(9);
        }

        if (lower_triangle) {
            if (row_idx < col_idx) std::swap(row_idx, col_idx);
            m_matrix[row_idx][col_idx] = value;
        } else {
            if (row_idx > col_idx) std::swap(row_idx, col_idx);
            m_matrix[row_idx][col_idx] = value;
        }
    }

    void SetValue(std::string& row_name, std::string& col_name, T value, bool lower_triangle) {
        auto row_idx_find = m_row_names_to_index.find(row_name);
        auto col_idx_find = m_col_names_to_index.find(col_name);

        if (row_idx_find == m_row_names_to_index.end()) {
            std::cerr << "row has no name " << row_name << std::endl;
            exit(9);
        }
        if (col_idx_find == m_col_names_to_index.end()) {
            std::cerr << "col has no name " << col_name << std::endl;
            exit(9);
        }
        SetValue(row_idx_find->second, col_idx_find->second, value, lower_triangle);
        any_set = true;
    }


    void PrintMatrix(std::ostream& os=std::cout, std::string const& sep="\t", size_t precision=4) const {
        os << std::fixed;
        for (auto col_idx = 0; col_idx < m_col_names.size(); col_idx++) {
            os << sep << GetName(m_col_names, col_idx);
        }
        os << std::endl;
        for (auto row_idx = 0; row_idx < m_row_names.size(); row_idx++) {
            os << GetName(m_row_names, row_idx);

            for (auto col_idx = 0; col_idx < m_col_names.size(); col_idx++) {
                os << sep << std::setprecision(precision) << m_matrix[row_idx][col_idx];
            }
            os << std::endl;
        }
    };

};

#endif //PROTAL_MATRIX_H
