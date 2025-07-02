#pragma once
#include <string>
#include <cctype>
#include <utility>
#include <tuple> 


class CigarIterator {
public:
    // value_type now holds: cpos, qpos, rpos, op
    using value_type = std::tuple<size_t, size_t, size_t, char>;

    class iterator {
    public:
        iterator(const std::string* cigar, size_t cpos, size_t qpos, size_t rpos, size_t cigar_pos, size_t op_len, char op_char, size_t op_idx, bool end)
            : m_cigar(cigar), m_cpos(cpos), m_qpos(qpos), m_rpos(rpos), m_cigar_pos(cigar_pos), m_op_len(op_len), m_op_char(op_char), m_op_idx(op_idx), m_end(end)
        {}

        value_type operator*() const {
            return {m_cpos, m_qpos, m_rpos, m_op_char};
        }

        iterator& operator++() {
            if (m_end) return *this;

            // Advance cpos always
            ++m_cpos;

            // Advance qpos and rpos depending on op
            if (m_op_char == 'D') {
                // Deletion: only reference advances
                ++m_rpos;
            } else if (m_op_char == 'I') {
                // Insertion: only query advances
                ++m_qpos;
            } else {
                // Match/Mismatch/Other: both advance
                ++m_qpos;
                ++m_rpos;
            }

            ++m_op_idx;
            if (m_op_idx >= m_op_len) {
                advance_cigar();
            }
            return *this;
        }

        bool operator==(const iterator& other) const {
            return m_end == other.m_end;
        }
        bool operator!=(const iterator& other) const {
            return !(*this == other);
        }

    private:
        void advance_cigar() {
            m_op_len = 0;
            m_op_char = 0;
            m_op_idx = 0;
            while (m_cigar_pos < m_cigar->size() && std::isdigit((*m_cigar)[m_cigar_pos])) {
                m_op_len = m_op_len * 10 + ((*m_cigar)[m_cigar_pos] - '0');
                ++m_cigar_pos;
            }
            if (m_cigar_pos < m_cigar->size() && std::isalpha((*m_cigar)[m_cigar_pos])) {
                m_op_char = (*m_cigar)[m_cigar_pos];
                ++m_cigar_pos;
            } else {
                m_end = true;
                return;
            }
            if (m_op_len == 0) {
                m_end = true;
            }
            m_op_idx = 0;
        }

        const std::string* m_cigar;
        size_t m_cpos;
        size_t m_qpos;
        size_t m_rpos;
        size_t m_cigar_pos;
        size_t m_op_len;
        char m_op_char;
        size_t m_op_idx;
        bool m_end;
        friend class CigarIterator;
    };

    CigarIterator(const std::string& cigar)
        : m_cigar(cigar)
    {}

    iterator begin() const {
        iterator it(&m_cigar, 0, 0, 0, 0, 0, 0, 0, false);
        it.advance_cigar();
        return it;
    }

    iterator end() const {
        return iterator(&m_cigar, 0, 0, 0, m_cigar.size(), 0, 0, 0, true);
    }

private:
    const std::string& m_cigar;
};