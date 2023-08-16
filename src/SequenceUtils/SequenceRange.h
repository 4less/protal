//
// Created by fritsche on 21/02/23.
//

#ifndef PROTAL_SEQUENCERANGE_H
#define PROTAL_SEQUENCERANGE_H


#include <cstdio>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>

using CoverageVec = std::vector<uint16_t>;

struct ReadInfo {
    size_t read_id;
    uint32_t start;
    uint32_t length;
    bool forward;

    std::string ToString() const {
        return "read_id: " + std::to_string(read_id) + " [" + std::to_string(start) + ", " + std::to_string(start+length) + "]:" + std::to_string(length) + " forward: " + std::to_string(forward);
    }

    bool Contains(size_t pos) const {
        return pos >= start && pos < (start + length);
    }
};

class SequenceRange {
    size_t m_start = 0;
    size_t m_end = 0;

    // Collect read evidence
    mutable std::vector<ReadInfo> m_read_info;

    bool Overlap(size_t start, size_t end) const;

    bool Overlap(SequenceRange const& range) const {
        return Overlap(range.m_start, range.m_end);
    }

    void Intersect(SequenceRange const& range) {
        assert(*this == range);
        m_start = std::max(m_start, range.m_start);
        m_end = std::min(m_end, range.m_end);

        auto end = std::remove_if(m_read_info.begin(),
                                  m_read_info.end(),
                                  [this](ReadInfo const& read_info) {
                                      return !this->Contains(read_info.start);    // remove odd numbers
                                  });

        m_read_info.erase(end, m_read_info.end());

    }

public:
    SequenceRange(size_t start, size_t end) : m_start(start), m_end(end) {};

    void AddReadInfo(ReadInfo const&& read_info) {
        m_read_info.emplace_back(read_info);
    }

    void AddReadInfo(ReadInfo const& read_info) {
        m_read_info.emplace_back(read_info);
    }

    const std::vector<ReadInfo> GetReadInfo() const {
        return m_read_info;
    }

    void SortReadInfo() const {
        std::sort(m_read_info.begin(), m_read_info.end(), [](ReadInfo const& a, ReadInfo const& b) {
            return a.start < b.start;
        });
    }

    CoverageVec CoverageVector() const {
        CoverageVec coverage;

        coverage.resize(m_end - m_start, 0);
        for (auto& read : m_read_info) {
            for (auto i = 0; i < read.length; i++) {
                // vector starts where m_start is, so m_start is the offset.
                coverage[i+(read.start-m_start)]++;
            }
        }

        return coverage;
    }

    bool operator < (const SequenceRange& range) const {
        return (m_end < range.m_start);
    }
    bool operator > (const SequenceRange& range) const {
        return (m_start > range.m_end);
    }
    bool operator == (const SequenceRange& range) const {
        return Overlap(range);
    }
    bool operator != (const SequenceRange& range) const {
        return !(*this == range);
    }
    bool operator <= (const SequenceRange& range) const {
        return *this < range || *this == range;
    }
    bool operator >= (const SequenceRange& range) const {
        return *this > range || *this == range;
    }

    size_t Length() const {
        return m_end - m_start;
    }

    std::string ToString() const {
        return "[" + std::to_string(m_start) + ',' + std::to_string(m_end) + "](" + std::to_string(m_read_info.size()) + ")";
    }

    std::string ToVerboseString() const {
        std::string output = "[" + std::to_string(m_start) + ',' + std::to_string(m_end) + "] {\n";
        for (auto& read_info : m_read_info) {
            output += read_info.ToString() + "\n";
        }
        return output;
    }

    size_t Start() const {
        return m_start;
    }

    size_t End() const {
        return m_end;
    }

    friend class SequenceRangeHandler;

    bool Contains(size_t pos) const {
        return pos >= m_start && pos < m_end;
    }

    bool Contains(size_t start, size_t end) const {
        return start >= m_start && end < m_end;
    }

    size_t Coverage(size_t pos) const {
        auto count = std::count_if(m_read_info.begin(), m_read_info.end(), [&pos](ReadInfo const& read_info) {
            return read_info.Contains(pos);
        });
//            std::cout << "count: " << count << " total: " << m_read_info.size() << std::endl;
        return count;
    }

    void Union(SequenceRange const& range) {
        assert(*this == range);
        m_start = std::min(m_start, range.m_start);
        m_end = std::max(m_end, range.m_end);

        // Copy elements from src to dest.
        // (Could move, but don't want to leave old object
        // in undefined state)
        m_read_info.insert(
                m_read_info.end(),
                range.m_read_info.begin(),
                range.m_read_info.end()
        );
    }
};


#endif //PROTAL_SEQUENCERANGE_H
