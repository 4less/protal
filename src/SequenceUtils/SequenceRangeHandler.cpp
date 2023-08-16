//
// Created by fritsche on 21/02/23.
//

#include <numeric>
#include "SequenceRangeHandler.h"


bool SequenceRange::Overlap(size_t start, size_t end) const {
    // MS###S#####ME Start between m_start and m_end or
    // MS##E######ME End between m_start and m_end
    // Start is inclusive, end is exclusive
    return start >= m_start && start <= m_end || end >= m_start && end <= m_end || m_start >= start && m_start <= end || m_end >= start && m_end <= end;
}



void SequenceRangeHandler::Add(size_t start, size_t end) {
    SequenceRange query_range(start, end);
    auto range_it = FindSequenceRange(start, end);

    if (range_it == m_ranges.end() || *range_it != query_range) {
        m_ranges.insert(range_it, query_range);

    } else if (*range_it == query_range) {
        range_it->Union(query_range);
        if (range_it+1 != m_ranges.end() && *(range_it+1) == *range_it) {
            auto del_it = range_it + 1;
            range_it->Union(*del_it);
            m_ranges.erase(del_it);
        }
    }
}

SequenceRangeHandler SequenceRangeHandler::Unify(const SequenceRangeHandler &handler) const {
    SequenceRangeHandler union_ranges;

    int this_i = 0;
    int other_i = 0;

    while (this_i < m_ranges.size() && other_i < handler.m_ranges.size()) {
        if (m_ranges[this_i] == handler.m_ranges[other_i]) {
            // create new copy
            auto union_range = m_ranges[this_i];
            // Intersect with other
            union_range.Union(handler.m_ranges[other_i]);
            // Push to new intersection handler
            union_ranges.m_ranges.emplace_back(union_range);
        }


        bool increment_this = m_ranges[this_i].m_end < handler.m_ranges[other_i].m_end;
        this_i += increment_this;
        other_i += !increment_this;

    }
    return union_ranges;
}

SequenceRangeList &SequenceRangeHandler::GetRanges() {
    return m_ranges;
}

const SequenceRangeList &SequenceRangeHandler::GetRanges() const {
    return m_ranges;
}


SequenceRangeHandler SequenceRangeHandler::Intersect(const SequenceRangeHandler &handler) const {
    SequenceRangeHandler intersection;

    int this_i = 0;
    int other_i = 0;

    while (this_i < m_ranges.size() && other_i < handler.m_ranges.size()) {
        if (m_ranges[this_i] == handler.m_ranges[other_i]) {
            // create new copy
            auto new_range = m_ranges[this_i];
            // Intersect with other
            new_range.Intersect(handler.m_ranges[other_i]);
            // Push to new intersection handler
            intersection.m_ranges.emplace_back(new_range);
        }


        bool increment_this = m_ranges[this_i].m_end < handler.m_ranges[other_i].m_end;
        this_i += increment_this;
        other_i += !increment_this;

    }
    return intersection;
}

//SequenceRangeHandler::SRIterator SequenceRangeHandler::FindSequenceRange(size_t start, size_t end) {
//    auto range_it = m_ranges.begin();
//    SequenceRange query_range(start, end);
//
//    // ranges are sorted in ascending order by their start.
//    // increment iterator until range is smaller
//    while (range_it != m_ranges.end() && *range_it < query_range) {
//        range_it++;
//    }
//    if (*range_it == query_range) return range_it;
//    return m_ranges.end();
//}

std::string SequenceRangeHandler::ToString() const {
    if (m_ranges.empty()) return "[]";
    std::string str = "[";
    str += m_ranges[0].ToString();
    for (int i = 1; i < m_ranges.size(); i++) {
        str += ", " + m_ranges[i].ToString();
    }
    str += "]";
    return str;
}

std::string SequenceRangeHandler::ToVerboseString() const {
    if (m_ranges.empty()) return "[]";
    std::string str = "";
    for (int i = 0; i < m_ranges.size(); i++) {
        str += m_ranges[i].ToVerboseString();
    }
    return str;
}


size_t SequenceRangeHandler::SequenceLength() const {
    return std::accumulate(m_ranges.begin(), m_ranges.end(), 0, [](int a, SequenceRange const& range){ return range.Length(); });
}

size_t SequenceRangeHandler::Size() const {
    return m_ranges.size();
}

SequenceRangeHandler::SRIterator SequenceRangeHandler::begin() noexcept {
    return m_ranges.begin();
}

SequenceRangeHandler::SRIterator SequenceRangeHandler::end() noexcept {
    return m_ranges.end();
}

SequenceRangeHandler::SRIterator
SequenceRangeHandler::GetOrAddSequenceRangeIterator(const size_t start, const size_t end) {
    auto range_it = m_ranges.begin();
    SequenceRange query_range(start, end);

    // ranges are sorted in ascending order by their start.
    // increment iterator until range is smaller
    while (range_it != m_ranges.end() && *range_it < query_range) {
        std::advance(range_it, 1);
    }
    if (range_it == m_ranges.end() || *range_it != query_range) {
        // Insert new range at end
        m_ranges.insert(range_it, query_range);
    } else if (*range_it == query_range) {
        range_it->Union(query_range);
        if (range_it+1 != m_ranges.end() && *(range_it+1) == *range_it) {
            auto del_it = range_it + 1;
            range_it->Union(*del_it);
            m_ranges.erase(del_it);
        }
    }
    return range_it;
}

SequenceRangeHandler::SRIterator
SequenceRangeHandler::FindSequenceRange(const size_t start, const size_t end) {
    auto range_it = m_ranges.begin();
    SequenceRange query_range(start, end);

    // ranges are sorted in ascending order by their start.
    // increment iterator until range is smaller
    while (range_it != m_ranges.end() && *range_it < query_range) {
//                std::cout << range_it->ToString() << " < " << query_range.ToString() << " " << (query_range > *range_it) << std::endl;
        range_it++;
    }
    return range_it;
}

const SequenceRange &SequenceRangeHandler::GetRange(size_t pos) const {
    auto find = std::find_if(m_ranges.begin(), m_ranges.end(), [&pos](SequenceRange const& range){
        return range.Contains(pos);
    });
    //TODO implement exception
    return *find;
}

SequenceRange &SequenceRangeHandler::GetRange(size_t pos) {
    auto find = std::find_if(m_ranges.begin(), m_ranges.end(), [&pos](SequenceRange const& range){
        return range.Contains(pos);
    });
    //TODO implement exception
    return const_cast<SequenceRange&>(*find);
}

bool SequenceRangeHandler::HasRange(size_t pos) const {
    auto find = std::find_if(m_ranges.begin(), m_ranges.end(), [&pos](SequenceRange const& range) {
        return range.Contains(pos);
    });
    //TODO implement exception
    return find != m_ranges.end();
}

void SequenceRangeHandler::Add(SequenceRange &range) {
    m_ranges.emplace_back(range);
}

void SequenceRangeHandler::Add(SequenceRange &&range) {
    m_ranges.emplace_back(range);
}

CoverageVec SequenceRangeHandler::CalculateCoverageVector() {
    for (auto &range : m_ranges) {
        if (range.m_start > m_cov.size()) {
            m_cov.resize(range.m_start, 0);
        }
        auto rcov = range.CoverageVector();
        m_cov.insert(m_cov.end(), rcov.begin(), rcov.end());
    }
    return m_cov;
}


CoverageVec SequenceRangeHandler::CalculateCoverageVector2() const {
    CoverageVec cov;
    for (auto &range : m_ranges) {
        if (range.m_start > cov.size()) {
            cov.resize(range.m_start, 0);
        }
        auto rcov = range.CoverageVector();
        cov.insert(cov.end(), rcov.begin(), rcov.end());
    }
    return cov;
}

size_t SequenceRangeHandler::CoveredPortion(uint16_t min_cov) {
    size_t count = 0;
    for (auto &range : m_ranges) {
        auto cov = range.CoverageVector();
        count += std::count_if(cov.begin(), cov.end(), [&min_cov](uint16_t const e) {
            return e >= min_cov;
        });
    }
    return count;
}


SequenceRangeHandler SequenceRangeHandler::FromCoverageVector(CoverageVec const& coverage_vector, size_t min_sequence_length) const {
    SequenceRangeHandler result{};

    auto start = -1;
    for (auto i = 0; i < coverage_vector.size(); i++) {
        auto cov = coverage_vector[i];
        if (cov == 0) {
            if (start != -1 && (i - start) > min_sequence_length) result.Add(start, i);
            start = -1;
        } else if (start == -1) start = i;
    }
    result.SetCoverageVector(coverage_vector);

    return result;
}


SequenceRangeHandler SequenceRangeHandler::Intersect(SequenceRangeHandler &handler, size_t min_coverage, size_t min_sequence_length) {
    auto cov1 = CalculateCoverageVector2();
    auto cov2 = handler.CalculateCoverageVector2();
    CoverageVec intersection = IntersectCoverageVectors(cov1, cov2, min_coverage);

    return FromCoverageVector(intersection, min_sequence_length);
}

void SequenceRangeHandler::SetCoverageVector(const CoverageVec &vector) {
    m_cov = vector;
}
