//
// Created by fritsche on 17/08/22.
//

#pragma once

#include <optional>
#include <string>
#include <vector>
#include "FastxReader.h"

namespace protal {
    enum DebugLevel {
        DEBUG_NONE,
        DEBUG_VERBOSE,
        DEBUG_EXTRAVERBOSE,
        DEBUG_BENCHMARK
    };

    struct NoBenchmark {
        void Print() {
            std::cout << "Dummy print" << std::endl;
        }
    };

    using TaxId = uint32_t;
    using GeneId = uint32_t;
    using GenePos = uint32_t;

    // Pair k-mer, readpos
    using KmerElement = std::pair<size_t, size_t>;
    using KmerList = std::vector<KmerElement>;

    struct LookupResult;
    using LookupList = std::vector<LookupResult>;
    using SeedList = LookupList;

    struct ChainAlignmentAnchor;
    using CAlignmentAnchor = ChainAlignmentAnchor;
    using AlignmentAnchorList = std::vector<CAlignmentAnchor>;

    struct AlignmentResult;
    using AlignmentResultList = std::vector<AlignmentResult>;
    using PairedAlignment = std::pair<AlignmentResult, AlignmentResult>;
    using PairedAlignmentResultList = std::vector<PairedAlignment>;

    struct InternalReadAlignment;
    using OptIRA = std::optional<InternalReadAlignment>;
    using OptIRAPair = std::pair<OptIRA, OptIRA>;

    using InternalReadAlignmentList = std::vector<InternalReadAlignment>;
    using InternalNonUniqueReadAlignmentList = std::vector<InternalReadAlignment>;
    using InternalNonUniqueReadOptIRAPairList = std::vector<OptIRAPair>;
    using InternalNonUniqueReadOptIRAPairListList = std::vector<InternalNonUniqueReadOptIRAPairList>;
    using InternalNonUniqueReadAlignmentListList = std::vector<InternalNonUniqueReadAlignmentList>;


    struct Sample;
    using SampleList = std::vector<Sample>;

    template <typename T>
    concept KmerIteratorConcept =
    requires(T t, size_t s, char* cs, KmerList klist) {
        { t(s) } -> std::same_as<bool>;
        { t.SetSequence(cs, s) } -> std::same_as<void>;
        { t.GetPos() } -> std::same_as<int64_t>;
        { t.KmerList(klist) } -> std::same_as<void>;
    };

    template <typename T>
    concept KmerHandlerConcept =
    requires(T t, std::string_view sv, std::string s, KmerList klist, size_t k) {
        { t(s, klist) } -> std::same_as<void>;
        { t.SetSequence(sv) } -> std::same_as<void>;
        { t.GetPos() } -> std::same_as<int64_t>;
        { t.Next(k) } -> std::same_as<bool>;
    };

    template <typename T>
    concept KmerStatisticsConcept =
    requires(T t) {
        { t.TotalKmers() } -> std::same_as<size_t>;
    };

    template <typename T>
    concept MinimizerStatisticsConcept =
    requires(T t) {
        { t.TotalMinimizers() } -> std::same_as<size_t>;
    };


    template <typename T>
    concept ContextFreeMinimizer =
    requires(T t, size_t s) {
        { t(s) } -> std::same_as<bool>;
    };


    template <typename T>
    concept KmerLookupConcept =
    requires(T t, LookupList l, std::istream is, uint32_t i, size_t k) {
        { t.Get(l, k, i) } -> std::same_as<void>;
//        { t.Load(is) } -> std::same_as<void>;
    };

    /*
     * Takes a list of k-mers and (k) and dumps a set of anchors into
     * AlignmentAnchorList
     * Check the using-type of AlignmentAnchorList to know what is expected
     */
    template <typename T>
    concept AnchorFinderConcept =
    requires(T t, KmerList k, SeedList s, AlignmentAnchorList a, std::string str) {
        { t(k, s, a, str) } -> std::same_as<void>;
    };

    /*
     * AlignmentHandlerConcept defines an interface that requires an
     * AlignmentAnchorList as potential targets and the query sequence as
     * input and gives an AlignmentResultList back.
     */
    template <typename T>
    concept AlignmentHandlerConcept =
    requires(T t, AlignmentAnchorList a, AlignmentResultList r, std::string s, size_t top) {
        { t(a, r, s, top, s) } -> std::same_as<void>;
    };

    /*
     * AlignmentOutputConcept defines an interface that requires an
     * AlignmentResultList as an input and performs some sort of output
     * which is up to the specific implementation of the concept.
     */
    template <typename T>
    concept AlignmentOutputConcept =
    requires(T t, AlignmentResultList r) {
        { t(r) } -> std::same_as<bool>;
    };



    template <typename T>
    concept KmerProcessorConcept =
    requires(T t, KmerList l, FastxRecord record) {
        { t(l, record) } -> std::same_as<void>;
    };
}