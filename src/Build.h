//
// Created by fritsche on 14/08/22.
//

#pragma once

#include "Options.h"
#include "SequenceUtils/SeqReader.h"
#include "SequenceUtils/KmerIterator.h"
#include "Statistics.h"
#include <iostream>
#include <fstream>
#include <omp.h>
#include "Constants.h"
#include "Hash/KmerPutter.h"
#include "robin_map.h"
#include "KmerUtils.h"

namespace protal::build {
    template<typename KmerHandler, typename KmerPutter, DebugLevel debug>
    static void Check(protal::Options const& options, KmerPutter& putter, KmerHandler& kmer_handler_global) {

        // Shared
        std::ifstream is(options.GetSequenceFilePath(), std::ios::in);
        size_t dummy = 0;
        int read_count = 0;

        // Set Thread Num
        omp_set_num_threads(1);

        size_t main_k = putter.GetMap().m_exact_k;
        size_t main_k_bits = putter.GetMap().m_main_bits;
        size_t flex_k = putter.GetMap().m_flex_k;
        size_t flex_k_bits = putter.GetMap().m_flex_k_bits;

        Statistics statistics;

        is = std::ifstream(options.GetSequenceFilePath(), std::ios::in);
        KmerLookupSM lookup_global(putter.GetMap());

#pragma omp parallel default(none) shared(std::cout, lookup_global, options, is, dummy, read_count, kmer_handler_global, statistics, putter, flex_k_bits, main_k_bits)
        {
            // Private variables
            FastxRecord record;

            // Extract variables from kmi_global
            KmerHandler kmer_handler(kmer_handler_global);
            SeqReader reader { is };
            size_t kmer = 0;
            Statistics thread_statistics;
            thread_statistics.thread_num = omp_get_thread_num();

            size_t taxid = 1;
            size_t geneid = 1;
            size_t genepos = 1;

            KmerList kmers;
            KmerLookupSM lookup(lookup_global);

            LookupList seeds;
            std::vector<Entry<20,20,20>*> max_sim_entries;
            uint32_t max_sim = 0;
            uint32_t best_possible_sim = putter.GetMap().m_flex_k;

            while (reader(record)) {
                kmer_handler.SetSequence(std::string_view(record.sequence));

                auto [taxonomic_id, gene_id] = KmerUtils::ExtractHeaderInformation(record.header);

                thread_statistics.reads++;

                // Retrieve kmers
                kmers.clear();
                seeds.clear();
                kmer_handler(std::string_view(record.sequence), kmers);

                for (auto pair : kmers) {
                    size_t pos = pair.second;
                    max_sim = 0;
                    // std::cout << pair.first << ", " << pair.second << std::endl;
                    max_sim_entries.clear();
                    lookup.GetFlex(pair.first, max_sim_entries, max_sim);
                    if (max_sim_entries.empty()) continue;
                }
            }


#pragma omp critical(statistics)
            statistics.Join(thread_statistics);
        }
    }

    template<typename KmerHandler, typename KmerPutter, DebugLevel debug>
    static Statistics Run(protal::Options const& options, KmerPutter& putter, KmerHandler& kmer_handler_global) {

        // Shared
        std::ifstream is(options.GetSequenceFilePath(), std::ios::in);
        size_t dummy = 0;
        int read_count = 0;

        // Set Thread Num
        omp_set_num_threads(options.GetThreads());

        Statistics statistics;

        size_t main_k = putter.GetMap().m_exact_k;
        size_t main_k_bits = putter.GetMap().m_main_bits;
        size_t flex_k = putter.GetMap().m_flex_k;
        size_t flex_k_bits = putter.GetMap().m_flex_k_bits;


        // TODO Fix multithreaded database building.
        if (options.GetThreads() > 1) {
            std::cerr << "Building db with multiple threads is currently broken" << std::endl;
            omp_set_num_threads(1);
        }

        std::cout << "Run Build" << std::endl;
        if constexpr(protal::HasFirstPut<KmerPutter>) {
#pragma omp parallel default(none) shared(std::cout, options, is, dummy, read_count, kmer_handler_global, statistics, putter, main_k_bits, flex_k_bits)
            {
                // Private variables
                FastxRecord record;

                // Extract variables from kmi_global
                KmerHandler kmer_handler(kmer_handler_global);
                SeqReader reader{ is };
                Statistics thread_statistics;
                thread_statistics.thread_num = omp_get_thread_num();

                KmerList kmers;

                std::cout << "Build: iterate records" << std::endl;

                while (reader(record)) {
                    thread_statistics.reads++;

                    // Retrieve kmers
                    kmers.clear();
                    kmer_handler(std::string_view(record.sequence), kmers);
                    for (auto pair : kmers) {
                        putter.FirstPut(pair.first);
                    }

                    if constexpr(KmerStatisticsConcept<KmerHandler>) {
                        thread_statistics.kmers_total += kmer_handler.TotalKmers();
                    }

                    if constexpr(KmerStatisticsConcept<KmerHandler>) {
                        thread_statistics.kmers_accepted += kmers.size();
                    }

                    if constexpr(debug == DEBUG_VERBOSE) {
                        thread_statistics.WriteStats(std::cout);
                    }
                    if constexpr(debug == DEBUG_EXTRAVERBOSE) {

                    }
                }

#pragma omp critical(statistics)
                statistics.Join(thread_statistics);
                std::cout << "minimizers: " << thread_statistics.kmers_accepted << std::endl;
            }

            // E.g. if k-mers are counted before they are inserted, call Initialize for put
            // To calculate the bucket sizes
            putter.InitializeForPut();
        }

        // Make sure input stream is reset to start
        is.clear();                 // clear fail and eof bits
        is.seekg(0, std::ios::beg); // back to the start!

#pragma omp parallel default(none) shared(std::cout, options, is, dummy, read_count, kmer_handler_global, statistics, putter, flex_k_bits, main_k_bits)
    {
        // Private variables
        FastxRecord record;

        // Extract variables from kmi_global
        KmerHandler kmer_handler(kmer_handler_global);
        SeqReader reader { is };
        size_t kmer = 0;
        Statistics thread_statistics;
        thread_statistics.thread_num = omp_get_thread_num();

        size_t taxid = 1;
        size_t geneid = 1;
        size_t genepos = 1;

        KmerList kmers;

        std::cout << "After first put" << std::endl;
        while (reader(record)) {
            kmer_handler.SetSequence(std::string_view(record.sequence));

            auto [taxonomic_id, gene_id] = KmerUtils::ExtractHeaderInformation(record.header);

            thread_statistics.reads++;

            // Retrieve kmers
            kmers.clear();
            kmer_handler(std::string_view(record.sequence), kmers);
            for (auto pair : kmers) {
                size_t pos = pair.second;
                putter.Put(pair.first, taxonomic_id, gene_id, pos);
            }

            if constexpr(KmerStatisticsConcept<KmerHandler>) {
                thread_statistics.kmers_total += kmer_handler.TotalKmers();
            }
            if constexpr(KmerStatisticsConcept<KmerHandler>) {
                thread_statistics.kmers_accepted += kmer_handler.TotalMinimizers();
            }

            if constexpr(debug == DEBUG_VERBOSE) {
                thread_statistics.WriteStats(std::cout);
            }
            if constexpr(debug == DEBUG_EXTRAVERBOSE) {

            }
        }


#pragma omp critical(statistics)
        statistics.Join(thread_statistics);
    }

        if (options.GetFullSequenceFilePath().empty()) {
            std::ofstream index_ostream(options.GetIndexFile(), std::ios::binary);
            putter.Save(index_ostream);
            index_ostream.close();
            return statistics;
        }

        std::cout << "Check Uniqueness: " << options.GetFullSequenceFilePath() << std::endl;

        omp_set_num_threads(options.GetThreads());
        is = std::ifstream(options.GetFullSequenceFilePath(), std::ios::in);
        KmerLookupSM lookup_global(putter.GetMap());

#pragma omp parallel default(none) shared(std::cout, lookup_global, options, is, dummy, read_count, kmer_handler_global, statistics, putter, flex_k_bits, main_k_bits)
    {
        // Private variables
        FastxRecord record;

        // Extract variables from kmi_global
        KmerHandler kmer_handler(kmer_handler_global);
        SeqReader reader { is };
        size_t kmer = 0;
        Statistics thread_statistics;
        thread_statistics.thread_num = omp_get_thread_num();

        size_t taxid = 1;
        size_t geneid = 1;
        size_t genepos = 1;

        KmerList kmers;
        KmerLookupSM lookup(lookup_global);

        LookupList seeds;
        std::vector<Entry<20,20,20>*> max_sim_entries;
        uint32_t max_sim = 0;
        uint32_t best_possible_sim = putter.GetMap().m_flex_k;

        while (reader(record)) {
            kmer_handler.SetSequence(std::string_view(record.sequence));

            auto [taxonomic_id, gene_id] = KmerUtils::ExtractHeaderInformation(record.header);

            thread_statistics.reads++;

            // Retrieve kmers
            kmers.clear();
            seeds.clear();
            kmer_handler(std::string_view(record.sequence), kmers);

            for (auto pair : kmers) {
                size_t pos = pair.second;
                max_sim = 0;
                // std::cout << pair.first << ", " << pair.second << std::endl;
                max_sim_entries.clear();
                lookup.GetFlex(pair.first, max_sim_entries, max_sim);
                if (max_sim_entries.empty()) continue;

                max_sim_entries.front()->Get(taxid, geneid, genepos);

                if (max_sim != best_possible_sim) {
                    continue;
                }
                if (max_sim_entries.size() == 1 && taxonomic_id == taxid) {
                    continue;
                }

                for (auto& entry : max_sim_entries) {
                    // std::cout << "SetNonUnique " << taxonomic_id << " != " << taxid << " entries: " << max_sim_entries.size() << " Isunique? " << entry->IsFlagUnique();
#pragma omp critical(SetNonUnique)
                    entry->SetFlagNonUnique();
                    // std::cout << " -> " << entry->IsFlagUnique() << std::endl;
                }
            }
        }


    #pragma omp critical(statistics)
        statistics.Join(thread_statistics);
    }

        std::cout << "Save unique kmer info: \n" << options.GetUniqueKmersFile() << std::endl;
        std::ofstream os(options.GetUniqueKmersFile());
        putter.GetMap().CountUniqueKmers(os, true, true);
        os.close();

        std::ofstream index_ostream(options.GetIndexFile(), std::ios::binary);
        putter.Save(index_ostream);
        index_ostream.close();


        return statistics;
    }


}