//
// Created by fritsche on 08/03/2021.
//

#pragma once

#include <cstdint>
#include <cstddef>
#include <iostream>
#include <err.h>
#include <sysexits.h>
#include <unordered_set>
#include <fstream>
#include "KmerUtils.h"
#include "compact_map_utils.h"
#include "ProgressBar.h"
#include <tgmath.h>

#include "robin_map.h"
#include "robin_set.h"
#include "sparse_map.h"
#include "sparse_set.h"

#define KEY_BITS 22
#define VALUE_BITS 42
#define MAX_VALUE ((1 << VALUE_BITS) - 1)

using Cell_T = uint64_t;

struct Cell {
    Cell_T data = 0llu;

public:
    inline bool empty(const size_t &value_bits) {
        return value(value_bits) == 0;
    }

    inline bool empty() {
        return data == 0;
    }

    inline uint64_t key(const size_t &value_bits) {
        return data >> value_bits;
    }

    inline uint64_t value(const size_t &value_bits) {
        return (data & ((1llu << value_bits) - 1));
    }

//    inline void paste(size_t &key, size_t &value, const size_t &value_bits) {
//        data = key << value_bits;
//        data |= value;
//    }

    inline void paste(size_t key, size_t value, const size_t &value_bits) {
        data = key << value_bits;
        data |= value;
    }

    inline void pasteSafe(size_t &key, size_t &value, const size_t &value_bits, const size_t &key_bits) {
        if (key_bits + value_bits != 64)
            errx(EX_SOFTWARE, "key size %llu and value size %llu must sum up to 64", (unsigned long long) key_bits, (unsigned long long) value_bits);
        if (key > ((1 << key_bits) - 1))
            errx(EX_SOFTWARE, "key size %llu is larger than max key size %llu", (unsigned long long) key, (unsigned long long) key_bits);
        if (value > ((1 << value_bits) - 1))
            errx(EX_SOFTWARE, "value size %llu is larger than max value size %llu", (unsigned long long) value, (unsigned long long) value_bits);
        paste(key, value, value_bits);
    }
};

class IndexedMap {
    size_t capacity_;
    size_t size_;

    Cell* map_ = nullptr;

    // debug
    size_t* psl_ = nullptr;

    IndexedMap(size_t capacity, size_t size, size_t key_bits,
               size_t value_bits, double load_factor, size_t offset_bits,
               size_t offset_size, uint64_t *offset, Cell *map);

//    IndexedMap(const IndexedMap& copy) :
//            key_bits_(copy.key_bits_),
//            value_bits_(copy.value_bits_),
//            load_factor_(copy.load_factor_),
//            offset_bits_(copy.offset_bits_),
//            buckets_count_(copy.buckets_count_),
//            capacity_(copy.capacity_),
//            size_(copy.size_),
//            offset_(copy.offset_),
//            map_(copy.map_),
//            psl_(copy.psl_) {}

    // No assignment operator

public:
    uint64_t* offset_;

    const size_t key_bits_ = 0;
    const size_t value_bits_ = 0;

    const double load_factor_ = 0.f;

    const size_t offset_bits_ = 0;
    const size_t buckets_count_ = 0;

    const bool multimap_ = false;

    using InsertKeyExists = std::function<size_t (const size_t, const size_t)>;
    InsertKeyExists function_;

    IndexedMap(size_t key_bits, size_t value_bits, size_t offset_bits, size_t* offset_table, double load_factor, bool multimap=false);
    IndexedMap(size_t key_bits, size_t value_bits, size_t offset_bits, std::string bucket_sizes, double load_factor, bool multimap=false);

    ~IndexedMap() {
        // Destructor
        delete[] offset_;
        delete[] map_;
        delete[] psl_;
    }

    void Insert(uint64_t key, uint64_t value);
    void InsertSafe(uint64_t key, uint64_t value);

    Cell* Find(uint64_t &key);

    Cell* begin() {
        return map_;
    }

    Cell* end() {
        return map_ + capacity_;
    }

    void Save(std::string file);
    static IndexedMap* Load(std::string file, bool silent=true);

    static inline uint64_t ComputeOffsetKey(uint64_t key, size_t const key_bits) {
        return key >> key_bits;
    }

    inline uint32_t GetBucketSize(const uint64_t key) {
        auto offset_key = ComputeOffsetKey(key, key_bits_);
        size_t bucket_index = offset_[offset_key];
        return offset_[offset_key + 1] - bucket_index;
    }

    static inline uint64_t ComputeInternalKey(uint64_t key, size_t const key_bits) {
        return key & ((1llu << key_bits) - 1);
    }

    size_t Capacity() const {
        return capacity_;
    }

    size_t K() {
        return (key_bits_ + offset_bits_) / 2;
    }

    size_t OffsetBits() const {
        return offset_bits_;
    }

    size_t KeyBits() const {
        return key_bits_;
    }

    size_t ValueBits() const {
        return value_bits_;
    }

    size_t Size() const {
        return size_;
    }

    inline uint32_t psl(uint64_t key, uint64_t pos, uint32_t bucket_size) {
        uint32_t hash = CompactMapUtils::Hash(key, bucket_size);
        return pos + ((hash > pos) * bucket_size) - hash;
    }

    void CheckBuckets (bool check_uniqueness = false) {
        ProgressBar bar(buckets_count_);
        size_t update_chunk = buckets_count_ / 1000;

        size_t min_size = numeric_limits<size_t>::max();
        size_t max_size = 0;
        double min_load = 1.0;
        double max_load = 0.0;

        size_t empty_buckets = 0;

        tsl::robin_set<size_t> key_set;
        size_t key_count = 0;
        for (int i = 0; i < buckets_count_; i++) {
            if (check_uniqueness)
                key_set.clear();
            key_count = 0;

            if ((i % update_chunk) == 0)
                bar.Update(i+1);

            auto bucket_index = offset_[i];
            auto bucket_size = offset_[i+1] - bucket_index;

            if (bucket_size == 0) ++empty_buckets;

            if (bucket_size > max_size) max_size = bucket_size;
            if (bucket_size < min_size) min_size = bucket_size;

            for (int j = 0; j < bucket_size; j++) {
                auto cell = map_[bucket_index + j];
                if (!cell.empty()) {
                    if (check_uniqueness)
                        key_set.insert(cell.key(value_bits_));
                    key_count++;
                }
            }

            auto load_factor = (double) key_count / bucket_size;
            if (load_factor > max_load) max_load = load_factor;
            if (load_factor < min_load) min_load = load_factor;

            if (check_uniqueness)
                assert(key_count == key_set.size());
        }

        std::cout << " " << std::endl;
        std::cout << "min_size:      \t" << min_size << std::endl;
        std::cout << "max_size:      \t" << max_size << std::endl;
        std::cout << "min_load:      \t" << min_load << std::endl;
        std::cout << "max_load:      \t" << max_load << std::endl;
        std::cout << "empty buckets: \t" << empty_buckets << std::endl;
    }

    void SaveBucketSizes(std::string file) {
        ofstream ofs(file, ios::binary);

        ofs.write((char *) &buckets_count_, sizeof(buckets_count_));

        size_t key_count = 0;

        for (int i = 0; i < buckets_count_; i++) {
            key_count = 0;

            auto& bucket_index = offset_[i];
            auto bucket_size = offset_[i+1] - bucket_index;
            for (int j = 0; j < bucket_size; j++) {
                auto& cell = map_[bucket_index + j];
                key_count += !cell.empty();
            }
            ofs.write((char *) &key_count, sizeof(key_count));
        }
        ofs.close();
    }

//    static pair<size_t*, size_t> LoadBucketSizes(std::string file) {
//        ifstream ifs(file, ios::binary);
//        if (!ifs.is_open()) {
//            std::cerr << "replace by exception" << std::endl;
//        }
//
//        size_t bucket_size;
//        ifs.read((char *) &bucket_size, sizeof(bucket_size));
//
//        size_t* buckets = new size_t[bucket_size];
//
//        ifs.read((char *) buckets, sizeof(buckets) * bucket_size);
//        ifs.close();
//
//        return { buckets, bucket_size };
//    }


    void LoadBucketSizes(std::string file) {
        if (map_ != nullptr) {
            delete[] map_;
        }
        ifstream ifs(file, ios::binary);
        if (!ifs.is_open()) {
            std::cerr << "replace by exception" << std::endl;
        }

        size_t buckets_count = 0;
        size_t bucket_size = 0;

        ifs.read((char *) &buckets_count, sizeof(buckets_count));

        offset_ = new size_t[buckets_count + 1];
        std::fill_n(offset_, buckets_count + 1, 0);

        size_t total = 0;
        offset_[0] = 0;
        for (size_t i = 1; i <= buckets_count; i++) {
            bucket_size = 0;
            ifs.read((char *) &bucket_size, sizeof(bucket_size));
            total += ComputeBucketCapacity2(bucket_size, load_factor_);
            offset_[i] += total;
        }
        capacity_ = total;
        map_ = new Cell[capacity_];
        memset(map_, 0, capacity_*8);
    }

    double ComputeBucketCapacity(size_t bucket_size, double load_factor) {
        return ceil((1/load_factor) * (double) bucket_size);
    }

    double ComputeBucketCapacity2(size_t bucket_size, double load_factor) {
        return (1/load_factor) * (double) bucket_size;
//        return 5 * (double) bucket_size;
    }

    void PrintMeta() {
        std::cout << " " << std::endl;
        std::cout << "### Meta Data for Map:\n\n";
        std::cout << "Capacity:        \t" << capacity_ << std::endl;
        std::cout << "Size:            \t" << size_ << std::endl;
        std::cout << "Offset Size:     \t" << buckets_count_ << std::endl;

        std::cout << "\nKey Bits:        \t" << key_bits_ << std::endl;
        std::cout << "Value Bits:      \t" << value_bits_ << std::endl;
        std::cout << "Offset Bits:     \t" << offset_bits_ << std::endl;
        std::cout << "Key Total Bits:  \t" << (offset_bits_ + key_bits_) << std::endl;
        std::cout << "\nLoad Factor:     \t" << load_factor_ << std::endl;
        std::cout << "Real Load Factor:\t" << ((double) size_/capacity_) << std::endl;

        auto mem_map = (capacity_ * 8) / (1024 * 1024);
        auto mem_offset = (buckets_count_ * 8) / (1024 * 1024);

        std::cout << "\nMemory:" << std::endl;
        std::cout << "Map (MB):        \t" << mem_map << std::endl;
        std::cout << "Map Offset (MB): \t" << mem_offset << std::endl;
        std::cout << "Total (MB):      \t" << (mem_offset + mem_map) << std::endl;
        std::cout << std::endl;
    }

    Cell* Map();

    /* #####################################################################
     * # Debug
     * #####################################################################
     */

    void printPSL() {
        for (int i = 0; i < 100; i++) {
            std::cout << i << ": " << psl_[i] << std::endl;
            if (!psl_[i]) break;
        }
    }

    void print() {
        size_t offset_index = 0;
        size_t offset_value_1 = 0;
        size_t offset_value_2 = 0;
        size_t size = offset_[offset_index] - offset_[offset_index - 1];
        for (int i = 0; i < capacity_; i++) {
            while (offset_value_2 == i) {
                std::cout << "--- " << offset_index << " --------------------------------------------------------------- size: ";
                offset_index++;
                offset_value_1 = offset_value_2;
                offset_value_2 = offset_[offset_index];
                size = offset_value_2 - offset_value_1;
                std::cout << size << std::endl;
            }
            auto c = map_[i];
            auto key = c.key(value_bits_);
            auto value = c.value(value_bits_);

            if (c.empty(value_bits_)) std::cout << (i - offset_value_1) << std::endl;
            else {
                auto hash = CompactMapUtils::Hash(key, (uint32_t) size);
                auto pos = (i - offset_value_1);
                std::cout << pos << ": " << KmerUtils::ToString(offset_index-1, offset_bits_) << " " << KmerUtils::ToString(key, key_bits_) << ": " << value << "\t\tHasholasho: " <<  hash << "\t\tpsl: " << psl(key, pos, size )  << std::endl;
            }

        }
    }

    void print(size_t pos, size_t length, bool print_empty=true) {
        std::cout << "print valuebits: " << value_bits_ << std::endl;
        size_t offset_index = 0;
        while (true) {
            if (offset_[offset_index] > pos) {
                offset_index--;
                break;
            }
            offset_index++;
        }
        size_t offset_value_1 = offset_[offset_index];
        size_t offset_value_2 = offset_[offset_index+1];
        size_t size = offset_value_2 - offset_value_1;

        for (int i = pos; i < capacity_ && length > 0; i++, length--) {

            while (offset_value_2 == i) {
                offset_index++;
                std::cout << "--- " << offset_index << " --------------------------------------------------------------- size: ";
                offset_value_1 = offset_value_2;
                offset_value_2 = offset_[offset_index+1];
                size = offset_value_2 - offset_value_1;
                std::cout << size << std::endl;
            }
            auto c = map_[i];
            auto key = c.key(value_bits_);
            auto value = c.value(value_bits_);

            if (c.empty(value_bits_)) std::cout << (i - offset_value_1) << std::endl;
            else {
                auto hash = CompactMapUtils::Hash(key, (uint32_t) size);
                auto pos = (i - offset_value_1);
                std::cout << pos << ": " << KmerUtils::ToString(offset_index, offset_bits_) << " " << key << " " << KmerUtils::ToString(key, key_bits_) << ": " << value << "\t\tHash: " <<  hash << "\t\tpsl: " << psl(key, pos, size ) << "\t\tdata: " << c.data << std::endl;
            }

        }
    }


    size_t Size2() const {
        size_t total = 0;
        for (int i = 0; i < capacity_; i++) {
            if (!map_[i].empty(value_bits_)) total++;
        }
        return total;
    }

    IndexedMap();

    void InsertPlus(uint64_t key, uint64_t value, function<size_t(const size_t, const size_t)> &function);

    size_t CalcSize();
};