//
// Created by fritsche on 17/02/2022.
//

#pragma once

#include <cstdint>
#include <cstddef>
#include <assert.h>
#include <string>
#include <bitset>
#include <iostream>
#include <tuple>
#include <fstream>
#include "Utils.h"
#include <bit>

namespace protal {
    template<uint64_t taxid_bits, uint64_t geneid_bits, uint64_t genepos_bits>
    struct Entry {
        uint64_t value = 0;

        inline void Put(uint64_t taxid, uint64_t geneid, uint64_t genepos) {
            value = 0;
            value |= taxid << (geneid_bits + genepos_bits);
            value |= geneid << genepos_bits;
            value |= genepos;
        }

        void Get(uint64_t &taxid, uint64_t &geneid, uint64_t &genepos) const {
            taxid = (value >> (geneid_bits + genepos_bits)) & ((1u << taxid_bits) - 1);
            geneid = (value >> (genepos_bits)) & ((1u << geneid_bits) - 1);
            genepos = value & ((1u << genepos_bits) - 1);
        }

        inline bool operator > (Entry const& other) const {
            auto [this_taxid, this_geneid, this_genepos] = Get();
            auto [other_taxid, other_geneid, other_genepos] = other.Get();
            if (this_taxid != other_taxid) {
                return this_taxid > other_taxid;
            }
            if (this_geneid != other_geneid) {
                return this_geneid > other_geneid;
            }
            return false;
        }

        inline bool operator < (Entry const& other) const {
            auto [this_taxid, this_geneid, this_genepos] = Get();
            auto [other_taxid, other_geneid, other_genepos] = other.Get();
            if (this_taxid != other_taxid) {
                return this_taxid < other_taxid;
            }
            if (this_geneid != other_geneid) {
                return this_geneid < other_geneid;
            }
            return false;
        };

//        inline bool operator > (LookupResult const& other) const = default;
        inline int operator == (Entry const& other) const {
            auto [this_taxid, this_geneid, this_genepos] = Get();
            auto [other_taxid, other_geneid, other_genepos] = other.Get();
            return this_taxid == other_taxid && this_geneid == other_geneid;
        }

        std::tuple<uint64_t, uint64_t, uint64_t> Get() const {
            auto taxid = (value >> (geneid_bits + genepos_bits)) & ((1u << taxid_bits) - 1);
            auto geneid = (value >> genepos_bits) & ((1u << geneid_bits) - 1);
            auto genepos = value & ((1u << genepos_bits) - 1);
            return { taxid, geneid, genepos };
        }

        std::string ToString() {
            std::string result;
            auto [taxid, geneid, genepos] = Get();
            result.reserve(60);
            result += "taxid:\t" + std::to_string(taxid) + "\t\t";
            result += "geneid:\t" + std::to_string(geneid) + "\t\t";
            result += "genepos:\t" + std::to_string(genepos);
            return result;
        }

        bool Empty() {
            return value == 0;
        }
    };

    class SeedmapUtils {
    public:
        template<size_t bits>
        static std::string BitString(uint64_t key) {
            std::bitset<bits> x(key);
            return x.to_string();
        }

        static inline void ParseHeader(std::string &header, size_t &internal_taxid, size_t &marker_gene_id) {
            auto tokens = Utils::split(header.substr(1, std::string::npos), "_");
            assert(tokens.size() > 0);
            internal_taxid = stoull(tokens[0]);
            marker_gene_id = (tokens.size() > 1) * stoull(tokens[1]);
        }
    };

//    │  values table                 │
//    │                               │
//    │                               │
//    ├───────────────┬───────────────┤ ─┐
//    │ 1  TCACACGTC  │ 2 GTCACACGC   │  │
//    ├───────────────┼───────────────┤  │   Here are the flexi-k values stored
//    │ 3  ATGCATGCT  │ 4 CGACTCGGC   │  │   Use a loop and similarity check
//    ├───────────────┼───────────────┤  │   Maybe sth with popcount to find
//    │ 5  ...        │ 6  ...        │  │   best matching k-mer in the value-section
//    ├───────────────┼───────────────┤  │
//    │ 7  ...        │ 8  ...        │  │
// ┌─ ├───────────────┴───────────────┤ ─┘
// │  │ 1 ( taxid, geneid, genepos )  │
// │  ├───────────────────────────────┤
// │  │ 2                             │
// │  ├───────────────────────────────┤
// │  │               .               │
// │  │               .               │
// │  │               .               │
// │  ├───────────────────────────────┤
// │  │ 7                             │
// │  ├───────────────────────────────┤
// │  │ 8                             │
// └─ ├───────────────────────────────┤
//    │                               │


    class Seedmap {
    public:
        // If a 15-mer has more than 8 locations, use flexi-k approach
        size_t m_exact_k = 15;
        size_t m_main_bits = m_exact_k * 2;
        size_t m_flex_threshold = 2;
        size_t m_flex_k = 16;
        size_t m_flex_k_half = m_flex_k/2;
        size_t m_flex_k_bits = m_flex_k*2;
        size_t m_main_key_mask = ((1llu << m_main_bits) - 1) << m_flex_k; // First mask, then shift
        size_t m_flex_key_mask_left = ((1llu << m_flex_k) - 1) << (m_flex_k + m_main_bits);
        size_t m_flex_key_mask_right = (1llu << m_flex_k) - 1;

    private:

//        using KeyMap_t = uint8_t;
        using KeyMap_t = uint16_t; // CHANGE THIS IN LATEST VERSION
        KeyMap_t* m_keymap = nullptr;
        uint64_t m_keymask = 0b0000000000000000000000000000000000111111111111111111111111111111;
        size_t keymap_size = 1u << m_main_bits;
        size_t keymap_max = m_keymask;

        // This variable defines how many keys are managed by one control block
        // Must be power of two
//        size_t m_keys_per_ctrl_block = 16; // Version 1
        size_t m_keys_per_ctrl_block = 8; // Version 2
//        size_t m_keys_per_ctrl_block = 32; // Flex 16-bit Cells

        // bitshift to find which control block a key is in.
        size_t ctrl_block_frequency_bitshift = log2(m_keys_per_ctrl_block);

        uint64_t ctrl_block_key_mask = (1 << ctrl_block_frequency_bitshift) - 1;

        // (1 << 8) = 256. This is how many fields an 8-bit number can index.
        // Spread the indexing power equally across the number of keys per block
        size_t max_key_ubiquity = (1 << (sizeof(KeyMap_t)*8)) / m_keys_per_ctrl_block;

//        size_t ctrl_block_frequency_bitshift = 4;

        // Bytespace a controlblock takes up
        size_t ctrl_block_byte_size = 4; // Byte
        size_t ctrl_block_cell_size = ctrl_block_byte_size/sizeof(KeyMap_t); // Byte
        size_t ctrl_block_byte_size_shift = log2(ctrl_block_cell_size);
//        size_t ctrl_block_byte_size_shift = 2;

        size_t keymap_size_total = keymap_size + (((keymap_size / m_keys_per_ctrl_block) + 1) * ctrl_block_cell_size);

        /*
         * The number of total value entries per block.
         * Each block holds indices for <m_keys_per_ctrl_block> keys
         * and the values for thes keys may not be more than <max_block_size>
         */
        size_t max_block_size = (1 << sizeof(KeyMap_t)*8); //  This is for block
        size_t max_key_multiplicity = 2048; //  This is for key, arbitrary value
        size_t values_size= 0;
        Entry<20,20,20>* m_map = nullptr;

        uint64_t m_found_counter = 0;

        void Print() {
            std::cout << "m_keys_per_ctrl_block:               " << m_keys_per_ctrl_block << std::endl;
            std::cout << "ctrl_block_frequency_bitshift:       " << ctrl_block_frequency_bitshift << std::endl;
            std::cout << "ctrl_block_byte_size:                " << ctrl_block_byte_size << std::endl;
            std::cout << "ctrl_block_cell_size:                " << ctrl_block_cell_size << std::endl;
            std::cout << "ctrl_block_byte_size_shift:          " << ctrl_block_byte_size_shift << std::endl;
            std::cout << "ctrl_block_key_mask:                 " << ctrl_block_key_mask << std::endl;
            std::cout << "keymap_size_total:                   " << keymap_size_total << std::endl;
            std::cout << "max_block_size:                      " << max_block_size << std::endl;
            std::cout << "values_size:                         " << values_size << std::endl;
        }

        Seedmap(Seedmap const& other) = delete;
    public:
        Seedmap() {
            m_keymap = new KeyMap_t[keymap_size_total];
        }

        Seedmap(std::string file) {
            Load(file);
        };

        ~Seedmap() {
            delete[] m_keymap;
            delete[] m_map;
        }

        inline uint64_t KeymapIndex(uint64_t key) const {
            return key + ((key >> ctrl_block_frequency_bitshift) << ctrl_block_byte_size_shift) + ctrl_block_cell_size;
        }

        inline uint64_t ControlBlockIndex(uint64_t key) const {
            return ((key >> ctrl_block_frequency_bitshift) << ctrl_block_frequency_bitshift) + ((key >> ctrl_block_frequency_bitshift) << (ctrl_block_byte_size_shift));
        }

        inline uint64_t BlockKey(uint64_t key) const {
            return key & ctrl_block_key_mask;
        }

        inline static uint64_t DivisionByTwoCeiling(uint64_t x) {
            return (x & 1llu) ? (x >> 1) + 1 : (x >> 1);
        }

        void Save(std::string file) {

//            SortForKeys();

            std::ofstream ofs(file, std::ios::binary);
            ofs.write((char *) &keymap_size, sizeof(keymap_size));
            ofs.write((char *) &keymap_size_total, sizeof(keymap_size_total));
            ofs.write((char *) &ctrl_block_byte_size, sizeof(ctrl_block_byte_size));
            ofs.write((char *) &m_keys_per_ctrl_block, sizeof(m_keys_per_ctrl_block));
            ofs.write((char *) &ctrl_block_frequency_bitshift, sizeof(ctrl_block_frequency_bitshift));

            ofs.write((char *) &values_size, sizeof(values_size));

            ofs.write((char *) m_keymap, sizeof(*m_keymap) * keymap_size_total);
            ofs.write((char *) m_map, sizeof(*m_map) * values_size);

            ofs.close();
        }

        void Save(std::ostream& ofs) {
//            SortForKeys();

            ofs.write((char *) &keymap_size, sizeof(keymap_size));
            ofs.write((char *) &keymap_size_total, sizeof(keymap_size_total));
            ofs.write((char *) &ctrl_block_byte_size, sizeof(ctrl_block_byte_size));
            ofs.write((char *) &m_keys_per_ctrl_block, sizeof(m_keys_per_ctrl_block));
            ofs.write((char *) &ctrl_block_frequency_bitshift, sizeof(ctrl_block_frequency_bitshift));

            ofs.write((char *) &values_size, sizeof(values_size));

            ofs.write((char *) m_keymap, sizeof(*m_keymap) * keymap_size_total);
            ofs.write((char *) m_map, sizeof(*m_map) * values_size);
        }

        void Load(std::istream &ifs) {
            ifs.read((char *) &keymap_size, sizeof(keymap_size));
            ifs.read((char *) &keymap_size_total, sizeof(keymap_size_total));
            ifs.read((char *) &ctrl_block_byte_size, sizeof(ctrl_block_byte_size));
            ifs.read((char *) &m_keys_per_ctrl_block, sizeof(m_keys_per_ctrl_block));
            ifs.read((char *) &ctrl_block_frequency_bitshift, sizeof(ctrl_block_frequency_bitshift));
            ifs.read((char *) &values_size, sizeof(values_size));

            delete[] m_keymap;
            m_keymap = new KeyMap_t[keymap_size_total];
            ifs.read((char *) m_keymap, sizeof(*m_keymap) * (keymap_size_total));

            m_map = new Entry<20,20,20>[values_size];
            ifs.read((char *) m_map, sizeof(*m_map) * (values_size));
        }

        void Load(std::string file) {

            std::ifstream ifs(file, std::ios::binary);

            ifs.read((char *) &keymap_size, sizeof(keymap_size));
            ifs.read((char *) &keymap_size_total, sizeof(keymap_size_total));
            ifs.read((char *) &ctrl_block_byte_size, sizeof(ctrl_block_byte_size));
            ifs.read((char *) &m_keys_per_ctrl_block, sizeof(m_keys_per_ctrl_block));
            ifs.read((char *) &ctrl_block_frequency_bitshift, sizeof(ctrl_block_frequency_bitshift));
            ifs.read((char *) &values_size, sizeof(values_size));

            m_keymap = new KeyMap_t[keymap_size_total];
            ifs.read((char *) m_keymap, sizeof(*m_keymap) * (keymap_size_total));

            m_map = new Entry<20,20,20>[values_size];
            ifs.read((char *) m_map, sizeof(*m_map) * (values_size));
        }


        void CountUpKey(uint64_t key) {
            assert(key < keymap_size);
            auto index = KeymapIndex(key);
            m_keymap[index] += m_keymap[index] < UINT16_MAX;
        }

        void CountUpKeyOMP(uint64_t key) {
            assert(key < keymap_size);
            auto index = KeymapIndex(key);

#pragma omp atomic
            m_keymap[index] += m_keymap[index] < UINT16_MAX;
        }

        uint64_t GetKey(uint64_t key) {
            return m_keymap[KeymapIndex(key)];
        }

        uint8_t Key(uint64_t key) {
            return m_keymap[KeymapIndex(key)];
        }

        void Get(uint64_t key, Entry<20, 20, 20>* &start, Entry<20, 20, 20>* &end, uint32_t* &flexblock_begin, uint32_t* &flexblock_end) {
            size_t main_key = MainKey(key);
            size_t flex_key = FlexKey(key);

            if (main_key > m_keymask) {
                std::cout << "Key is larger than keymask" << key << " > " << m_keymask << std::endl;
            }
            uint64_t block_start_idx = ControlBlockIndex(main_key);
            uint64_t block_end_idx = block_start_idx + m_keys_per_ctrl_block + ctrl_block_cell_size;

            uint64_t block_value_start_idx = *((uint32_t *) (m_keymap + block_start_idx));
            uint64_t block_value_end_idx = *((uint32_t *) (m_keymap + block_end_idx));
            uint64_t block_value_size = block_value_end_idx - block_value_start_idx;

            if (!block_value_size) {
                start = nullptr;
                end = nullptr;
                return;
            }

            size_t key_index = block_start_idx + ctrl_block_cell_size + BlockKey(main_key);

            auto key_value_start = block_value_start_idx + m_keymap[key_index];
            auto key_value_end = (key_index + 1) == block_end_idx ? block_value_end_idx : block_value_start_idx +
                                                                                          m_keymap[key_index + 1];
            auto key_value_block_size = key_value_end - key_value_start;

            if (!key_value_block_size) {
                start = nullptr;
                end = nullptr;
                return;
            }

            flexblock_begin = nullptr;
            flexblock_end = nullptr;
            if (key_value_block_size >= m_flex_threshold) {
                flexblock_begin = (uint32_t*) (m_map + key_value_start);
                flexblock_end = flexblock_begin + key_value_block_size - FlexBlockSize(key_value_block_size);
                key_value_start += FlexBlockSize(key_value_block_size);
            }

            start = m_map + key_value_start;
            end =  m_map + key_value_end;
        }

        void Get(uint64_t key, Entry<20, 20, 20>* &start, Entry<20, 20, 20>* &end) {
            size_t main_key = MainKey(key);
            size_t flex_key = FlexKey(key);

            if (main_key > m_keymask) {
                std::cout << "Key is larger than keymask" << key << " > " << m_keymask << std::endl;
            }
            uint64_t block_start_idx = ControlBlockIndex(main_key);
            uint64_t block_end_idx = block_start_idx + m_keys_per_ctrl_block + ctrl_block_cell_size;

            uint64_t block_value_start_idx = *((uint32_t *) (m_keymap + block_start_idx));
            uint64_t block_value_end_idx = *((uint32_t *) (m_keymap + block_end_idx));
            uint64_t block_value_size = block_value_end_idx - block_value_start_idx;

            if (!block_value_size) {
                start = nullptr;
                end = nullptr;
                return;
            }

            size_t key_index = block_start_idx + ctrl_block_cell_size + BlockKey(main_key);

            auto key_value_start = block_value_start_idx + m_keymap[key_index];
            auto key_value_end = (key_index + 1) == block_end_idx ? block_value_end_idx : block_value_start_idx +
                                                                                          m_keymap[key_index + 1];
            auto key_value_block_size = key_value_end - key_value_start;

            if (!key_value_block_size) {
                start = nullptr;
                end = nullptr;
                return;
            }

            if (key_value_block_size >= m_flex_threshold) {
                key_value_start += FlexBlockSize(key_value_block_size);
            }
//            if (key_value_block_size >= m_flex_threshold) {
//                std::cout << key_value_block_size << std::endl;
//                PrintBlock(main_key);
//
//                uint32_t* flexmer_ptr = (uint32_t*) (m_map + key_value_start);
//                // Update.. also give out positions of extended keys.
//                key_value_start += FlexBlockSize(key_value_block_size);
//
//                auto values_size = key_value_end - key_value_start;
//
//                std::cout << "Flexmer: " << KmerUtils::ToString(flex_key, m_flex_k_bits) << " " << SeedmapUtils::BitString<64>(flex_key) << std::endl;
//                for (auto i = 0; i < values_size; i++) {
//                    std::cout << " -> " << KmerUtils::ToString(flexmer_ptr[i], m_flex_k_bits) << " " << Similarity(flex_key, flexmer_ptr[i]) << " " << SeedmapUtils::BitString<64>(flexmer_ptr[i]) << std::endl;
//                }
//
//                Utils::Input();
//            }


            if (key_value_start > key_value_end) {
                std::cout << "Block end offset: " << m_keys_per_ctrl_block + ctrl_block_cell_size << std::endl;
                std::cout << "block_value_start_idx: " << block_value_start_idx << std::endl;
                std::cout << "block_value_end_idx: " << block_value_end_idx << std::endl;
                std::cout << "block_start_idx: " << block_start_idx << std::endl;
                std::cout << "key_index: " << key_index << std::endl;
                std::cout << "BlockKey(main_key): " << BlockKey(main_key) << std::endl;
                std::cout << main_key << " " << KmerUtils::ToString(main_key, m_main_bits) << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;
                std::cout << "key_value_start: " << key_value_start << std::endl;
                std::cout << "key_value_end: " << key_value_end << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;
//                PrintBlock(main_key);
                std::cout << "Seedmap: 240" << std::endl;
                Print();
                exit(71);
            }

            start = m_map + key_value_start;
            end =  m_map + key_value_end;
        }

        bool Put(uint64_t &key, uint64_t &taxid, uint64_t &geneid, uint64_t &genepos) {
            if (m_map == nullptr) {
                std::cerr << "You need to initialize the values with BuildValuePointers()" << std::endl;
                exit(8);
            }
            if (key > keymap_max) {
                std::cerr << "Key too large (" << key << " > " << keymap_max << ")" << std::endl;
                exit(8);
            }


            uint64_t block_start_idx = ControlBlockIndex(key);
            uint64_t block_end_idx = block_start_idx + m_keys_per_ctrl_block + ctrl_block_cell_size;

            uint64_t block_value_start_idx = *((uint32_t*) (m_keymap + block_start_idx ) );
            uint64_t block_value_end_idx = *((uint32_t*) (m_keymap + block_end_idx) );
            uint64_t block_value_size = block_value_end_idx - block_value_start_idx;

            if (!block_value_size) {
                return false;
            }

            size_t key_index = block_start_idx + ctrl_block_cell_size + BlockKey(key);
            auto key_value_start = block_value_start_idx + m_keymap[key_index];
            auto key_value_end = (key_index + 1) == block_end_idx ? block_value_end_idx : block_value_start_idx + m_keymap[key_index + 1];
            auto key_value_block_size = key_value_end - key_value_start;

            if (!key_value_block_size) {
                return false;
            }

            auto index = key_value_start;

            if (key_value_start > key_value_end) {
                std::cout << "block_value_start_idx: " << block_value_start_idx << std::endl;
                std::cout << "block_value_end_idx: " << block_value_end_idx << std::endl;
                std::cout << "block_start_idx: " << block_start_idx << std::endl;
                std::cout << "key_index: " << key_index << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;
                std::cout << "key_value_start: " << key_value_start << std::endl;
                std::cout << "key_value_end: " << key_value_end << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;


//                PrintBlock(key);
                std::cout << "Key value start (" << key_value_start << ") can not be larger than key value end (" << key_value_end << ")" << std::endl;
                exit(72);
            }

            for (; !m_map[index].Empty() && index < key_value_end; index++);

            if (index == key_value_end) {
                std::cout << "bucket too large, not present " << key_value_start << " size: " << (key_value_end - key_value_start) << std::endl;
                return false;
            }
            m_map[index].Put(taxid, geneid, genepos);

            Entry<20,20,20> *begin, *end;
            Get(key, begin, end);

            if constexpr(false) {
//            std::cout << "put index: " << index << " pointer: " << &m_map[index] << " begin pointer " << begin  << " " << taxid << " " << genepos << " " << geneid << std::endl;
                bool found = false;
                for (auto it = begin; it < end; it++) {
//                std::cout << it->ToString() << std::endl;
                    auto [t, g, p] = it->Get();
                    if (taxid == t && geneid == g && genepos == p) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    std::cout << "____searched: " << taxid << ", " << geneid << ", " << genepos << std::endl;
                    for (auto it = begin; it < end; it++) {
                        std::cout << it->ToString() << std::endl;
                    }
                    exit(9);
                } else {
                    m_found_counter++;
                    if (m_found_counter % 10'000'000 == 0)
                        std::cout << "found counter: " << m_found_counter << std::endl;
                }
            }

//            std::string stop;
//            std::cin >> stop;
            return true;
        }

        uint64_t MainKey(uint64_t key) {
            return (key & m_main_key_mask) >> m_flex_k ; // shift by half of flexbits which is flex_k
        }

        uint64_t FlexKey(uint64_t key) {
            return ((key & m_flex_key_mask_left) >> (m_main_bits)) | (key & m_flex_key_mask_right) ; // shift by half of flexbits which is flex_k
        }

        uint64_t FlexBlockSize(uint64_t block_size) {
            return (block_size + 2) / 3;
        }

//        static uint32_t Similarity(uint32_t a, uint32_t b) {
//            return std::popcount(((~(a ^ b) >> 1) & ~(a ^ b)) & 0b00000000010101010101010101010101);
//        }
        static uint32_t Similarity(uint32_t a, uint32_t b) {
            return std::popcount(((~(a ^ b) >> 1) & ~(a ^ b)) & 0b01010101010101010101010101010101);
        }

        bool PutOMP(uint64_t &key, uint64_t &taxid, uint64_t &geneid, uint64_t &genepos) {
            uint64_t main_key = MainKey(key);
            uint64_t flex_key = FlexKey(key);

            static size_t count_flex = 0;
            if (m_map == nullptr) {
                std::cerr << "You need to initialize the values with BuildValuePointers()" << std::endl;
                exit(8);
            }
            if (main_key > keymap_max) {
                std::cerr << "Key too large (" << main_key << " > " << keymap_max << ")" << std::endl;
                exit(8);
            }


            uint64_t block_start_idx = ControlBlockIndex(main_key);
            uint64_t block_end_idx = block_start_idx + m_keys_per_ctrl_block + ctrl_block_cell_size;

            uint64_t block_value_start_idx = *((uint32_t*) (m_keymap + block_start_idx ) );
            uint64_t block_value_end_idx = *((uint32_t*) (m_keymap + block_end_idx) );
            uint64_t block_value_size = block_value_end_idx - block_value_start_idx;

            if (!block_value_size) {
                return false;
            }


            size_t key_index = block_start_idx + ctrl_block_cell_size + BlockKey(main_key);
            auto key_value_start = block_value_start_idx + m_keymap[key_index];
            auto key_value_end = (key_index + 1) == block_end_idx ? block_value_end_idx : block_value_start_idx + m_keymap[key_index + 1];
            auto key_value_block_size = key_value_end - key_value_start;

            if (!key_value_block_size) {
                return false;
            }

            bool is_flex = key_value_block_size >= m_flex_threshold;
            size_t flex_block_size = FlexBlockSize(key_value_block_size);
            auto value_index = key_value_start + (is_flex * flex_block_size);

            if (value_index >= values_size) {
                std::cout << "key: " << key << std::endl;
                std::cout << "key: " << SeedmapUtils::BitString<64>(key) << std::endl;
                std::cout << "key: " << SeedmapUtils::BitString<64>(main_key) << std::endl;
                std::cout << "key: " << SeedmapUtils::BitString<64>(flex_key) << std::endl;
                std::cout << "block_value_start_idx: " << block_value_start_idx << std::endl;
                std::cout << "block_value_end_idx: " << block_value_end_idx << std::endl;
                std::cout << "block_start_idx: " << block_start_idx << std::endl;
                std::cout << "block_end_idx: " << block_end_idx << std::endl;
                std::cout << "key_index: " << key_index << std::endl;
                std::cout << "BlockKey(main_key): " << BlockKey(main_key) << std::endl;
                std::cout << main_key << " " << KmerUtils::ToString(main_key, m_main_bits) << std::endl;
                std::cout << flex_key << " " << KmerUtils::ToString(flex_key, m_flex_k_bits) << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;
                std::cout << "is_flex: " << is_flex << std::endl;
                std::cout << "flex_block_size: " << flex_block_size << std::endl;
                std::cout << "key_value_start: " << key_value_start << std::endl;
                std::cout << "key_value_end: " << key_value_end << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;
                std::cout << "Value index " << value_index << "/" << values_size << std::endl;

//                PrintBlock(main_key);
                exit(73);
            }

            // Find next empty slot
            for (; value_index < key_value_end && !m_map[value_index].Empty(); value_index++);
            if (value_index == key_value_end) {
                std::cout << "block_value_start_idx: " << block_value_start_idx << std::endl;
                std::cout << "block_value_end_idx: " << block_value_end_idx << std::endl;
                std::cout << "block_start_idx: " << block_start_idx << std::endl;
                std::cout << "key_index: " << key_index << std::endl;
                std::cout << "BlockKey(main_key): " << BlockKey(main_key) << std::endl;
                std::cout << main_key << " " << KmerUtils::ToString(main_key, m_main_bits) << std::endl;
                std::cout << flex_key << " " << KmerUtils::ToString(flex_key, m_flex_k_bits) << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;
                std::cout << "is_flex: " << is_flex << std::endl;
                std::cout << "flex_block_size: " << flex_block_size << std::endl;
                std::cout << "key_value_start: " << key_value_start << std::endl;
                std::cout << "key_value_end: " << key_value_end << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;
                std::cout << "Value index " << value_index << "/" << values_size << std::endl;
//                PrintBlock(main_key);
                exit(74);
            }

#pragma omp critical(put)
{
            while(value_index < key_value_end && !m_map[value_index].Empty()) value_index++;
            if (value_index == key_value_end) {
                std::cout << "block_value_start_idx: " << block_value_start_idx << std::endl;
                std::cout << "block_value_end_idx: " << block_value_end_idx << std::endl;
                std::cout << "block_start_idx: " << block_start_idx << std::endl;
                std::cout << "key_index: " << key_index << std::endl;
                std::cout << "BlockKey(main_key): " << BlockKey(main_key) << std::endl;
                std::cout << main_key << " " << KmerUtils::ToString(main_key, m_main_bits) << std::endl;
                std::cout << flex_key << " " << KmerUtils::ToString(flex_key, m_flex_k_bits) << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;
                std::cout << "is_flex: " << is_flex << std::endl;
                std::cout << "flex_block_size: " << flex_block_size << std::endl;
                std::cout << "key_value_start: " << key_value_start << std::endl;
                std::cout << "key_value_end: " << key_value_end << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;
                std::cout << "Value index " << value_index << "/" << values_size << std::endl;
//                PrintBlock(main_key);
                exit(75);
            }

//            if (is_flex) {
//                count_flex++;
//                uint64_t flexpos = value_index - key_value_start - flex_block_size;
//                uint32_t* flex_ptr = (uint32_t*) (m_map + key_value_start);
//                std::cout << "GRAND INDEX: " << key_value_start << "  + FLEXPOS: " << flexpos << " INSERT SIDEKEY " << KmerUtils::ToString(flex_key, 24) << std::endl;
//                std::cout << "sidekey: " << flex_key << std::endl;
//                std::cout << "before flex_ptr[flexpos]: " << flex_ptr[flexpos] << std::endl;
//                flex_ptr[flexpos] = static_cast<uint32_t>(flex_key);
//                std::cout << "after  flex_ptr[flexpos]: " << flex_ptr[flexpos] << std::endl;
//
////                if (flexpos == key_value_block_size - flex_block_size - 1) {
//                    std:cout << "Before Put" << std::endl;
//                    std::cout << "Put: " << taxid << ", " << geneid << ", " << genepos << "  at: " << value_index << std::endl;
//                    std::cout << "Started search at index : " << (key_value_start + (is_flex * flex_block_size)) << " flex start : " << key_value_start << std::endl;
//                    PrintBlock(main_key);
////                }
//
////                if ((count_flex % 1000) == 0) {
////                    std::cout << count_flex << std::endl;
////                }
//            }

            // put in at value_index:
            // Need to add m_flex_k to genepos because the main exact match sits in the middle
            m_map[value_index].Put(taxid, geneid, genepos + m_flex_k_half);

            if (is_flex) {
                count_flex++;
                uint64_t flexpos = value_index - key_value_start - flex_block_size;
                uint32_t* flex_ptr = (uint32_t*) &(m_map[key_value_start]);

//                assert(flex_ptr[flexpos] == 0);
                if (flex_ptr[flexpos] != 0) {
                    std::cout << "m_map[" << key_value_start << "]: " << m_map[key_value_start].value << std::endl;
                    std::cout << ((uint32_t *) (m_map + key_value_start))[0] << " "
                              << ((uint32_t *) (m_map + key_value_start))[1] << std::endl;
                    std::cout << KmerUtils::ToString(((uint32_t *) (m_map + key_value_start))[0], m_flex_k_bits) << " "
                              << KmerUtils::ToString(((uint32_t *) (m_map + key_value_start))[1], m_flex_k_bits) << std::endl;
                    std::cout << "m_map[key_value_start]empty?: " << m_map[key_value_start].Empty() << std::endl;

                    if (flexpos == key_value_block_size - flex_block_size - 1) {
                        std::cout << "Main key:          " << KmerUtils::ToString(main_key, m_main_bits) << std::endl;
                        std::cout << "Side key:          " << KmerUtils::ToString(flex_key, m_flex_k_bits) << std::endl;
                        std::cout << "At:                " << flexpos + 1 << "/"
                                  << key_value_block_size - flex_block_size << std::endl;
                        std::cout << "With Flexpos-size: " << flex_block_size << std::endl;
                        std::cout << "Block size:        " << key_value_block_size - flex_block_size << std::endl;
                        std::cout << "key_value_start:   " << key_value_start << std::endl;
                        std::cout << "flex_start:        " << (key_value_start + (is_flex * flex_block_size))
                                  << std::endl;
                    }
//                    PrintBlock(main_key);
                    Utils::Input();
                }
                flex_ptr[flexpos] = static_cast<uint32_t>(flex_key);
//                auto tuple = m_map[key_value_start + (flexpos >> 1)].Get();
//                if (get<0>(tuple) == 38055 && get<1>(tuple) == 479241 && get<2>(tuple) == 304245) {
//                    PrintBlock(main_key);
//                    std::cout << key_value_start << " flexpos: " << flexpos << std::endl;
//                    std::cout << m_map[key_value_start + (flexpos >> 1)].ToString() << std::endl;
//                    exit(55);
//                }
            }
}

            Entry<20,20,20> *begin, *end;
            Get(main_key, begin, end);

            return true;
        }


        struct subkey {
            uint32_t key = 0;
            uint32_t count = 0;

            void Reset() {
                key = 0; count = 0;
            }
            void Print() {
                std::cout << "key: " << key << "  count: " << count << std::endl;
            }
            std::string ToString() const {
                return std::to_string(key) + '\t' + std::to_string(count);
            }
        };

        void BuildValuePointersFirstIteration(uint64_t& total_count, subkey* subkey, int* kmer_frequencies,
                                              uint64_t block_start_index) {
            // Set ubiquity of each key in subkey
            for (auto key_pos = 0; key_pos < m_keys_per_ctrl_block; key_pos++) {
                subkey[key_pos].key = key_pos;
                subkey[key_pos].count = m_keymap[block_start_index + key_pos];
                kmer_frequencies[subkey[key_pos].count]++;
                total_count += subkey[key_pos].count;
            }

            // Sort so keys are sorted by ubiquity.
            std::sort(subkey, subkey + m_keys_per_ctrl_block, [](const struct subkey &a, const struct subkey &b) {
                return a.count < b.count;
            });
        }

        void BuildValuePointersSecondIteration(subkey* subkey, size_t& count_failed_demand, size_t& count_total_stored,
                                               size_t& count_total_demand, uint64_t& max_keyblock_size) {
            uint64_t current_sum = 0;
            uint64_t total_demand = 0;

            for (auto key_pos = 0; key_pos < m_keys_per_ctrl_block; key_pos++) {
                if (subkey[key_pos].count == 0) continue;

                // Given the Flex-K approach block demand quantifies the number of value cells
                // a given key needs.
                size_t key_value_space_demand = subkey[key_pos].count >= m_flex_threshold ?
                                                subkey[key_pos].count + DivisionByTwoCeiling(subkey[key_pos].count) : subkey[key_pos].count;

                total_demand += key_value_space_demand;
//                    subkey[key_pos].Print();
                // This if decides whether keys are taken or not !!!
                // Important to get the key_value_space demand correct otherwise the offsets will be wrong
                // and it will come to SEGMENTATION FAULTS
                if (subkey[key_pos].count < max_key_multiplicity && key_value_space_demand < max_keyblock_size && (current_sum + key_value_space_demand) < max_block_size) {
                    current_sum += key_value_space_demand;
                } else {
                    subkey[key_pos].count = 0;
                }
            }
            // Only statistics.
            count_total_demand += total_demand;
            count_total_stored += current_sum;
            count_failed_demand += (current_sum < total_demand);

            std::sort(subkey, subkey + m_keys_per_ctrl_block, [](const struct subkey &a, const struct subkey &b) {
                return a.key < b.key;
            });
        }

        void BuildValuePointersBlockThirdIteration(subkey* subkey, uint64_t block_start_index, uint32_t& global_position) {
            uint64_t local_position = 0;
            for (auto key_pos = 0; key_pos < m_keys_per_ctrl_block; key_pos++) {
                m_keymap[block_start_index + key_pos] = local_position;
                size_t key_value_space_demand = subkey[key_pos].count >= m_flex_threshold ?
                                                subkey[key_pos].count + DivisionByTwoCeiling(subkey[key_pos].count) : subkey[key_pos].count;

                local_position += key_value_space_demand;
                global_position += key_value_space_demand;
            }
        }

        void BuildValuePointersBlock(uint64_t block, uint32_t& global_position,
                                     size_t& count_failed_demand, size_t& count_total_stored, size_t& count_total_demand,
                                     int* kmer_frequencies, uint64_t& max_keyblock_size, subkey* subkey) {

            uint32_t* control_ptr = nullptr;
            auto block_index = block * (m_keys_per_ctrl_block + ctrl_block_cell_size);
            auto block_start_index = block_index + ctrl_block_cell_size;

            uint64_t total_count = 0;

            control_ptr = (uint32_t*) (m_keymap + block_index);
            *control_ptr = global_position;


            //####################################################################################################
            // first iteration
            // - After this keys in block are sorted by their frequency
            // second iteration
            // - After this, keys in block are resorted by their key.
            // Third iteration
            // - Here, the global offset values and local offset values are set.
            BuildValuePointersFirstIteration(total_count, subkey, kmer_frequencies, block_start_index);
            if (!total_count) return;
            BuildValuePointersSecondIteration(subkey, count_failed_demand, count_total_stored, count_total_demand,
                                              max_keyblock_size);
            BuildValuePointersBlockThirdIteration(subkey, block_start_index, global_position);
        }


        void BuildValuePointers() {
            uint32_t control_idx = 0;
            uint32_t* control_ptr = (uint32_t*) m_keymap;

            uint32_t global_position = 0;
            uint64_t num_ctrl_blocks = keymap_size >> ctrl_block_frequency_bitshift;
            std::cout << "number ctrl blocks: " << num_ctrl_blocks << std::endl;
            uint64_t max_keyblock_size = max_key_ubiquity*2;


            subkey subkey[m_keys_per_ctrl_block];

            // Understand kmer frequencies better to improve sensitivity of alignment.
            constexpr int kmer_freq_size = 1 << (sizeof(KeyMap_t) * 8);
            int kmer_frequencies[kmer_freq_size] = { 0 };
            size_t count_failed_demand = 0;
            size_t count_total_stored = 0;
            size_t count_total_demand = 0;

            // Iterate through blocks. Each block manages <m_keys_per_ctrl_block> 15-mers
            for (auto block = 0; block < num_ctrl_blocks; block++) {
                BuildValuePointersBlock(block, global_position, count_failed_demand, count_total_stored,
                                        count_total_demand, kmer_frequencies, max_keyblock_size, subkey);
            }

            // set last control pointer
            control_ptr = (uint32_t*) (m_keymap + keymap_size_total - ctrl_block_cell_size);
            *control_ptr = global_position;
            values_size = global_position;
            std::cout << "ctrl_block_cell_size: " << ctrl_block_cell_size << std::endl;

            m_map = new Entry<20, 20, 20>[values_size];

            constexpr bool verbose = true;
            if constexpr(verbose) {
                std::cout << "Values size: " << values_size << std::endl;
                std::cout << "space requirements" << std::endl;
                double key_mem = (double) keymap_size_total / (1024 * 1024 * 1024);
                double val_mem = (double) (values_size * 8) / (1024 * 1024 * 1024);
                std::cout << "keys:   " << key_mem << " GB" << std::endl;
                std::cout << "values: " << val_mem << " GB" << std::endl;
                std::cout << "total:  " << key_mem + val_mem << " GB" << std::endl;


                for (auto i = 0; i < kmer_freq_size; i++) {
                    std::cout << i << '\t' << kmer_frequencies[i] << std::endl;
                }
                std::cout << "Total stored: " << count_total_stored << std::endl;
                std::cout << "Total demand: " << count_total_demand << std::endl;
                std::cout << "Stored "
                          << (100 * static_cast<double>(count_total_stored) / static_cast<double>(count_total_demand))
                          << " of total values." << std::endl;
                std::cout << "Failed to meet demands? " << count_failed_demand << std::endl;
            }
        }


        void SortForKeys() {
            std::cout << "Sort for keys" << std::endl;
            size_t row_index = 0;
            std::string line = "";
            size_t non_empty_block_count = 0;
            for (auto key = 0; key < keymap_max; key++) {
                auto ctrl_block = ControlBlockIndex(key);
                auto value_block_idx_start = *((uint32_t*) (m_keymap + ctrl_block));
                auto value_block_idx_end = *((uint32_t*) (m_keymap + ctrl_block + m_keys_per_ctrl_block + ctrl_block_cell_size));
                auto value_block_size = value_block_idx_end - value_block_idx_start;

                std::sort(m_map+value_block_idx_start,m_map+value_block_idx_end);
            }
        }

// ctrl_block_keys_index
//  │                       Values array
//  │   │          │       │             │
//  │   ├──────────┤       ├─────────────┤
//  └──►│CTRL-BLOCK├──────►│             │<- ctrl_block_values_begin
//      ├──────────┤       │             │
//      │Key1      │       │             │
//      ├──────────┤       │             │
//      │Key2      │       │             │
//      ├──────────┤       │             │
//      │Key3      │       │             │
//      ├──────────┤       │             │    m_keys_per_ctrl_block
//      │Key4      │       │             │
//      ├──────────┤       │             │    #keys in this range
//      │Key5      │       │             │
//      ├──────────┤       │             │
//      │Key6      │       │             │
//      ├──────────┤       │             │
//      │Key7      │       │             │
//      ├──────────┤       │             │
//      │Key8      │       │             │
//      ├──────────┤       │             │
//      │CTRL-Block├───┐   │             │
//      ├──────────┤   │   │             │
//      │          │   │   ├─────────────┤
//      │          │   └──►│             │<- ctrl_block_values_end
//      ├──────────┤       │             │
//      │          │

        void PrintBlock(size_t key) {
            std::cout << "Print block key: " << KmerUtils::ToString(key, m_main_bits) << std::endl;
            auto ctrl_block_keys_index = ControlBlockIndex(key);
            auto ctrl_block_values_begin = *((uint32_t*)(m_keymap + ctrl_block_keys_index));
            auto ctrl_block_values_end = *((uint32_t*)(m_keymap + ctrl_block_keys_index + m_keys_per_ctrl_block + ctrl_block_cell_size));
            auto ctrl_block_values_size = ctrl_block_values_end - ctrl_block_values_begin;

            std::cout << "\n#######################################\n## KEY ARRAY ###########################~\n" << std::endl;
            std::cout << "-- Ctrl Block --- Keys: " << ctrl_block_keys_index << " -- Values: " << ctrl_block_values_begin << std::endl;
            std::cout << "Idx in values: " << ctrl_block_values_begin << " - " << ctrl_block_values_end << std::endl;
            std::cout << "----------------" << std::endl;
            for (int i = ctrl_block_cell_size; i < ctrl_block_cell_size + m_keys_per_ctrl_block; i++) {
                std::cout << (i-ctrl_block_cell_size) << ": " << (uint32_t) m_keymap[ctrl_block_keys_index + i] << ", ";
            }
            std::cout << std::endl;

            std::cout << "-- Next Ctrl Block --- Keys: " << ctrl_block_keys_index + m_keys_per_ctrl_block + ctrl_block_cell_size << " -- Values: " << ctrl_block_values_end << std::endl;
            std::cout << "   .... " << std::endl;
            std::cout << "\n## VALUE ARRAY ###########################~" << std::endl;

            Utils::Input();

            for (int j = ctrl_block_cell_size; j < ctrl_block_cell_size + m_keys_per_ctrl_block; j++) {
                size_t key_index = ctrl_block_keys_index + j;
                auto key_value_start = ctrl_block_values_begin + m_keymap[key_index];
                auto key_value_end = j == ctrl_block_cell_size + m_keys_per_ctrl_block - 1 ?
                        ctrl_block_values_end : ctrl_block_values_begin + m_keymap[key_index + 1];
                auto key_value_block_size = key_value_end - key_value_start;
                auto i = 0;
                std::cout << "-- Value index: " << key_value_start << " -------------" << (j - ctrl_block_cell_size) << " From, To: " << key_value_start << " - " << key_value_end << " (Size: " << FlexBlockSize(key_value_block_size) << ")      ";
                if (key_value_block_size >= m_flex_threshold) {
                    std::cout << " Flexi-K block (" << key_value_start << ")" << std::endl;
                    auto flex_block_size = FlexBlockSize(key_value_block_size);
                    for (; i < flex_block_size; i++) {
                        auto* cell = m_map + key_value_start + i;
                        std::cout << 2*i << ": " << KmerUtils::ToString(((uint32_t*)cell)[0], m_flex_k_bits) << "  " << (2*i + 1) << ": " << KmerUtils::ToString(((uint32_t*)cell)[1], m_flex_k_bits) << std::endl;
                    }
                    std::cout << "-------------";
                }
                std::cout << " Value block" << std::endl;
                for (; i < key_value_block_size; i++) {
                    std::cout << key_value_start + i << ": " << m_map[key_value_start + i].ToString() << std::endl;
                }
            }
            std::cout << "-------------";
            std::cout << "End printing block \n##########################" << std::endl;
        }

        void PrintKeyMap(size_t print_n, size_t wrap_around, bool skip_empty=true) {
            size_t row_index = 0;
            std::string line = "";
            size_t non_empty_block_count = 0;
            for (auto row = 0; row_index < keymap_size_total; row++, row_index = row * wrap_around) {
                bool line_empty = true;
                // check if there is a non empty cell
                for (auto col = 0; col < wrap_around && row_index + col < keymap_size_total; col++) {
                    if (m_keymap[row_index + col]) {
//                        std::cout << "line not empty(" << row + col << ", " << row << ", " << col << ") " << (uint32_t) keymap[row_index + col] << std::endl;
                        line_empty = false;
                        break;
                    }
                }
                // print if there is non-empty cell
                if (line_empty) continue;
                non_empty_block_count++;
                if (non_empty_block_count > print_n) return;

                std::cout << row_index << " ";
                for (auto col = 0; col < wrap_around; col++) {
                    std::cout << SeedmapUtils::BitString<8>(m_keymap[row_index + col]) << " ";
                }
                std::cout << '\n';
            }
        }
    };
}