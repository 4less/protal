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

namespace protal {
    template<uint64_t taxid_bits, uint64_t geneid_bits, uint64_t genepos_bits>
    struct Entry {
        uint64_t value;

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



    class Seedmap {
    private:
        using KeyMap_t = uint8_t;
        KeyMap_t* m_keymap = nullptr;
        uint64_t m_keymask = 0b0000000000000000000000000000000000111111111111111111111111111111;
        size_t keymap_size = 1u << 30u;
        size_t keymap_max = m_keymask;


        // This variable defines how many keys are managed by one control block
        // Must be power of two
        size_t ctrl_block_frequency = 16; // Version 1
//        size_t ctrl_block_frequency = 8; // Version 2

        // bitshift to find which control block a key is in.
        size_t ctrl_block_frequency_bitshift = log2(ctrl_block_frequency);

        uint64_t ctrl_block_key_mask = (1 << ctrl_block_frequency_bitshift) - 1;



        size_t max_key_ubiquity = (1 << 8) / ctrl_block_frequency;

//        size_t ctrl_block_frequency_bitshift = 4;

        // Bytespace a controlblock takes up
        size_t ctrl_block_size = 4; // Byte
        size_t ctrl_block_size_shift = log2(ctrl_block_size);
//        size_t ctrl_block_size_shift = 2;

        size_t keymap_size_total = keymap_size + (((keymap_size/ctrl_block_frequency)+1) * ctrl_block_size);

        /*
         * The number of total value entries per block.
         * Each block holds indices for <ctrl_block_frequency> keys
         * and the values for thes keys may not be more than <max_block_size>
         */
        size_t max_block_size = (1 << sizeof(KeyMap_t)*8);

        size_t values_size= 0;
        Entry<20,20,20>* m_map = nullptr;

        uint64_t m_found_counter = 0;

        void Print() {
            std::cout << "ctrl_block_frequency:             " << ctrl_block_frequency << std::endl;
            std::cout << "ctrl_block_frequency_bitshift:    " << ctrl_block_frequency_bitshift << std::endl;
            std::cout << "ctrl_block_size:                  " << ctrl_block_size << std::endl;
            std::cout << "ctrl_block_size_shift:            " << ctrl_block_size_shift << std::endl;
            std::cout << "ctrl_block_key_mask:              " << ctrl_block_key_mask << std::endl;
            std::cout << "keymap_size_total:                " << keymap_size_total << std::endl;
            std::cout << "max_block_size:                   " << max_block_size << std::endl;
            std::cout << "values_size:                      " << values_size << std::endl;
        }

        Seedmap(Seedmap const& other) = delete;
    public:
        Seedmap() {
            Print();
            m_keymap = new uint8_t[keymap_size_total];
        }

        Seedmap(std::string file) {
            Load(file);
        };

        ~Seedmap() {
            delete[] m_keymap;
            delete[] m_map;
        }

        inline uint64_t KeymapIndex(uint64_t key) const {
            return key + ((key >> ctrl_block_frequency_bitshift) << ctrl_block_size_shift) + ctrl_block_size;
        }

        inline uint64_t ControlBlockIndex(uint64_t key) const {
            return ((key >> ctrl_block_frequency_bitshift) << ctrl_block_frequency_bitshift) + ((key >> ctrl_block_frequency_bitshift) << (ctrl_block_size_shift));
        }

        inline uint64_t BlockKey(uint64_t key) const {
            return key & ctrl_block_key_mask;
        }

        void Save(std::string file) {

//            SortForKeys();

            std::ofstream ofs(file, std::ios::binary);
            ofs.write((char *) &keymap_size, sizeof(keymap_size));
            ofs.write((char *) &keymap_size_total, sizeof(keymap_size_total));
            ofs.write((char *) &ctrl_block_size, sizeof(ctrl_block_size));
            ofs.write((char *) &ctrl_block_frequency, sizeof(ctrl_block_frequency));
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
            ofs.write((char *) &ctrl_block_size, sizeof(ctrl_block_size));
            ofs.write((char *) &ctrl_block_frequency, sizeof(ctrl_block_frequency));
            ofs.write((char *) &ctrl_block_frequency_bitshift, sizeof(ctrl_block_frequency_bitshift));

            ofs.write((char *) &values_size, sizeof(values_size));

            ofs.write((char *) m_keymap, sizeof(*m_keymap) * keymap_size_total);
            ofs.write((char *) m_map, sizeof(*m_map) * values_size);
        }

        void Load(std::istream &ifs) {
            ifs.read((char *) &keymap_size, sizeof(keymap_size));
            ifs.read((char *) &keymap_size_total, sizeof(keymap_size_total));
            ifs.read((char *) &ctrl_block_size, sizeof(ctrl_block_size));
            ifs.read((char *) &ctrl_block_frequency, sizeof(ctrl_block_frequency));
            ifs.read((char *) &ctrl_block_frequency_bitshift, sizeof(ctrl_block_frequency_bitshift));
            ifs.read((char *) &values_size, sizeof(values_size));

            delete[] m_keymap;
            m_keymap = new uint8_t[keymap_size_total];
            ifs.read((char *) m_keymap, sizeof(*m_keymap) * (keymap_size_total));

            m_map = new Entry<20,20,20>[values_size];
            ifs.read((char *) m_map, sizeof(*m_map) * (values_size));
        }

        void Load(std::string file) {

            std::ifstream ifs(file, std::ios::binary);

            ifs.read((char *) &keymap_size, sizeof(keymap_size));
            ifs.read((char *) &keymap_size_total, sizeof(keymap_size_total));
            ifs.read((char *) &ctrl_block_size, sizeof(ctrl_block_size));
            ifs.read((char *) &ctrl_block_frequency, sizeof(ctrl_block_frequency));
            ifs.read((char *) &ctrl_block_frequency_bitshift, sizeof(ctrl_block_frequency_bitshift));
            ifs.read((char *) &values_size, sizeof(values_size));

            m_keymap = new uint8_t[keymap_size_total];
            ifs.read((char *) m_keymap, sizeof(*m_keymap) * (keymap_size_total));

            m_map = new Entry<20,20,20>[values_size];
            ifs.read((char *) m_map, sizeof(*m_map) * (values_size));
        }


        void CountUpKey(uint64_t key) {
            assert(key < keymap_size);
            auto index = KeymapIndex(key);
            m_keymap[index] += m_keymap[index] < UINT8_MAX;
        }

        void CountUpKeyOMP(uint64_t key) {
            assert(key < keymap_size);
            auto index = KeymapIndex(key);

#pragma omp atomic
            m_keymap[index] += m_keymap[index] < UINT8_MAX;
        }

        uint64_t GetKey(uint64_t key) {
            return m_keymap[KeymapIndex(key)];
        }

        uint8_t Key(uint64_t key) {
            return m_keymap[KeymapIndex(key)];
        }

        void Get(uint64_t key, Entry<20, 20, 20>* &start, Entry<20, 20, 20>* &end) {
            if (key > m_keymask) {
                std::cout << "Key is larger than keymask" << key << " > " << m_keymask << std::endl;
            }
            uint64_t block_start_idx = ControlBlockIndex(key);
            uint64_t block_end_idx = block_start_idx + ctrl_block_frequency + ctrl_block_size;

            uint64_t block_value_start_idx = *((uint32_t*) (m_keymap + block_start_idx) );
            uint64_t block_value_end_idx = *((uint32_t*) (m_keymap + block_end_idx) );
            uint64_t block_value_size = block_value_end_idx - block_value_start_idx;

            if (!block_value_size) {
                start = nullptr;
                end = nullptr;
                return;
            }

            size_t key_index = block_start_idx + ctrl_block_size + BlockKey(key);

            auto key_value_start = block_value_start_idx + m_keymap[key_index];
            auto key_value_end = (key_index + 1) == block_end_idx ? block_value_end_idx : block_value_start_idx + m_keymap[key_index + 1];
            auto key_value_block_size = key_value_end - key_value_start;

            if (!key_value_block_size) {
                start = nullptr;
                end = nullptr;
                return;
            }

            if (key_value_start > key_value_end) {
                PrintBlock(key);
                std::cout << "Seedmap: 240" << std::endl;
                exit(9);
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
            uint64_t block_end_idx = block_start_idx + ctrl_block_frequency + ctrl_block_size;

            uint64_t block_value_start_idx = *((uint32_t*) (m_keymap + block_start_idx ) );
            uint64_t block_value_end_idx = *((uint32_t*) (m_keymap + block_end_idx) );
            uint64_t block_value_size = block_value_end_idx - block_value_start_idx;

            if (!block_value_size) {
                return false;
            }

            size_t key_index = block_start_idx + ctrl_block_size + BlockKey(key);
            auto key_value_start = block_value_start_idx + m_keymap[key_index];
            auto key_value_end = (key_index + 1) == block_end_idx ? block_value_end_idx : block_value_start_idx + m_keymap[key_index + 1];
            auto key_value_block_size = key_value_end - key_value_start;

            if (!key_value_block_size) {
                return false;
            }

            auto index = key_value_start;

            if (key_value_start > key_value_end) {
                PrintBlock(key);
                std::cout << "Key value start (" << key_value_start << ") can not be larger than key value end (" << key_value_end << ")" << std::endl;
                exit(9);
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


        bool PutOMP(uint64_t &key, uint64_t &taxid, uint64_t &geneid, uint64_t &genepos) {
            if (m_map == nullptr) {
                std::cerr << "You need to initialize the values with BuildValuePointers()" << std::endl;
                exit(8);
            }
            if (key > keymap_max) {
                std::cerr << "Key too large (" << key << " > " << keymap_max << ")" << std::endl;
                exit(8);
            }


            uint64_t block_start_idx = ControlBlockIndex(key);
            uint64_t block_end_idx = block_start_idx + ctrl_block_frequency + ctrl_block_size;

            uint64_t block_value_start_idx = *((uint32_t*) (m_keymap + block_start_idx ) );
            uint64_t block_value_end_idx = *((uint32_t*) (m_keymap + block_end_idx) );
            uint64_t block_value_size = block_value_end_idx - block_value_start_idx;

            if (!block_value_size) {
                return false;
            }

            size_t key_index = block_start_idx + ctrl_block_size + BlockKey(key);
            auto key_value_start = block_value_start_idx + m_keymap[key_index];
            auto key_value_end = (key_index + 1) == block_end_idx ? block_value_end_idx : block_value_start_idx + m_keymap[key_index + 1];
            auto key_value_block_size = key_value_end - key_value_start;

            if (!key_value_block_size) {
                return false;
            }

            auto index = key_value_start;

            if (key_value_start > key_value_end) {
                PrintBlock(key);
                std::cout << "Seedmap: 363" << std::endl;
                exit(9);
            }

            for (; !m_map[index].Empty() && index < key_value_end; index++);

            if (index == key_value_end) {
                std::cout << "bucket too large, not present " << key_value_start << " size: " << (key_value_end - key_value_start) << std::endl;
                std::cout << "block_end_idx:  " << block_end_idx << " size: " << (key_value_end - key_value_start) << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;
                return false;
            }
#pragma omp critical(put)
{
            while(!m_map[index].Empty()) index++;
            m_map[index].Put(taxid, geneid, genepos);
}

            Entry<20,20,20> *begin, *end;
            Get(key, begin, end);

            return true;
        }


        void BuildValuePointers() {
            uint32_t control_idx = 0;
            uint32_t* control_ptr = (uint32_t*) m_keymap;

            uint64_t local_position = 0;

            uint32_t global_position = 0;

            uint64_t num_ctrl_blocks = keymap_size >> ctrl_block_frequency_bitshift;

            std::cout << "number ctrl blocks: " << num_ctrl_blocks << std::endl;

            uint64_t max_keyblock_size = max_key_ubiquity*2;


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

            subkey subkey[ctrl_block_frequency];
            auto subkey_begin = subkey;
            auto subkey_end = subkey + ctrl_block_frequency;

            // Understand kmer frequencies better to improve sensitivity of alignment.
            constexpr int kmer_freq_size = 1 << (sizeof(KeyMap_t) * 8);
            int kmer_frequencies[kmer_freq_size] = { 0 };
            size_t count_failed_demand = 0;
            size_t count_total_stored = 0;
            size_t count_total_demand = 0;


            for (auto block = 0; block < num_ctrl_blocks; block++) {
//                std::cout << "new block" << std::endl;
                auto block_index = block * (ctrl_block_frequency + ctrl_block_size);
                auto block_start_index = block_index + ctrl_block_size;


                control_ptr = (uint32_t*) (m_keymap + block_index);
                *control_ptr = global_position;

                // first iteration
                uint64_t total_count = 0;
                for (auto key_pos = 0; key_pos < ctrl_block_frequency; key_pos++) {
                    subkey[key_pos].key = key_pos;
                    subkey[key_pos].count = m_keymap[block_start_index + key_pos];

                    kmer_frequencies[subkey[key_pos].count]++;

                    total_count += subkey[key_pos].count;
                }
                if (!total_count) continue;
                // sort by count

                std::sort(subkey_begin, subkey_end, [](const struct subkey &a, const struct subkey &b) {
                    return a.count < b.count;
                });

                // second iteration
                uint64_t current_sum = 0;
                uint64_t total_demand = 0;

                for (auto key_pos = 0; key_pos < ctrl_block_frequency; key_pos++) {
                    total_demand += subkey[key_pos].count;
//                    subkey[key_pos].Print();
                    if (subkey[key_pos].count < max_keyblock_size && current_sum + subkey[key_pos].count < max_block_size) {
                        current_sum += subkey[key_pos].count;
                    } else {
                        subkey[key_pos].count = 0;
                    }
                }
                count_total_demand += total_demand;
                count_total_stored += current_sum;
                count_failed_demand += (current_sum < total_demand);



                // sort by key
//                std::cout << "----" << std::endl;
// WHY? Figure out
                std::sort(subkey_begin, subkey_end, [](const struct subkey &a, const struct subkey &b) {
                    return a.key < b.key;
                });

                local_position = 0;
                for (auto key_pos = 0; key_pos < ctrl_block_frequency; key_pos++) {
                    m_keymap[block_start_index + key_pos] = local_position;

                    local_position += subkey[key_pos].count;
                    global_position += subkey[key_pos].count;
//                    subkey[key_pos].Print();
                }
//                if (total_demand > 200) {
//                    std::cout << "current_sum: " << current_sum << "  total_demand: " << total_demand << std::endl;
//                    std::string stop;
//                    std::cin >> stop;
//                }
            }
//            exit(9);

//            for (auto key = 0; key < keymap_size; key++) {
//                if (key % 100'000'000 == 0) std::cout << key << std::endl;
//
//                auto key_index = KeymapIndex(key);
//                auto key_count = m_keymap[key_index];
//
//                auto current_control_idx = ControlBlockIndex(key);
//                if (current_control_idx != control_idx) {
//                    control_idx = current_control_idx;
//                    control_ptr = (uint32_t*) (m_keymap + control_idx);
//                    *control_ptr = global_position;
//
//                    local_position = 0;
//                }
//
//                m_keymap[key_index] = local_position;
//
//                if (key_count <= 16) {
//                    local_position += key_count;
//                    global_position += key_count;
//                }
//
//
//                if (current_control_idx >= keymap_size_total) {
//                    std::cout << current_control_idx << " >= " << keymap_size_total << std::endl;
//                    exit(12);
//                }
//            }

            // set last control pointer
            control_ptr = (uint32_t*) (m_keymap + keymap_size_total - 4);
            *control_ptr = global_position;
            values_size = global_position;

            std::cout << "space requirements" << std::endl;
            double key_mem = (double) keymap_size_total/(1024*1024*1024);
            double val_mem = (double) (values_size * 8)/(1024*1024*1024);
            std::cout << "keys:   " << key_mem << " GB" << std::endl;
            std::cout << "values: " << val_mem << " GB" << std::endl;
            std::cout << "total:  " << key_mem+val_mem << " GB" << std::endl;

//            exit(9);
            m_map = new Entry<20, 20, 20>[values_size];
//            std::string stop;
//            std::cin >> stop;

            for (auto i = 0; i < kmer_freq_size; i++) {
                std::cout << i << '\t' << kmer_frequencies[i] << std::endl;
            }
            std::cout << "Total stored: " << count_total_stored << std::endl;
            std::cout << "Total demand: " << count_total_demand << std::endl;
            std::cout << "Stored " << (100*static_cast<double>(count_total_stored)/static_cast<double>(count_total_demand)) << " of total values." << std::endl;
            std::cout << "Failed to meet demands? " << count_failed_demand << std::endl;
        }


        void SortForKeys() {
            std::cout << "Sort for keys" << std::endl;
            size_t row_index = 0;
            std::string line = "";
            size_t non_empty_block_count = 0;
            for (auto key = 0; key < keymap_max; key++) {
                auto ctrl_block = ControlBlockIndex(key);
                auto value_block_idx_start = *((uint32_t*) (m_keymap + ctrl_block));
                auto value_block_idx_end = *((uint32_t*) (m_keymap + ctrl_block + ctrl_block_frequency + ctrl_block_size));
                auto value_block_size = value_block_idx_end - value_block_idx_start;

                std::sort(m_map+value_block_idx_start,m_map+value_block_idx_end);
            }
        }

        void PrintBlock(size_t key) {
            auto ctrl_block = ControlBlockIndex(key);
            auto value_block_idx_start = *((uint32_t*) (m_keymap + ctrl_block));
            auto value_block_idx_end = *((uint32_t*) (m_keymap + ctrl_block + ctrl_block_frequency + ctrl_block_size));
            auto value_block_size = value_block_idx_end - value_block_idx_start;

            std::cout << "-- Ctrl Block (" << value_block_size << ") -- ctrl block idx : " << ctrl_block << std::endl;
            std::cout << value_block_idx_start << std::endl;
            std::cout << "----------------" << std::endl;
            for (int i = 4; i < 20; i++) {
                std::cout << (uint32_t) m_keymap[ctrl_block + i] << std::endl;
            }
            auto last_index = m_keymap[ctrl_block + 19];

            std::cout << "-- Ctrl Block (+1) -- ctrl block idx : " << ctrl_block << std::endl;
            std::cout << value_block_idx_end << std::endl;
            std::cout << "----------------" << std::endl;

            if (value_block_idx_end < value_block_idx_start) {
                std::cerr << value_block_idx_end << " < " << value_block_idx_start << std::endl;
            }
            if (value_block_idx_start + last_index > value_block_idx_end) {
                std::cerr << value_block_idx_start << " + " << (uint32_t) last_index << " > " << value_block_idx_end << std::endl;
                std::string stop;
                std::cin >> stop;
            }
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