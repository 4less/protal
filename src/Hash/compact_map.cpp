//
// Created by fritsche on 08/03/2021.
//

#include "compact_map.h"
#include "KmerUtils.h"
#include <tgmath.h>
#include <cstring>
#include <fstream>
#include "omp.h"
#include "assert.h"

#define USE_0MP 1

IndexedMap::IndexedMap(size_t key_bits, size_t value_bits, size_t offset_bits, size_t* offset_sizes, double load_factor, bool multimap) :
        size_(0),
        key_bits_(key_bits),
        value_bits_(value_bits),
        offset_bits_(offset_bits),
        buckets_count_(1llu << offset_bits_),
        load_factor_(load_factor),
        multimap_(multimap) {

    offset_ = new uint64_t[buckets_count_ + 1];

    size_t total = 0;
    for (auto i = 0; i < buckets_count_; i++) {
        offset_[i] = total;
        total += ComputeBucketCapacity(offset_sizes[i], load_factor_);
    }
    offset_[buckets_count_] = total;
    capacity_ = total;

    map_ = new Cell[capacity_];
    memset(map_, 0, capacity_);

    // analyse psls
    psl_ = new size_t[100];
    std::fill_n(psl_, 100, 0llu);
}


IndexedMap::IndexedMap(size_t key_bits, size_t value_bits, size_t offset_bits, std::string bucket_sizes, double load_factor, bool multimap) :
        size_(0),
        key_bits_(key_bits),
        value_bits_(value_bits),
        offset_bits_(offset_bits),
        buckets_count_(1llu << offset_bits_),
        load_factor_(load_factor),
        multimap_(multimap)  {

    LoadBucketSizes(bucket_sizes);

    psl_ = new size_t[100];
    std::fill_n(psl_, 100, 0llu);
}



Cell* IndexedMap::Find(uint64_t &key) {
    uint64_t bucket_key = key >> key_bits_;
    uint64_t internal_key = key & ((1llu << key_bits_) - 1);

    assert(key < (1llu << (key_bits_ + offset_bits_)));
    assert(bucket_key < buckets_count_);

    size_t bucket_index = offset_[bucket_key];
    uint32_t bucket_size = offset_[bucket_key + 1] - bucket_index;

    // turn off when checking the db against reference to see if all k-mers are contained
    if (bucket_size == 0) return nullptr;

    assert(bucket_size != 0);

    size_t cur_pos = CompactMapUtils::Hash(internal_key, bucket_size);
    Cell* cell = map_ + bucket_index + cur_pos;
    auto probes = 0;

    for (int i = 0; i < 1; i++) {
        if (probes == bucket_size) return nullptr;
        if (cell->key(value_bits_) == internal_key && !cell->empty()) {
            return cell;
        }

        ++cell;
//        cell -= (++cur_pos == bucket_size) * bucket_size;
//        cur_pos -= (cur_pos == bucket_size) * bucket_size;
        if (++cur_pos == bucket_size) {
            cell = map_ + bucket_index;
            cur_pos = 0;
        }
        ++probes;
    }

    size_t c_psl = psl(cell->key(value_bits_), cur_pos, bucket_size);

    while (probes <= c_psl && !cell->empty()) {
        if (cell->key(value_bits_) == internal_key) {
//            psl_[c_psl]++;
            return cell;
        }

        ++cell;
        ++probes;
        if (++cur_pos == bucket_size) {
            cell = map_ + bucket_index;
            cur_pos = 0;
        }
        c_psl = psl(cell->key(value_bits_), cur_pos, bucket_size);
    }

    return nullptr;
}

size_t IndexedMap::CalcSize() {
    size_t size = 0;
    for (auto& cell : *this) {
        if (!cell.empty()) {
            size++;
        }
    }
    return size;
}

void IndexedMap::Insert(uint64_t key, uint64_t value) {

    // Within 0 < key < pre_key_bits
    // Key that is being implicitly stored in offset_
    uint64_t bucket_key = key >> key_bits_;

    // Key that is being explicitly stored in map_
    uint64_t main_key = key & ((1llu << key_bits_) - 1);

    size_t bucket_index = offset_[bucket_key];
    uint32_t bucket_size = offset_[bucket_key + 1] - bucket_index;

    assert(GetBucketSize(key) == bucket_size);
    assert(bucket_size != 0);
    assert(bucket_key < buckets_count_);
    assert(key < (1llu << (key_bits_ + offset_bits_)));

    // new ntry psl and resident entry psl
    size_t cur_pos = CompactMapUtils::Hash(main_key, bucket_size);
    size_t ins_psl = 0;
    size_t res_psl = 0;

    // navigate to bucket in map_
    Cell* cell = map_ + bucket_index + cur_pos;


    // probing until found
    while (!cell->empty()) {
        if (cell->key(value_bits_) == main_key) {
            size_--;
            break;
        }

//        #if DEBUG
//        if (ins_psl >= bucket_size) {
//            print(bucket_index, bucket_size);
//
//            std::cout << "key: " << key << std::endl;
//            std::cout << "bucket_key: " << bucket_key << std::endl;
//            std::cout << "bucket_size: " << bucket_size << std::endl;
//            std::cout << "main_key: " << main_key << std::endl;
//            std::cout << "Hash(main_key): " << CompactMapUtils::Hash(main_key, bucket_size) << std::endl;
//        }
//        #endif


        if (ins_psl >= bucket_size) {
            std::cout << "bucket_size: " << bucket_size << std::endl;
            std::cout << "insert psl: " << ins_psl << std::endl;
            std::cout << "bucket_key: " << bucket_key << std::endl;
            std::cout << "key: " << key << std::endl;
            std::cout << "mainkey: " << main_key << std::endl;
            std::cout << "keybits: " << key_bits_ << std::endl;
            print(bucket_index, bucket_size, true);
            exit(9);
        }

        assert(ins_psl < bucket_size);
        // TEST NEW FEATURE
        if (!multimap_ && cell->key(value_bits_) == main_key) {
            size_--;
            exit(9);
            break;
        }

        res_psl = psl(cell->key(value_bits_), cur_pos, bucket_size);

        if (ins_psl > res_psl) {
            // find
            // Swap entries by copying cell data into tmp key and value
            auto tmp_key = cell->key(value_bits_);
            auto tmp_value = cell->value(value_bits_);

            // update cell values with current key and value
            cell->paste(main_key, value, value_bits_);

            // continue with swapped key
            main_key = tmp_key;
            value = tmp_value;

            ins_psl = res_psl;
        }

        cur_pos++;
        ins_psl++;
        cell++;

        if (cur_pos == bucket_size) {
            cell = map_ + bucket_index;
            cur_pos = 0;
        }
    }

    size_++;

    if (!cell->empty()) {
        std::cerr << "cell should be empty" << std::endl;
        exit(9);
    }

    cell->paste(main_key, value, value_bits_);

//    if (CalcSize() != size_) {
//        std::cerr << "calc: " << CalcSize() << " vs " << size_ << std::endl;
//        exit(9);
//    }
}

using InsertKeyExists = std::function<size_t (const size_t, const size_t)>;
void IndexedMap::InsertPlus(uint64_t key, uint64_t value, InsertKeyExists &function) {
    // Within 0 < key < pre_key_bits
    // Key that is being implicitly stored in offset_
    uint64_t bucket_key = key >> key_bits_;

    // Key that is being explicitly stored in map_
    uint64_t main_key = key & ((1llu << key_bits_) - 1);

    size_t bucket_index = offset_[bucket_key];
    uint32_t bucket_size = offset_[bucket_key + 1] - bucket_index;

    assert(GetBucketSize(key) == bucket_size);

    assert(bucket_size != 0);
    assert(bucket_key < buckets_count_);
    assert(key < (1llu << (key_bits_ + offset_bits_)));

    // new ntry psl and resident entry psl
    size_t cur_pos = CompactMapUtils::Hash(main_key, bucket_size);
    size_t ins_psl = 0;
    size_t res_psl = 0;

    // navigate to bucket in map_
    Cell* cell = map_ + bucket_index + cur_pos;

    // probing until found
    while (!cell->empty()) {
#if DEBUG
        if (ins_psl >= bucket_size) {
            print(bucket_index, bucket_size);

            std::cout << "key: " << key << std::endl;
            std::cout << "bucket_key: " << bucket_key << std::endl;
            std::cout << "bucket_size: " << bucket_size << std::endl;
            std::cout << "main_key: " << main_key << std::endl;
            std::cout << "Hash(main_key): " << CompactMapUtils::Hash(main_key, bucket_size) << std::endl;
        }
#endif

        assert(ins_psl < bucket_size);
        if (cell->key(value_bits_) == main_key) {
//            std::cout << "Key exists" << std::endl;

            // This is the "plus" functionality and it is experimental
            cell->paste(main_key, function(cell->value(value_bits_), value), value_bits_);
            cell = nullptr;
            size_--; break;
        }

        res_psl = psl(cell->key(value_bits_), cur_pos, bucket_size);

        if (ins_psl > res_psl) {
            // find
            // Swap entries by copying cell data into tmp key and value
            auto tmp_key = cell->key(value_bits_);
            auto tmp_value = cell->value(value_bits_);

            // update cell values with current key and value
            cell->paste(main_key, value, value_bits_);

            // continue with swapped key
            main_key = tmp_key;
            value = tmp_value;

            ins_psl = res_psl;
        }

        cur_pos++;
        ins_psl++;
        cell++;

        if (cur_pos == bucket_size) {
            cell = map_ + bucket_index;
            cur_pos = 0;
        }
    }

    size_++;

    if (cell)
        cell->paste(main_key, value, value_bits_);
}

void IndexedMap::InsertSafe(uint64_t key, uint64_t value) {
    if (key > (1llu << (offset_bits_ + key_bits_)))
        errx(EX_SOFTWARE, "key %u must be smaller than offset_bits_ %u + key_size %llu = %llu", (unsigned long long) key, (unsigned int) offset_bits_, (unsigned long long) (1 < key), (unsigned long long) (1llu << (offset_bits_ + key_bits_)));
    if (value > (1llu << value_bits_))
        errx(EX_SOFTWARE, "value %u must be smaller than value_bits %u", (unsigned int) value, (unsigned int) value_bits_);

    Insert(key, value);
}

void IndexedMap::Save(std::string file) {
    ofstream ofs(file, ofstream::binary);
    ofs.write((char *) &capacity_, sizeof(capacity_));
    ofs.write((char *) &size_, sizeof(size_));
    ofs.write((char *) &load_factor_, sizeof(load_factor_));

    ofs.write((char *) &key_bits_, sizeof(key_bits_));
    ofs.write((char *) &value_bits_, sizeof(value_bits_));

    ofs.write((char *) &offset_bits_, sizeof(offset_bits_));
    ofs.write((char *) &buckets_count_, sizeof(buckets_count_));

    ofs.write((char *) offset_, sizeof(offset_) * (buckets_count_ + 1));
    ofs.write((char *) map_, sizeof(map_) * capacity_);
    ofs.close();
}

IndexedMap* IndexedMap::Load(std::string file, bool silent) {
    size_t capacity, size, offset_bits, offset_size, key_bits, value_bits;
    double load_factor;

    std::ifstream ifs(file, std::ifstream::binary);

    ifs.read((char *) &capacity, 8);
    ifs.read((char *) &size, 8);
    ifs.read((char *) &load_factor, 8);

    ifs.read((char *) &key_bits, 8);
    ifs.read((char *) &value_bits, 8);

    ifs.read((char *) &offset_bits, 8);
    ifs.read((char *) &offset_size, 8);

    if (!silent) {
        std::cout << "capacity:   \t" << capacity << std::endl;
        std::cout << "size:       \t" << size << std::endl;
        std::cout << "load_factor:\t" << load_factor << std::endl;
        std::cout << "key_bits:   \t" << key_bits << std::endl;
        std::cout << "value_bits: \t" << value_bits << std::endl;
        std::cout << "offset_bits:\t" << offset_bits << std::endl;
        std::cout << "offset_size:\t" << offset_size << std::endl;
    }


    uint64_t* offsets = new uint64_t[offset_size + 1];
    ifs.read((char *) offsets, sizeof(offsets) * (offset_size + 1));

//    for (int i = 0; i < offset_size+1; i++) {
//        std::cout << i << ": " << offsets[i] << std::endl;
//    }

    Cell* map;
    try {
        map = new Cell[capacity];
    } catch (std::bad_alloc &ex) {
        std::cerr << "Failed attempt to allocate " << (sizeof(*map) * capacity) << "bytes;\n"
                  << "you may not have enough free memory to load this database.\n"
                  << "If your computer has enough RAM, perhaps reducing memory usage from\n"
                  << "other programs could help you load this database?" << std::endl;
        errx(EX_OSERR, "unable to allocate hash table memory");
    }
    ifs.read((char *) map, sizeof(map) * capacity);

//    for (int i = 0; i < capacity; i++) {
//        std::cout << map[i].data << std::endl;
//    }

    if (!std::char_traits<char>::eof() == ifs.get())
        errx(EX_IOERR, "Failed to load database file %s", file.c_str());

    IndexedMap* imap = new IndexedMap(capacity, size, key_bits, value_bits, load_factor, offset_bits, offset_size, offsets, map);
    return imap;
}

IndexedMap::IndexedMap(size_t capacity, size_t size, const size_t key_bits,
                       const size_t value_bits, const double load_factor, const size_t offset_bits,
                       const size_t offset_size, uint64_t *offset, Cell *map) : capacity_(capacity),
                                                                                size_(size),
                                                                                key_bits_(key_bits),
                                                                                value_bits_(value_bits),
                                                                                load_factor_(load_factor),
                                                                                offset_bits_(offset_bits),
                                                                                buckets_count_(offset_size),
                                                                                offset_(offset), map_(map) {}

IndexedMap::IndexedMap() {}

Cell *IndexedMap::Map() {
    return map_;
}
