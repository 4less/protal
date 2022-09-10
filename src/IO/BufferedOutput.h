//
// Created by fritsche on 28/11/2021.
//

#ifndef VARKIT_BUFFEREDOUTPUT_H
#define VARKIT_BUFFEREDOUTPUT_H

#include <bits/stdc++.h>
#include <fmt/core.h>
#include <fmt/os.h>

using namespace std;

enum BufferedOutputUnit {
    N, MB, GB
};

template<typename T>
class BufferedOutput {
    T* buffer_ = nullptr;

    size_t buffer_fill_ = 0;

    size_t buffer_size_;
    static constexpr size_t entry_size_ = sizeof(T);


public:
    BufferedOutput(size_t buffer_size) {
        buffer_size_ = buffer_size;
        buffer_ = new T[buffer_size_];
    }

    BufferedOutput(size_t buffer_size, BufferedOutputUnit bsu) {
        switch (bsu) {
            case N:
                buffer_size_ = buffer_size;
                break;
            case MB:
                buffer_size_ = (buffer_size*1024*1024) / sizeof(T);
                break;
            case GB:
                buffer_size_ = (buffer_size*1024*1024*1024) / sizeof(T);
                break;
        }
//        std::cout << "Allocate with: " << buffer_size_ << " -> " << (sizeof(T) * buffer_size_)/(1024*1024) << " MB" << std::endl;

        buffer_ = new T[buffer_size_];
    }

    ~BufferedOutput() {
        delete[] buffer_;
    }

    inline bool Write(T &item) {
        assert(buffer_fill_ + entry_size_ < buffer_size_);
        memcpy(buffer_ + buffer_fill_, &item, entry_size_);
        buffer_fill_++;
        return buffer_fill_ < buffer_size_;
    }

    inline void Write(std::ostream &ofs) {
        if (buffer_fill_ != 0) {
//            std::cout << "write snps" << buffer_fill_ << std::endl;

            ofs.write(reinterpret_cast<char *>(buffer_), buffer_fill_ * sizeof(T));
            buffer_fill_ = 0;
        }
    }
};

class BufferedStringOutput {
    char *m_buffer = nullptr;

    constexpr static double DEFAULT_SOFT_CAPACITY_RATIO = 0.8;

    size_t m_buffer_size = 0;

    size_t m_buffer_capacity;
    size_t m_buffer_soft_capacity;
    double m_buffer_soft_capacity_ratio;

    void EmergencyIncreaseBuffer(size_t by) {
        std::cout << "EmergencyIncreaseBuffer from " << m_buffer_capacity;
        m_buffer_capacity += by;
        std::cout << " to " << m_buffer_capacity << std::endl;
        char* new_buffer = new char[m_buffer_capacity];
        std::memcpy(new_buffer, m_buffer, m_buffer_size);
        delete[] m_buffer;
        m_buffer = new_buffer;
        new_buffer = nullptr;

        exit(91);
    }

public:
    explicit BufferedStringOutput(size_t capacity, double soft_capacity_ratio=DEFAULT_SOFT_CAPACITY_RATIO) :
            m_buffer_capacity(capacity),
            m_buffer_soft_capacity_ratio(soft_capacity_ratio),
            m_buffer_soft_capacity(soft_capacity_ratio * capacity) {
        m_buffer = new char[m_buffer_capacity];
    }

    BufferedStringOutput(BufferedStringOutput const& other) {
        std::cout << "copy constructor" << std::endl;
        exit(9);
    }

    ~BufferedStringOutput() {
        std::cout << "Destructor bufferedStringOutput" << std::endl;
        delete[] m_buffer;
        std::cout << "Destructor bufferedStringOutput after" << std::endl;
    }

    inline BufferedStringOutput& operator << (std::string const&& str) {
        if (m_buffer_size + str.length() > m_buffer_capacity) {
            EmergencyIncreaseBuffer(str.length());
        }
        std::memcpy(m_buffer + m_buffer_size, str.c_str(), str.length());
        m_buffer_size += str.length();
        return *this;
    }

    inline BufferedStringOutput& operator << (std::string const& str) {
        if (m_buffer_size + str.length() > m_buffer_capacity) {
            EmergencyIncreaseBuffer(str.length());
        }
        std::memcpy(m_buffer + m_buffer_size, str.c_str(), str.length());
        m_buffer_size += str.length();
        return *this;
    }

    inline size_t Capacity() const {
        return m_buffer_capacity;
    }

    [[nodiscard]] bool Full() const {
        return m_buffer_size >= m_buffer_soft_capacity;
    }

    inline void Write(std::ostream &os) {
        if (m_buffer_size != 0) {
            os.write(reinterpret_cast<char *>(m_buffer), m_buffer_size);
            m_buffer_size = 0;
        }
    }
};





#endif //VARKIT_BUFFEREDOUTPUT_H
