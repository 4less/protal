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
    explicit BufferedOutput(size_t buffer_size) {
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
    std::string m_line_buffer{};

    size_t m_buffer_size = 0;
    size_t m_buffer_capacity;

public:
    explicit BufferedStringOutput(size_t capacity) :
            m_buffer_capacity(capacity) {
        m_buffer = new char[m_buffer_capacity];
    }

    BufferedStringOutput(size_t buffer_capacity, BufferedOutputUnit bsu) {
        switch (bsu) {
            case N:
                m_buffer_capacity = buffer_capacity;
                break;
            case MB:
                m_buffer_capacity = (buffer_capacity*1024lu*1024lu);
                break;
            case GB:
                m_buffer_capacity = (buffer_capacity*1024lu*1024lu*1024lu);
                break;
        }
        m_buffer = new char[m_buffer_capacity];
    }

    ~BufferedStringOutput() {
        delete[] m_buffer;
    }

    inline size_t Capacity() const {
        return m_buffer_capacity;
    }

    inline void Write(std::ostream &os) {
        if (m_buffer_size != 0) {
            os.write(reinterpret_cast<char *>(m_buffer), m_buffer_size);
            os.write(m_line_buffer.c_str(), m_line_buffer.length());
            m_line_buffer.clear();
            m_buffer_size = 0;
        }
    }

    inline bool Write(std::string& line) {
        if (m_buffer_size + line.length() > m_buffer_capacity) {
            m_line_buffer = line;
            return false;
        }
        std::memcpy(m_buffer + m_buffer_size, line.c_str(), line.length());
        m_buffer_size += line.length();
    }
};





#endif //VARKIT_BUFFEREDOUTPUT_H
