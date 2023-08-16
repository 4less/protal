//
// Created by fritsche on 28/11/2021.
//

#ifndef VARKIT_BUFFEREDOUTPUT_H
#define VARKIT_BUFFEREDOUTPUT_H

#include <bits/stdc++.h>
//#include <fmt/core.h>
//#include <fmt/os.h>

using namespace std;

enum BufferedOutputUnit {
    N, MB, GB
};

template<typename T>
class BufferedOutput {
    T* m_buffer = nullptr;

    size_t m_buffer_size = 0;

    size_t m_buffer_capacity;
    static constexpr size_t m_entry_size = sizeof(T);


public:
    explicit BufferedOutput(size_t buffer_size) {
        m_buffer_capacity = buffer_size;
        m_buffer = new T[m_buffer_capacity];
    }

    inline size_t Capacity() const {
        return m_buffer_capacity;
    }

    BufferedOutput(size_t buffer_size, BufferedOutputUnit bsu) {
        switch (bsu) {
            case N:
                m_buffer_capacity = buffer_size;
                break;
            case MB:
                m_buffer_capacity = (buffer_size * 1024 * 1024) / sizeof(T);
                break;
            case GB:
                m_buffer_capacity = (buffer_size * 1024 * 1024 * 1024) / sizeof(T);
                break;
        }
//        std::cout << "Allocate with: " << m_buffer_capacity << " -> " << (sizeof(T) * m_buffer_capacity)/(1024*1024) << " MB" << std::endl;

        m_buffer = new T[m_buffer_capacity];
    }

    ~BufferedOutput() {
        delete[] m_buffer;
    }

    inline bool Write(T &item) {
        if (m_buffer_size >= m_buffer_capacity) {
            std::cout << "m_buffer_size " << m_buffer_size << std::endl;
            std::cout << "m_entry_size " << m_entry_size << std::endl;
            std::cout << "m_buffer_capacity " << m_buffer_capacity << std::endl;
        }
        assert(m_buffer_size < m_buffer_capacity);
        memcpy(m_buffer + m_buffer_size, &item, m_entry_size);
        m_buffer_size++;
        return m_buffer_size < m_buffer_capacity;
    }

    inline void Write(std::ostream &ofs) {
        if (m_buffer_size != 0) {
//            std::cout << "write snps" << m_buffer_size << std::endl;

            ofs.write(reinterpret_cast<char *>(m_buffer), m_buffer_size * sizeof(T));
            m_buffer_size = 0;
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
        return true;
    }

    inline bool Write(std::string&& line) {
        if (m_buffer_size + line.length() + 1 > m_buffer_capacity) {
            m_line_buffer = line + '\n';
            return false;
        }
        std::memcpy(m_buffer + m_buffer_size, line.c_str(), line.length());
        m_buffer_size += line.length();
        m_buffer[m_buffer_size++] = '\n';
        return true;
    }
};





#endif //VARKIT_BUFFEREDOUTPUT_H
