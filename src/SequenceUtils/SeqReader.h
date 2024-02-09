//
// Created by fritsche on 14/08/22.
//

#pragma once

#include "FastxReader.h"
#include "omp.h"

namespace protal {

    class SeqReader {
    private:
//        const size_t m_block_size = (10 * 1024 * 1024);
        const size_t m_block_size = (1024);
        BufferedFastxReader m_reader;
        bool m_valid_fragment = false;
        bool m_valid_block = true;

        std::istream& m_is;

    public:
        SeqReader(std::istream& is) :
                m_is(is) {};

        SeqReader(SeqReader const& other) :
                m_is(other.m_is) {};

        inline void LoadBlockOMP(std::istream& is) {
#pragma omp critical(reader)
            {
                m_valid_block = m_reader.LoadBlock(is, m_block_size);
            }
        }

        bool operator() (FastxRecord &record) {
            m_valid_fragment = m_reader.NextSequence(record);
            if (m_valid_fragment) return true;

            LoadBlockOMP(m_is);
            if (!m_valid_block) return false;

            m_valid_fragment = m_reader.NextSequence(record);
            return m_valid_fragment;
        }
    };

    class SeqReaderPE {
        enum SeqReaderPEError {
            NO_ERROR,DIFF_LENGTH_SEQUENCE, DIFFERENT_LENGTH_BLOCK
        };
    private:
//        const size_t m_block_size = (10 * 1024 * 1024);
        const size_t m_block_size = (1024);
        const size_t m_record_count = 32;
        BufferedFastxReader m_reader_1;
        BufferedFastxReader m_reader_2;

        bool m_valid_fragment_1 = false;
        bool m_valid_fragment_2 = false;

        bool m_valid_block_1 = true;
        bool m_valid_block_2 = true;

        std::istream& m_is1;
        std::istream& m_is2;

        SeqReaderPEError m_error_code = NO_ERROR;
        bool m_success = true;

    public:
        SeqReaderPE(std::istream& is1, std::istream& is2) :
                m_is1(is1),
                m_is2(is2) {};

        SeqReaderPE(SeqReaderPE const& other) :
                m_is1(other.m_is1),
                m_is2(other.m_is2) {};


        auto& GetFirstStream() {
            return m_is1;
        }

        inline void LoadBlockOMP() {
#pragma omp critical(reader)
            {
//                m_valid_block_1 = m_reader_1.LoadBlock(m_is1, m_block_size);
//                m_valid_block_2 = m_reader_2.LoadBlock(m_is2, m_block_size);
                m_valid_block_1 = m_reader_1.LoadBatch(m_is1, m_record_count);
                m_valid_block_2 = m_reader_2.LoadBatch(m_is2, m_record_count);
            }
        }

        inline void LoadBlockOMP(bool first) {
#pragma omp critical(reader)
            {
                if (first)
                    m_valid_block_1 = m_reader_1.LoadBlock(m_is1, m_block_size);
                else
                    m_valid_block_2 = m_reader_2.LoadBlock(m_is2, m_block_size);
            }
        }

        bool Success() const {
            return m_success;
        }

        void UpdateSuccess(SeqReaderPE const& other) {
            m_success &= other.m_success;
        }

        bool operator() (FastxRecord &record1, FastxRecord &record2) {
            m_valid_fragment_1 = m_reader_1.NextSequence(record1);
            m_valid_fragment_2 = m_reader_2.NextSequence(record2);

            if (m_reader_1.Error() || m_reader_2.Error()) {
                m_success = false;
                return false;
            }

            if (m_valid_fragment_1 && m_valid_fragment_2) {
                return true;
            }

            // Paired end states must always be identical.
            if (m_valid_fragment_1 != m_valid_fragment_2) {
                std::cerr << omp_get_thread_num() << "Error with next sequence. Paired end file streams are not of the same length. Abort. (R1: " << m_valid_fragment_1 << ", R2: " << m_valid_fragment_2 << ")" << std::endl;
                m_error_code = DIFF_LENGTH_SEQUENCE;
                m_success = false;
                return false;
            }

            LoadBlockOMP();
            if (!m_valid_block_1 && !m_valid_block_2) {
                return false;
            }

            // Paired end states must always be identical.
            if (m_valid_block_1 != m_valid_block_2) {
                std::cerr << omp_get_thread_num() << " Error after load block. Paired end file streams are not of the same length. Abort." << std::endl;
                m_error_code = DIFFERENT_LENGTH_BLOCK;
                m_success = false;
                return false;
            }

            m_valid_fragment_1 = m_reader_1.NextSequence(record1);
            m_valid_fragment_2 = m_reader_2.NextSequence(record2);

            if (m_reader_1.Error() || m_reader_2.Error()) {
                m_success = false;
                return false;
            }

            // Paired end states must always be identical.
            if (m_valid_fragment_1 != m_valid_fragment_2) {
                std::cerr << "2 Error with next sequence. Paired end file streams are not of the same length. Abort." << std::endl;
                m_error_code = DIFF_LENGTH_SEQUENCE;
                m_success = false;
                return false;
            }

            return m_valid_fragment_1 && m_valid_fragment_2;
        }
    };

}
