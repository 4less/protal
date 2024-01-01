#include <iostream>
#include "RunProtal.h"
#include "Test.h"

//#include <htslib/sam.h>
//#include "Test.h"

int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);
    std::cout << "Test" << std::endl;
//
//    std::string first2 = "/media/fritsche/Extreme_SSD/data/reads/postdrama/test/SRR7280802_1c.fastq";
//    std::string second2 = "/media/fritsche/Extreme_SSD/data/reads/postdrama/test/SRR7280802_2c.fastq";
//    std::string first = "/media/fritsche/Extreme_SSD/data/reads/postdrama/test/sub1c.fastq";
//    std::string second = "/media/fritsche/Extreme_SSD/data/reads/postdrama/test/sub2c.fastq";
//
//    igzstream is1 { first.c_str() };
//    igzstream is2 { second.c_str() };
//    std::ifstream is1 { first.c_str() };
//    std::ifstream is2 { second.c_str() };
//    protal::SeqReaderPE reader{is1, is2};
//    FastxRecord record1;
//    FastxRecord record2;
//
//    for (std::string line; std::getline(is1, line);) {
//        std::cout << line << std::endl;
//    }
//
//    exit(9);
//    BufferedFastxReader readerx;

//    while (true) {
//        bool ok = false;
//
//        ok = readerx.LoadBatch(is1, 2048);
//        if (!ok) break;
//
//        // Read records from datablock
//        while (true) {
//            // Read from read1
//            auto valid_fragment = readerx.NextSequence(record1);
//            if (!valid_fragment) break;
//
//            std::cout << record1.to_string() << std::endl;
//        }
//    }
//    exit(12);
//
//    while (reader(record1, record2)) {
//
//        std::cout << record1.sequence.length() << " " << record1.quality.length() << std::endl;
//        std::cout << record2.sequence.length() << " " << record2.quality.length() << std::endl;
//
//        std::cout << "-----" << std::endl;
//        std::cout << record1.to_string() << std::endl;
//        std::cout << record2.to_string() << std::endl;
//    }
//    exit(90);

    // test::TestEndsFree();

    // test::SIMDSimilarity();
    protal::Run(argc, argv);

    return 0;
}
