#include <iostream>
#include "RunProtal.h"
#include <htslib/sam.h>
#include "Test.h"

int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);


//
//    std::string sam_path = "/usr/users/QIB_fr017/fritsche/Projects/data/protal/database/test/gprev/short.sam";
//
//    samFile *sam_in = sam_open(sam_path.c_str(), "r");
//    bam_hdr_t *header = sam_hdr_read(sam_in); //read header
//    bam1_t *aln = bam_init1();
//
//    std::cout << "Header: " << header << std::endl;
//
//    std::cout << header->ref_count << std::endl;
//    std::cout << header->n_targets << std::endl;
//    std::cout << header->text << std::endl;
//    std::cout << header->l_text << std::endl;
//    std::cout << sam_in->state << std::endl;
//    std::cout << sam_in->line.l << std::endl;
//
//    while (sam_read1(sam_in, header, aln) >= 0) {
//
//        int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
//        char *chr = header->target_name[aln->core.tid] ; //contig name (chromosome)
//        uint32_t len = aln->core.l_qseq; //length of the read.
//
//        uint8_t *q = bam_get_seq(aln); //quality string
//        uint32_t q2 = aln->core.qual ; //mapping quality
//
//
//        char *qseq = (char *)malloc(len);
//
//        std::cout << "Malloc: " << len << std::endl;
//        for(int i=0; i< len ; i++){
//            qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
//        }
//
//        std::cout << bam_get_qname(aln) << std::endl;
//        printf("%s\t%d\t%d\t%s\t%s\t%d\n",chr,pos,len,qseq,q,q2);
//
//    }
//
//    bam_destroy1(aln);
//    sam_close(sam_in);
//
//    exit(0);

//    test::TestEndsFree();
    protal::Run(argc, argv);

    return 0;
}
