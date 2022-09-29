//
// Created by fritsche on 03/02/2022.
//

#pragma once

#include <cstddef>
#include <ostream>
#include <iostream>

struct BinaryClassifierEvaluator {
    std::string name = "";

    size_t tp = 0;
    size_t tn = 0;
    size_t fp = 0;
    size_t fn = 0;

    void SetName(std::string name) {
        this->name = name;
    }

    double Sensitivity() {
        if (tp + fn == 0) return -1;
        return (double) tp/(tp + fn);
    }

    double Specificity() {
        if (tn + fp == 0) return -1;
        return (double) tn/(tn + fp);
    }

    double Precision() {
        if (tp + fp == 0) return -1;
        return (double) tp/(tp + fp);
    }

    double Accuracy() {
        if (tp + fp + tn + fn == 0) return -1;
        return (double) (tp + tn)/(tp + tn + fp + fn);
    }

    double F1() {
        if (tp + fp + tn + fn == 0) return -1;
        return (double) 2*tp/(2 * tp + fp + fn);
    }


    void WriteStats(std::ostream &os=std::cout) {
        os << "Evaluate binary classifier: " << std::endl;
        os << "tp: " << tp << " fp: " << fp << " tn: " << tn << " fn: " << fn << std::endl;
        os << "Sensitivity:\t\t" << Sensitivity() << std::endl;
        os << "Specificity:\t\t" << Specificity() << std::endl;
        os << "Precision:  \t\t" << Precision() << std::endl;
        os << "Accuracy:   \t\t" << Accuracy() << std::endl;
        os << "F1:         \t\t" << F1() << std::endl;
    }

    void WriteRowHeader(std::ostream &os=std::cout, std::string name_col="Name")  {
        os << name_col << "\tTP\tFP\tTN\tFN\tSensitivity\tSpecificity\tPrecision\tAccuracy\tF1" << std::endl;
    }

    void Join(BinaryClassifierEvaluator const& other) {
        this->tp += other.tp;
        this->tn += other.tn;
        this->fp += other.fp;
        this->fn += other.fn;
    }

    void WriteRowStats(std::ostream &os=std::cout) {
        os << name << '\t';
        os << tp << '\t';
        os << fp << '\t';
        os << tn << '\t';
        os << fn << '\t';
        os << Sensitivity() << '\t';
        os << Specificity() << '\t';
        os << Precision() << '\t';
        os << Accuracy() << '\t';
        os << F1() << std::endl;
    }
};