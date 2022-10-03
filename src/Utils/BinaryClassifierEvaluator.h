//
// Created by fritsche on 03/02/2022.
//

#pragma once

#include <cstddef>
#include <ostream>
#include <iostream>

struct BinaryClassifierEvaluator {
    std::string name {};

    size_t tp = 0;
    size_t tn = 0;
    size_t fp = 0;
    size_t fn = 0;

    BinaryClassifierEvaluator() = default;
    explicit BinaryClassifierEvaluator(std::string name) :
        name(name) {};

    void SetName(std::string name) {
        this->name = name;
    }

    double Sensitivity() const {
        if (tp + fn == 0) return -1;
        return (double) tp/static_cast<double>(tp + fn);
    }

    double Specificity() const {
        if (tn + fp == 0) return -1;
        return (double) tn/static_cast<double>(tn + fp);
    }

    double Precision() const {
        if (tp + fp == 0) return -1;
        return (double) tp/static_cast<double>(tp + fp);
    }

    double Accuracy() const {
        if (tp + fp + tn + fn == 0) return -1;
        return (double) (tp + tn)/static_cast<double>(tp + tn + fp + fn);
    }

    double F1() const {
        if (tp + fp + tn + fn == 0) return -1;
        return (double) 2*tp/static_cast<double>(2 * tp + fp + fn);
    }


    void WriteStats(std::ostream &os=std::cout) const {
        os << "Evaluate binary classifier: " << std::endl;
        os << "tp: " << tp << " fp: " << fp << " tn: " << tn << " fn: " << fn << std::endl;
        os << "Sensitivity:\t\t" << Sensitivity() << std::endl;
        os << "Specificity:\t\t" << Specificity() << std::endl;
        os << "Precision:  \t\t" << Precision() << std::endl;
        os << "Accuracy:   \t\t" << Accuracy() << std::endl;
        os << "F1:         \t\t" << F1() << std::endl;
    }

    static void WriteRowHeader(std::ostream &os=std::cout, std::string name_col="Name") {
        os << name_col << "\tTP\tFP\tTN\tFN\tSensitivity\tSpecificity\tPrecision\tAccuracy\tF1" << std::endl;
    }

    void Join(BinaryClassifierEvaluator const& other) {
        this->tp += other.tp;
        this->tn += other.tn;
        this->fp += other.fp;
        this->fn += other.fn;
    }

    void WriteRowStats(std::ostream &os=std::cout) const {
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