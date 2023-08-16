//
// Created by fritsche on 21/02/23.
//

#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <numeric>


enum VariantType {
    INS, DEL, SNP
};

using VariantID = size_t;
using Orientation = bool;
using Base = char;
using FirstRead = bool;
using ReadId = size_t;
using SSize = uint16_t;
using Qual = uint8_t;
using QualList = std::vector<Qual>;
using VariantPos = uint32_t;

class Variant {

    VariantID variant_id = 0;   // size_t
    VariantType variant_type;   // enum
    VariantPos position = UINT32_MAX;        // uint32_t
    Base reference = 'X';             // char
    Base variant = 'X';               // char
    SSize structural_size = 1;  // uint16_t
    uint16_t observations_fwd = 0;
    uint16_t observations_rev = 0;
    bool is_valid = true;
    bool is_major = false;
    std::string* structural = nullptr;
    QualList quals;

public:
    Variant(VariantPos pos, Base snp, Base ref) :
            variant_type(VariantType::SNP), position(pos), variant(snp), reference(ref) {};

    Variant(VariantType type, VariantPos pos, Base ref, std::string& structural) :
            variant_type(type), position(pos), reference(ref), structural(new std::string(structural)),
            structural_size(structural.length()) {};

    size_t Observations() const {
        return observations_fwd + observations_rev;
    }

    auto& GetQualList() {
        return quals;
    }
    auto GetQualListCopy() const {
        return quals;
    }

    auto IsUninitialized() const {
        return variant == 'X' && reference == 'X';
    }

    size_t QualitySum() const {
        if (IsUninitialized()) {
            exit(12);
        }
        if (IsReference()) return Observations() * 40;
        return std::accumulate(quals.begin(), quals.end(), 0, [](size_t acc, const uint16_t q) {
            return acc + q;
        });
    }
    size_t MeanQuality() const {
        return static_cast<size_t>(static_cast<double>(QualitySum())/Observations());
    }


    void SetValid(bool is_valid) {
        this->is_valid = is_valid;
    }

    bool GetValid() const {
        return is_valid;
    }

    void SetMajorAllele(bool is_major) {
        this->is_major = is_major;
    }

    bool IsMajorAllele() const {
        return is_major;
    }

    ~Variant() {
        if (!structural) {
            delete[] structural;
        }
    }


    std::string ToString() const {
        std::string str;
        str += '{';
        str += std::to_string(is_valid) + ' ';
        if (variant_type == VariantType::SNP) {
            str += (IsReference() ? "REF " : "SNP ");
        } else if (variant_type == VariantType::INS) {
            str += "INS ";
        } else if (variant_type == VariantType::DEL) {
            str += "DEL ";
        }

        if (variant_type == VariantType::SNP) {
            str += (IsReference() ? std::string(1, reference) : std::string(1, reference) + "->" + std::string(1, variant)) + " ";
//            str += std::to_string(quality) + '\t';
        } else if (variant_type == VariantType::INS || variant_type == VariantType::DEL) {
            str += *structural + '\t';
            str += std::to_string(structural_size);
        }

        str += std::to_string(Observations()) + " (" + std::to_string(observations_fwd) + '/' + std::to_string(observations_rev) + ") ";
        str += std::to_string(position);
        str += " " + std::to_string(MeanQuality());
        str += '}';

        return str;
    }

    std::string ToMinimalString() const {
        std::string str;

        if (variant_type == VariantType::SNP) {
            str += std::string(1, variant) + "(";
        } else if (variant_type == VariantType::INS || variant_type == VariantType::DEL) {
            str += *structural + "(";
        }

        str += std::to_string(Observations());
        str += ", ";

        if (variant_type == VariantType::SNP) {
            str += (IsReference() ? "REF)" : "SNP)");
        } else if (variant_type == VariantType::INS) {
            str += "INS)";
        } else if (variant_type == VariantType::DEL) {
            str += "DEL)";
        }

        return str;
    }

    void AddObservation(Qual quality, bool from_forward) {
        observations_fwd += from_forward;
        observations_rev += !from_forward;
        if (Observations() == 65535) exit(9);
        quals.emplace_back(quality);
    }

    bool IsSNP() const {
        return variant_type == VariantType::SNP;
    }

    char Reference() const {
        return reference;
    }

    VariantPos Position() const {
        return position;
    }

    bool IsReference() const {
        return reference == variant;
    }

    bool HasFwdAndRev() const {
        return observations_fwd > 0 && observations_rev > 0;
    }

    void SetObservations(size_t obs) {
        observations_fwd = obs;
    }

    size_t GetStructuralSize() const {
        return structural_size;
    }

    std::string GetStructural() const {
        return *structural;
    }

    char GetVariant() const {
        return variant;
    }

    bool IsINDEL() const {
        return variant_type == VariantType::INS || variant_type == VariantType::DEL;
    }

    bool IsINS() const {
        return variant_type == VariantType::INS;
    }

    bool IsDEL() const {
        return variant_type == VariantType::DEL;
    }

    bool Match(VariantPos pos, Base base, Base ref) const {
        return pos == position &&
            base == variant && ref == reference;
    }

    bool Match(VariantType type, VariantPos pos, std::string& structural) const {
        return type == variant_type && pos == position &&
                *this->structural == structural;
    }

    bool Match(Variant const& other) const {
        return variant_type == other.variant_type &&
               position == other.position &&
                (variant_type == VariantType::SNP && variant == other.variant ||
                 variant_type == VariantType::INS && *structural == *other.structural);
    };

    void SetStructural(std::string&& structural_string) {
        structural_size = structural_string.length();
        structural = new std::string(structural_string);
    }
    void SetStructural(std::string& structural_string) {
        structural_size = structural_string.length();
        structural = new std::string(structural_string);
    }

    Variant() {}
};