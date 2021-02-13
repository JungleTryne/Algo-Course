#include <algorithm>
#include <iostream>
#include <vector>

class SuffixArray {
private:
    size_t alpha_size_ = 128;
    char terminal_char_ = '$';

    std::string text_;
    std::vector<size_t> suffix_;
    std::vector<int64_t> lcp_;

    std::vector<size_t> classes_;
    std::vector<size_t> new_classes_;
    size_t current_length_ = 0;

    inline size_t GetSuffixLength(size_t suf_ind);
    void BuildInitialBuckets();
    void BucketSortLoop();
    void MergeClasses();
    void LCP();
    void BuildSuffixArray();
public:
    explicit SuffixArray(std::string text);
    int64_t GetDistinctStringsNumber();
};

SuffixArray::SuffixArray(std::string text) : text_(std::move(text)) {
    BuildSuffixArray();
    LCP();
}

/*
 * Dividing the string into classes
 * in order to perform bucket sort
 */
void SuffixArray::BuildInitialBuckets() {
    text_ += terminal_char_;

    std::vector<size_t> counter((text_.size() > alpha_size_) ? text_.size() : alpha_size_, 0);
    for(char ch : text_) {
        ++counter[static_cast<size_t>(ch)];
    }

    for(size_t i = 1; i < counter.size(); ++i) {
        counter[i] += counter[i-1];
    }

    suffix_.resize(text_.size());
    for(size_t i = 0; i < text_.size(); ++i) {
        suffix_[--counter[static_cast<size_t>(text_[i])]] = i;
    }

    classes_.resize(text_.size(), 0);
    size_t class_number = 0;
    for(unsigned long suff : suffix_) {
        if(text_[suff] != terminal_char_) {
            terminal_char_ = text_[suff];
            ++class_number;
        }
        classes_[suff] = class_number;
    }
}

void SuffixArray::BuildSuffixArray() {
    BuildInitialBuckets();
    new_classes_.resize(text_.size(), 0);
    current_length_ = 1;

    while(current_length_ <= text_.size()) {
        BucketSortLoop();
    }
}

void SuffixArray::BucketSortLoop() {
    std::vector<size_t> suffix_second(text_.size());
    for(size_t i = 0; i < suffix_second.size(); ++i) {
        auto biased_second_suffix = static_cast<int64_t>(suffix_[i] - current_length_);
        suffix_second[i] = biased_second_suffix < 0 ?
                           text_.size() + static_cast<size_t>(biased_second_suffix) :
                           static_cast<size_t>(biased_second_suffix);
    }

    std::vector<size_t> counter(classes_.size());

    std::for_each(suffix_second.begin(), suffix_second.end(), [this, &counter](size_t pos) {
        ++counter[classes_[pos]];
    });

    for(size_t i = 1; i < counter.size(); ++i) {
        counter[i] += counter[i-1];
    }

    std::vector<size_t> new_suffix(suffix_.size());
    for(auto i = static_cast<int64_t>(suffix_second.size() - 1); i > -1; --i) {
        auto index = static_cast<size_t>(i);
        new_suffix[--counter[classes_[suffix_second[index]]]] = suffix_second[index];
    }
    suffix_ = std::move(new_suffix);
    MergeClasses();
}

void SuffixArray::MergeClasses() {
    new_classes_.clear();
    new_classes_.resize(text_.size());
    size_t new_classes_number = 0;
    new_classes_[0] = 0;
    for(size_t i = 1; i < text_.size(); ++i) {
        size_t mid_1 = (suffix_[i] + current_length_) % text_.size();
        size_t mid_2 = (suffix_[i - 1] + current_length_) % text_.size();
        if(classes_[suffix_[i]] != classes_[suffix_[i - 1]] || classes_[mid_1] != classes_[mid_2]) {
            ++new_classes_number;
        }
        new_classes_[suffix_[i]] = new_classes_number;
    }
    classes_ = new_classes_;
    current_length_ *= 2;
}

void SuffixArray::LCP() {
    std::vector<size_t> pos(suffix_.size());
    lcp_.resize(suffix_.size());

    for(size_t i = 0; i < suffix_.size(); ++i) {
        pos[suffix_[i]] = i;
    }

    int64_t current_bias = 0;
    for(size_t i = 0; i < text_.size() - 1; ++i) {
        if(current_bias > 0) --current_bias;
        if(pos[i] == text_.size() - 1) {
            lcp_[text_.size() - 1] = -1;
            current_bias = 0;
        } else {
            size_t j = suffix_[pos[i] + 1];
            int64_t index_one = static_cast<int64_t>(i) + current_bias;
            int64_t index_two = static_cast<int64_t>(j) + current_bias;

            while(std::max(index_one, index_two) < static_cast<int64_t>(text_.size()) &&
                  text_[static_cast<size_t>(index_one)] == text_[static_cast<size_t>(index_two)]) {
                ++current_bias;
                index_one = static_cast<int64_t>(i) + current_bias;
                index_two = static_cast<int64_t>(j) + current_bias;
            }
            lcp_[pos[i]] = current_bias;
        }
    }
}

int64_t SuffixArray::GetDistinctStringsNumber() {
    auto count = static_cast<int64_t>(GetSuffixLength(0));
    for(size_t i = 1; i < suffix_.size(); ++i) {
        count += (static_cast<int64_t>(GetSuffixLength(i)) - lcp_[i]);
    }
    return count - static_cast<int64_t>(text_.size()) - 1;
}

size_t SuffixArray::GetSuffixLength(size_t suf_ind) {
    return text_.size() - suffix_[suf_ind];
}

int main() {
    std::string text;
    std::cin >> text;
    SuffixArray sf(text);
    std::cout << sf.GetDistinctStringsNumber() << std::endl;
    return 0;
}
