#include <iostream>
#include <vector>
#include <map>

class SuffixArray {
private:
    std::string text;
    std::vector<size_t> suffix;
    std::vector<int64_t> lcp;

    inline size_t GetSuffixLength(size_t suf_ind);
    void LCP();
    void BuildSuffixArray();
public:
    explicit SuffixArray(std::string text);
    int64_t GetDistinctStringsNumber();
};

SuffixArray::SuffixArray(std::string text) : text(std::move(text)) {
    BuildSuffixArray();
    LCP();
}

void SuffixArray::BuildSuffixArray() {
    constexpr size_t ALPHA_SIZE = 128;
    text += '$';
    std::vector<size_t> counter((text.size() > ALPHA_SIZE) ? text.size() : ALPHA_SIZE, 0);
    std::for_each(text.begin(), text.end(), [&counter](char ch) {
        ++counter[static_cast<size_t>(ch)];
    });

    for(size_t i = 1; i < counter.size(); ++i) {
        counter[i] += counter[i-1];
    }

    suffix.resize(text.size());
    for(size_t i = 0; i < text.size(); ++i) {
        suffix[--counter[static_cast<size_t>(text[i])]] = i;
    }

    std::vector<size_t> classes(text.size(), 0);
    size_t class_number = 0;
    char last_char = '$';
    for(unsigned long suff : suffix) {
        if(text[suff] != last_char) {
            last_char = text[suff];
            ++class_number;
        }
        classes[suff] = class_number;
    }

    std::vector<size_t> new_classes(text.size(), 0);
    size_t current_length = 1;

    while(current_length <= text.size()) {
        std::vector<size_t> suffix_second(text.size());
        for(size_t i = 0; i < suffix_second.size(); ++i) {
            auto biased_second_suffix = static_cast<int64_t>(suffix[i] - current_length);
            suffix_second[i] = biased_second_suffix < 0 ?
                    text.size() + static_cast<size_t>(biased_second_suffix):
                    static_cast<size_t>(biased_second_suffix);
        }

        counter.clear();
        counter.resize(classes.size());

        std::for_each(suffix_second.begin(), suffix_second.end(), [&classes, &counter](size_t pos) {
            ++counter[classes[pos]];
        });

        for(size_t i = 1; i < counter.size(); ++i) {
            counter[i] += counter[i-1];
        }

        std::vector<size_t> new_suffix(suffix.size());
        for(auto i = static_cast<int64_t>(suffix_second.size() - 1); i > -1; --i) {
            auto index = static_cast<size_t>(i);
            new_suffix[--counter[classes[suffix_second[index]]]] = suffix_second[index];
        }
        suffix = std::move(new_suffix);

        new_classes.clear();
        new_classes.resize(text.size());
        size_t new_classes_number = 0;
        new_classes[0] = 0;
        for(size_t i = 1; i < text.size(); ++i) {
            size_t mid_1 = (suffix[i] + current_length) % text.size();
            size_t mid_2 = (suffix[i-1] + current_length) % text.size();
            if(classes[suffix[i]] != classes[suffix[i-1]] || classes[mid_1] != classes[mid_2]) {
                ++new_classes_number;
            }
            new_classes[suffix[i]] = new_classes_number;
        }
        classes = new_classes;
        current_length *= 2;
    }
}

void SuffixArray::LCP() {
    std::vector<size_t> pos(suffix.size());
    lcp.resize(suffix.size());

    for(size_t i = 0; i < suffix.size(); ++i) {
        pos[suffix[i]] = i;
    }

    int64_t current_bias = 0;
    for(size_t i = 0; i < text.size() - 1; ++i) {
        if(current_bias > 0) --current_bias;
        if(pos[i] == text.size() - 1) {
            lcp[text.size() - 1] = -1;
            current_bias = 0;
        } else {
            size_t j = suffix[pos[i] + 1];
            int64_t index_one = static_cast<int64_t>(i) + current_bias;
            int64_t index_two = static_cast<int64_t>(j) + current_bias;

            while(std::max(index_one, index_two) < static_cast<int64_t>(text.size()) &&
                  text[static_cast<size_t>(index_one)] == text[static_cast<size_t>(index_two)]) {
                ++current_bias;
                index_one = static_cast<int64_t>(i) + current_bias;
                index_two = static_cast<int64_t>(j) + current_bias;
            }
            lcp[pos[i]] = current_bias;
        }
    }
}

int64_t SuffixArray::GetDistinctStringsNumber() {
    auto count = static_cast<int64_t>(GetSuffixLength(0));
    for(size_t i = 1; i < suffix.size(); ++i) {
        count += (static_cast<int64_t>(GetSuffixLength(i)) - lcp[i]);
    }
    return count - static_cast<int64_t>(text.size()) - 1;
}

size_t SuffixArray::GetSuffixLength(size_t suf_ind) {
    return text.size() - suffix[suf_ind];
}

int main() {
    std::string text;
    std::cin >> text;
    SuffixArray sf(text);
    std::cout << sf.GetDistinctStringsNumber() << std::endl;
    return 0;
}
