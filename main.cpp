#include <algorithm>
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
    explicit SuffixArray(std::string&& first, std::string&& second);
    std::string GetKthCommonSubstring(size_t k);
};

SuffixArray::SuffixArray(std::string&& first, std::string&& second) : text(first + "#" + second) {
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

size_t SuffixArray::GetSuffixLength(size_t suf_ind) {
    return text.size() - suffix[suf_ind];
}

std::string SuffixArray::GetKthCommonSubstring(size_t k) {
    size_t substring_counter = 0;
    size_t first_word_size = suffix[0];
    bool on_first_word = false;
    for(size_t i = 1; i < text.size(); ++i) {
        size_t j = i + 1;
        int64_t min_lcp = lcp[i-1];
        while(on_first_word ? (suffix[j] < first_word_size) : (suffix[j] >= first_word_size)) {
            min_lcp = min_lcp > lcp[j-1] ? lcp[j-1] : min_lcp;
            ++j;
            if(j >= text.size()) return "";
        }
        min_lcp = min_lcp < lcp[j-1] ? min_lcp : lcp[j-1];
        //size_t word_size = !on_first_word ? first_word_size : suffix.size() - 1;
        size_t update = static_cast<size_t>(lcp[j - 1]) - static_cast<size_t>(min_lcp);
        if(update + substring_counter >= k) {
            size_t first = suffix[j];
            size_t length = static_cast<size_t>(min_lcp) + k - substring_counter;
            return text.substr(first, length);
        } else {
            on_first_word = !on_first_word;
            substring_counter += update;
            i = j - 1;
        }
    }
    return "";
}

int main() {
    std::string str1;
    std::string str2;
    std::cin >> str1 >> str2;
    SuffixArray sf(std::move(str1), std::move(str2));
    size_t k = 0;
    std::cin >> k;
    std::string result = sf.GetKthCommonSubstring(k);
    if (result.empty()) {
        std::cout << -1;
    } else {
        std::cout << result;
    }
    std::cout << '\n';


    return 0;
}
