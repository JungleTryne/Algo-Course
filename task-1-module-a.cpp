#include <cstdio>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

template <typename OutputIterator>
class SubstringFinder {
private:
    long long text_char_counter_ = -1;
    bool buffer_ready_ = false;

    std::string pattern_;
    std::string buffer_;
    std::vector<size_t> prefix_func_;
    OutputIterator& out_stream_;

    void ProceedPrefixFuncLoop(size_t last_char_index, size_t new_char_index, char new_char);
    void FindPrefixFunc();
    void UpdatePrefixFunc(char new_char);
    static inline size_t GetPrefixIndex(size_t text_index, size_t pattern_size);
public:
    explicit SubstringFinder(const std::string& pattern, OutputIterator& out_stream);
    void AddTextCharacter(char new_char);
    void SavePatternPosition(size_t pos) ;
};

template <typename OutputIterator>
SubstringFinder<OutputIterator>::SubstringFinder(const std::string& pattern, OutputIterator& out_stream):
                                                                                        pattern_(pattern),
                                                                                        buffer_(pattern + "$"),
                                                                                        out_stream_(out_stream)
{}

template <typename OutputIterator>
inline size_t SubstringFinder<OutputIterator>::GetPrefixIndex(size_t text_index, size_t pattern_size) {
    return (text_index % pattern_size) + pattern_size + 1;
}

template <typename OutputIterator>
void SubstringFinder<OutputIterator>::AddTextCharacter(char new_char) {
    ++text_char_counter_;
    //Filling up the buffer
    while (static_cast<size_t>(text_char_counter_) < pattern_.size() - 1) {
        buffer_ += new_char;
        return;
    }

    //As long as the buffer is full, we can start calculation of prefix function
    if (!buffer_ready_) {
        buffer_ready_ = true;
        prefix_func_.resize(buffer_.size());
        FindPrefixFunc();
        if (prefix_func_.back() == pattern_.size()) {
            *out_stream_ = 0;
        }
    }
    UpdatePrefixFunc(new_char);
}

template <typename OutputIterator>
void SubstringFinder<OutputIterator>::ProceedPrefixFuncLoop(size_t last_char_index,size_t new_char_index,
                                                            char new_char) {
    size_t prefix_index = prefix_func_[last_char_index];

    if (buffer_[prefix_index] == new_char) {
        prefix_func_[new_char_index] = prefix_index + 1;
        if (prefix_func_[new_char_index] == pattern_.size()) {
            size_t pattern_start_pos = static_cast<size_t>(text_char_counter_) - pattern_.size() + 1;
            SavePatternPosition(pattern_start_pos);
        }
        return;
    }

    while (prefix_index > 0 && buffer_[prefix_index] != new_char) {
        prefix_index = prefix_func_[prefix_index - 1];
    }

    if (buffer_[prefix_index] == new_char) {
        ++prefix_index;
    }

    prefix_func_[new_char_index] = prefix_index;
    if (prefix_func_[new_char_index] == pattern_.size()) {
        size_t pattern_start_pos = static_cast<size_t>(text_char_counter_) - pattern_.size() + 1;
        SavePatternPosition(pattern_start_pos);
    }
}

template <typename OutputIterator>
void SubstringFinder<OutputIterator>::UpdatePrefixFunc(char new_char) {
    size_t last_char_index = GetPrefixIndex(static_cast<size_t>(text_char_counter_) - 1,pattern_.size());
    size_t new_prefix_index = GetPrefixIndex(static_cast<size_t>(text_char_counter_), pattern_.size());

    ProceedPrefixFuncLoop(last_char_index, new_prefix_index, new_char);
}

template <typename OutputIterator>
void SubstringFinder<OutputIterator>::FindPrefixFunc() {
    for (size_t i = 1; i < prefix_func_.size(); ++i) {
        ProceedPrefixFuncLoop(i - 1, i, buffer_[i]);
    }
}

template <typename OutputIterator>
void SubstringFinder<OutputIterator>::SavePatternPosition(size_t pos) {
    *out_stream_ = pos;
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::string pattern;
    std::vector<size_t> pattern_positions;

    std::cin >> pattern;
    auto out_iter = std::ostream_iterator<size_t>(std::cout, " ");
    SubstringFinder finder(pattern, out_iter);

    char input;
    while (std::cin >> input) {
        finder.AddTextCharacter(input);
    }
    std::cout << std::endl;
    return 0;
}
