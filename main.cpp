#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

class SubstringFinder {
private:
    size_t text_char_counter = -1;
    bool buffer_ready = false;

    std::string pattern;
    std::string buffer;
    std::vector<size_t> prefix_func;
    std::vector<size_t> pattern_positions;

    void FindPrefixFunc();
    static inline size_t GetPrefixIndex(size_t text_index, size_t patter_size);
public:
    explicit SubstringFinder(const std::string& pattern);
    void AddTextCharacter(char new_char);
    std::vector<size_t> get_pattern_positions();
};

inline size_t SubstringFinder::GetPrefixIndex(size_t text_index, size_t patter_size) {
    return (text_index % patter_size) + patter_size + 1;
}

void SubstringFinder::FindPrefixFunc() {
    for(size_t i = 1; i < prefix_func.size(); ++i) {
        size_t recursive_index = prefix_func[i - 1];

        if(buffer[recursive_index] == buffer[i]) {
            prefix_func[i] = prefix_func[i - 1] + 1;
            continue;
        }
        while(recursive_index > 0 && buffer[recursive_index] != buffer[i]) {
            recursive_index = prefix_func[recursive_index - 1];
        }
        if(buffer[recursive_index] == buffer[i]) {
            ++recursive_index;
        }
        prefix_func[i] = recursive_index;
    }
}

SubstringFinder::SubstringFinder(const std::string& pattern) : pattern(pattern) {
    buffer = (pattern + "$");
}

void SubstringFinder::AddTextCharacter(char new_char) {
    ++text_char_counter;
    while(text_char_counter < pattern.size() - 1) {
        buffer += new_char;
        return;
    }

    if(!buffer_ready) {
        buffer_ready = true;
        prefix_func.resize(buffer.size());
        FindPrefixFunc();
        size_t prefix_index = GetPrefixIndex(pattern.size() - 1, pattern.size());
        if(prefix_func[prefix_index] == pattern.size()) {
            pattern_positions.push_back(0);
        }
    }

    size_t prefix_index = GetPrefixIndex(text_char_counter - 1, pattern.size());
    prefix_index = prefix_func[prefix_index];

    if(buffer[prefix_index] == new_char) {
        prefix_func[GetPrefixIndex(text_char_counter, pattern.size())] = prefix_index + 1;
        if(prefix_func[GetPrefixIndex(text_char_counter, pattern.size())] == pattern.size()) {
            pattern_positions.push_back(text_char_counter - pattern.size() + 1);
        }
        return;
    }

    while(prefix_index > 0 && buffer[prefix_index] != new_char) {
        prefix_index = prefix_func[prefix_index - 1];
    }
    if(buffer[prefix_index] == new_char) {
        ++prefix_index;
    }
    prefix_func[GetPrefixIndex(text_char_counter, pattern.size())] = prefix_index;
    if(prefix_func[GetPrefixIndex(text_char_counter, pattern.size())] == pattern.size()) {
        pattern_positions.push_back(text_char_counter - pattern.size() + 1);
    }

}

std::vector<size_t> SubstringFinder::get_pattern_positions() {
    return pattern_positions;
}

int main() {
    std::string pattern;
    std::vector<size_t> pattern_positions;

    std::cin >> pattern;
    SubstringFinder finder(pattern);

    getchar();
    char next_char = static_cast<char>(getchar());
    do {
        finder.AddTextCharacter(next_char);
        next_char = static_cast<char>(getchar());
    } while(next_char > 15);

    for(size_t pos : finder.get_pattern_positions()) {
        std::cout << pos << ' ';
    }
    std::cout << std::endl;

    return 0;
}
