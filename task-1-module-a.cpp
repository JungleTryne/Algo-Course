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
    std::ostream& out_stream;

    void ProceedPrefixFuncLoop(size_t last_char_index, size_t new_char_index, char new_char);
    void FindPrefixFunc();
    void UpdatePrefixFunc(char new_char);
    static inline size_t GetPrefixIndex(size_t text_index, size_t pattern_size);
public:
    explicit SubstringFinder(const std::string& pattern, std::ostream& out_stream);
    void AddTextCharacter(char new_char);
    void SavePatternPosition(size_t pos) ;
};

SubstringFinder::SubstringFinder(const std::string& pattern, std::ostream& out_stream): pattern(pattern),
                                                                                        out_stream(out_stream),
                                                                                        buffer(pattern + "$")
{}

/* Функия получения циклического индекса во время обработки символа текста под индексом text_index */
inline size_t SubstringFinder::GetPrefixIndex(size_t text_index, size_t pattern_size) {
    return (text_index % pattern_size) + pattern_size + 1;
}

void SubstringFinder::AddTextCharacter(char new_char) {
    ++text_char_counter;
    //Заполняем буффер
    while(text_char_counter < pattern.size() - 1) {
        buffer += new_char;
        return;
    }

    //Буффер заполнился и можно начинать вычислять префикс функцию
    if(!buffer_ready) {
        buffer_ready = true;
        prefix_func.resize(buffer.size());
        FindPrefixFunc();
        if (prefix_func.back() == pattern.size()) {
            std::cout << 0 << ' ';
        }
    }

    UpdatePrefixFunc(new_char);

}

void SubstringFinder::ProceedPrefixFuncLoop(size_t last_char_index, size_t new_char_index, char new_char) {
    size_t prefix_index = prefix_func[last_char_index];

    if(buffer[prefix_index] == new_char) {
        prefix_func[new_char_index] = prefix_index + 1;
        if(prefix_func[new_char_index] == pattern.size()) {
            size_t pattern_start_pos = text_char_counter - pattern.size() + 1;
            SavePatternPosition(pattern_start_pos);
        }
        return;
    }

    while(prefix_index > 0 && buffer[prefix_index] != new_char) {
        prefix_index = prefix_func[prefix_index - 1];
    }

    if(buffer[prefix_index] == new_char) {
        ++prefix_index;
    }

    prefix_func[new_char_index] = prefix_index;
    if(prefix_func[new_char_index] == pattern.size()) {
        size_t pattern_start_pos = text_char_counter - pattern.size() + 1;
        SavePatternPosition(pattern_start_pos);
    }
}

void SubstringFinder::UpdatePrefixFunc(char new_char) {
    size_t last_char_index = GetPrefixIndex(text_char_counter - 1, pattern.size());
    size_t new_prefix_index = GetPrefixIndex(text_char_counter, pattern.size());

    ProceedPrefixFuncLoop(last_char_index, new_prefix_index, new_char);
}

void SubstringFinder::FindPrefixFunc() {
    for(size_t i = 1; i < prefix_func.size(); ++i) {
        ProceedPrefixFuncLoop(i - 1, i, buffer[i]);
    }
}

void SubstringFinder::SavePatternPosition(size_t pos) {
    out_stream << pos << ' ';
}

int main() {
    std::string pattern;
    std::vector<size_t> pattern_positions;

    std::cin >> pattern;
    SubstringFinder finder(pattern, std::cout);

    getchar();
    char next_char = static_cast<char>(getchar());
    do {
        finder.AddTextCharacter(next_char);
        next_char = static_cast<char>(getchar());
    } while(next_char > 15);

    std::cout << std::endl;

    return 0;
}
