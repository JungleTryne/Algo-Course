#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <memory>

class TemplateFinder {
private:
    /* Вершина бора */
    struct Node {
        bool terminal = false;
        size_t word_size = 0;
        char parent_char = 0;

        std::shared_ptr<Node> parent;
        std::shared_ptr<Node> suffix;
        std::shared_ptr<Node> shrink_suffix;

        std::vector<size_t> word_bias; //Смещение подшаблона. Подшаблоны могут повторяться -> массив
        std::unordered_map<char, std::shared_ptr<Node>> transitions;
        std::unordered_map<char, std::shared_ptr<Node>> delta_function;
    };

    size_t subpattern_count = 0;
    std::string pattern;

    std::shared_ptr<Node> root;
    std::vector<std::shared_ptr<Node>> nodes;

    void AddSubTemplate(const std::string& subtemplate, size_t word_bias);
    void ProcessShrunk(const std::shared_ptr<Node>& current_p, size_t char_pos, std::vector<size_t>& pattern_entries);

    auto GetSuffix(const std::shared_ptr<TemplateFinder::Node>& current_p) -> std::shared_ptr<Node>;
    auto GoDelta(const std::shared_ptr<Node>& current_p, char c) -> std::shared_ptr<Node> ;
    auto GetShrunkSuffix(const std::shared_ptr<Node>& current_p) -> std::shared_ptr<Node>;

    static auto split(const std::string& text, char splitter) -> std::pair<std::vector<std::string>, std::vector<size_t>>;
public:
    explicit TemplateFinder(const std::string& pattern);
    void FindEntries(const std::string& text, std::ostream& out);
};

/* Функция добавление подшаблона в бор */
void TemplateFinder::AddSubTemplate(const std::string &subtemplate, size_t word_bias) {
    auto p_current = root;
    for (const auto& c : subtemplate) {
        if (p_current->transitions.find(c) == p_current->transitions.end()) {
            p_current->transitions[c] = std::make_shared<Node>();
            p_current->transitions[c]->parent = p_current;
            p_current->transitions[c]->parent_char = c;
            nodes.push_back(p_current->transitions[c]);
        }
        p_current = p_current->transitions[c];
    }
    p_current->terminal = true;
    p_current->word_bias.push_back(word_bias);
    p_current->word_size = subtemplate.size();
    ++subpattern_count;
}

TemplateFinder::TemplateFinder(const std::string& pattern) : pattern(pattern) {
    root = std::make_shared<Node>();
    nodes.push_back(root);
    auto [split_text, bias] = split(pattern, '?');
    for (size_t i = 0; i < split_text.size(); ++i) {
        AddSubTemplate(split_text[i], bias[i]);
    }
}

/* Функция разбиентия шаблона на подшаблоны */
auto TemplateFinder::split(const std::string &text, char splitter) -> std::pair<std::vector<std::string>, std::vector<size_t>>  {
    std::vector<std::string> split_text;
    std::vector<size_t> bias; //Позиция подшаблонов в шаблоне
    std::string buffer;

    size_t counter = 0;
    for (const auto& c : text) {
        if (c == splitter && !buffer.empty()) {
            bias.push_back(counter - buffer.size());
            split_text.push_back(buffer);
            buffer = "";
        } else if (c != splitter) {
            buffer += c;
        }
        ++counter;
    }
    if (!buffer.empty()) {
        bias.push_back(counter - buffer.size());
        split_text.push_back(buffer);
    }
    return std::make_pair(split_text, bias);
}

/* Получение суффиксной ссылки у вершины */
auto TemplateFinder::GetSuffix(const std::shared_ptr<TemplateFinder::Node>& current_p) -> std::shared_ptr<TemplateFinder::Node> {
    if (!current_p->suffix) {
        if (current_p == root || current_p->parent == root) {
            current_p->suffix = root;
        } else {
            current_p->suffix = GoDelta(GetSuffix(current_p->parent), current_p->parent_char);
        }
    }
    return current_p->suffix;
}

/* Функция перехода в автомате */
auto TemplateFinder::GoDelta(const std::shared_ptr<TemplateFinder::Node>& current_p, char c) -> std::shared_ptr<TemplateFinder::Node> {
    if (current_p->delta_function.find(c) == current_p->delta_function.end()) {
        if (current_p->transitions.find(c) != current_p->transitions.end()) {
            current_p->delta_function[c] = current_p->transitions[c];
        } else if (current_p == root) {
            current_p->delta_function[c] = root;
        } else {
            current_p->delta_function[c] = GoDelta(GetSuffix(current_p), c);
        }
    }
    return current_p->delta_function[c];
}

/* Функция получения сжатой суффиксной ссылки */
auto TemplateFinder::GetShrunkSuffix(const std::shared_ptr<TemplateFinder::Node>& current_p) -> std::shared_ptr<TemplateFinder::Node> {
    if (!current_p->shrink_suffix) {
        std::shared_ptr<Node> suffix_link = GetSuffix(current_p);
        if (suffix_link->terminal) {
            current_p->shrink_suffix = suffix_link;
        } else if (suffix_link == root) {
            current_p->shrink_suffix = root;
        } else {
            current_p->shrink_suffix = GetShrunkSuffix(suffix_link);
        }
    }
    return current_p->shrink_suffix;
}

/* Основная функция алгоритма - найти шаблон в заданном тексте */
void TemplateFinder::FindEntries(const std::string &text, std::ostream& out) {
    std::shared_ptr<Node> current_p = root;
    std::vector<size_t> pattern_entries(text.size());
    
    for (size_t char_pos = 0; char_pos < text.size(); ++char_pos) {
        current_p = GoDelta(current_p, text[char_pos]);
        ProcessShrunk(current_p, char_pos, pattern_entries);

        if (current_p->terminal) {
            auto update_entries = [current_p, char_pos, &pattern_entries](size_t bias) {
                size_t pattern_pos = char_pos - bias - current_p->word_size + 1;
                if (pattern_pos >= 0 && pattern_pos < pattern_entries.size()) {
                    ++pattern_entries[pattern_pos];
                }
            };
            std::for_each(current_p->word_bias.begin(), current_p->word_bias.end(), update_entries);
        }
    }

    for (size_t char_pos = 0; char_pos < pattern_entries.size(); ++char_pos) {
        if (pattern_entries[char_pos] == subpattern_count && char_pos + pattern.size() < text.size() + 1) {
            out << char_pos << ' ';
        }
    }
}

/* Фукнция прохождения по сжатым суффиксным ссылкам */
void TemplateFinder::ProcessShrunk(const std::shared_ptr<Node>& current_p, size_t char_pos, std::vector<size_t> &pattern_entries) {
    for (auto shrunk_p = GetShrunkSuffix(current_p); shrunk_p != root; shrunk_p = GetShrunkSuffix(shrunk_p)) {
        for (auto bias : shrunk_p->word_bias) {
            size_t pattern_pos = char_pos - bias - shrunk_p->word_size + 1;
            if (pattern_pos >= 0 && pattern_pos < pattern_entries.size()) {
                ++pattern_entries[pattern_pos];
            }
        }
    }
}

int main() {
    std::string templ;
    std::string text;
    std::cin >> templ >> text;
    TemplateFinder finder(templ);
    finder.FindEntries(text, std::cout);
    std::cout << '\n';
    return 0;
}
