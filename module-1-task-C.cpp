#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <memory>

class TemplateFinder {
private:
    /* Trie node */
    struct Node {
        bool terminal_ = false;
        size_t word_size_ = 0;
        char parent_char_ = 0;

        std::shared_ptr<Node> parent_;
        std::shared_ptr<Node> suffix_;
        std::shared_ptr<Node> shrink_suffix_;

        std::vector<size_t> word_bias_; //Subtemplate bias. Subtemplates can be repeated -> several biases
        std::unordered_map<char, std::shared_ptr<Node>> transitions_;
        std::unordered_map<char, std::shared_ptr<Node>> delta_function_;
    };

    size_t subpattern_count_ = 0;
    std::string pattern_;

    std::shared_ptr<Node> root_;
    std::vector<std::shared_ptr<Node>> nodes_;

    void AddSubTemplate(const std::string& subtemplate, size_t word_bias);
    void ProcessShrunk(const std::shared_ptr<Node>& current_p, size_t char_pos, std::vector<size_t>& pattern_entries);

    auto GetSuffix(const std::shared_ptr<TemplateFinder::Node>& current_p) -> std::shared_ptr<Node>;
    auto GoDelta(const std::shared_ptr<Node>& current_p, char c) -> std::shared_ptr<Node> ;
    auto GetShrunkSuffix(const std::shared_ptr<Node>& current_p) -> std::shared_ptr<Node>;

    static auto UpdateEntries(const std::shared_ptr<Node>& current_p, size_t char_position,
                              std::vector<size_t>& pattern_entries) -> void;

    static auto split(const std::string& text, char splitter)
        -> std::pair<std::vector<std::string>, std::vector<size_t>>;
public:
    explicit TemplateFinder(const std::string& pattern);
    void FindEntries(const std::string& text, std::ostream& out);
};

/* Функция добавление подшаблона в бор */
void TemplateFinder::AddSubTemplate(const std::string &subtemplate, size_t word_bias) {
    auto p_current = root_;
    for (const auto& c : subtemplate) {
        if (p_current->transitions_.find(c) == p_current->transitions_.end()) {
            p_current->transitions_[c] = std::make_shared<Node>();
            p_current->transitions_[c]->parent_ = p_current;
            p_current->transitions_[c]->parent_char_ = c;
            nodes_.push_back(p_current->transitions_[c]);
        }
        p_current = p_current->transitions_[c];
    }
    p_current->terminal_ = true;
    p_current->word_bias_.push_back(word_bias);
    p_current->word_size_ = subtemplate.size();
    ++subpattern_count_;
}

TemplateFinder::TemplateFinder(const std::string& pattern) : pattern_(pattern) {
    root_ = std::make_shared<Node>();
    nodes_.push_back(root_);
    auto [split_text, bias] = split(pattern, '?');
    for (size_t i = 0; i < split_text.size(); ++i) {
        AddSubTemplate(split_text[i], bias[i]);
    }
}

/* Splitting the template to subtemplates */
auto TemplateFinder::split(const std::string &text, char splitter)
    -> std::pair<std::vector<std::string>, std::vector<size_t>>
{
    std::vector<std::string> split_text;
    std::vector<size_t> bias; //Position of subtemplates in the template
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

/* Getting suffix link of the node */
auto TemplateFinder::GetSuffix(const std::shared_ptr<TemplateFinder::Node>& current_p)
    -> std::shared_ptr<TemplateFinder::Node>
{
    if (!current_p->suffix_) {
        if (current_p == root_ || current_p->parent_ == root_) {
            current_p->suffix_ = root_;
        } else {
            current_p->suffix_ = GoDelta(GetSuffix(current_p->parent_), current_p->parent_char_);
        }
    }
    return current_p->suffix_;
}

/* Delta function of automata */
auto TemplateFinder::GoDelta(const std::shared_ptr<TemplateFinder::Node>& current_p, char c)
    -> std::shared_ptr<TemplateFinder::Node>
{
    if (current_p->delta_function_.find(c) == current_p->delta_function_.end()) {
        if (current_p->transitions_.find(c) != current_p->transitions_.end()) {
            current_p->delta_function_[c] = current_p->transitions_[c];
        } else if (current_p == root_) {
            current_p->delta_function_[c] = root_;
        } else {
            current_p->delta_function_[c] = GoDelta(GetSuffix(current_p), c);
        }
    }
    return current_p->delta_function_[c];
}

/* Getting shrunk suffix link of the node */
auto TemplateFinder::GetShrunkSuffix(const std::shared_ptr<TemplateFinder::Node>& current_p)
    -> std::shared_ptr<TemplateFinder::Node>
{
    if (!current_p->shrink_suffix_) {
        std::shared_ptr<Node> suffix_link = GetSuffix(current_p);
        if (suffix_link->terminal_) {
            current_p->shrink_suffix_ = suffix_link;
        } else if (suffix_link == root_) {
            current_p->shrink_suffix_ = root_;
        } else {
            current_p->shrink_suffix_ = GetShrunkSuffix(suffix_link);
        }
    }
    return current_p->shrink_suffix_;
}

/* Main algorithm function - finding pattern in the text  */
void TemplateFinder::FindEntries(const std::string &text, std::ostream& out) {
    std::shared_ptr<Node> current_p = root_;
    std::vector<size_t> pattern_entries(text.size());
    
    for (size_t char_pos = 0; char_pos < text.size(); ++char_pos) {
        current_p = GoDelta(current_p, text[char_pos]);
        ProcessShrunk(current_p, char_pos, pattern_entries);

        if (current_p->terminal_) {
            UpdateEntries(current_p, char_pos, pattern_entries);
        }
    }

    for (size_t char_pos = 0; char_pos < pattern_entries.size(); ++char_pos) {
        if (pattern_entries[char_pos] == subpattern_count_ && char_pos + pattern_.size() < text.size() + 1) {
            out << char_pos << ' ';
        }
    }
}

/* Shrunk suffix traversal */
auto TemplateFinder::ProcessShrunk(const std::shared_ptr<Node>& current_p, size_t char_pos,
                                   std::vector<size_t> &pattern_entries) -> void
{
    for (auto shrunk_p = GetShrunkSuffix(current_p); shrunk_p != root_; shrunk_p = GetShrunkSuffix(shrunk_p)) {
        UpdateEntries(shrunk_p, char_pos, pattern_entries);
    }
}

auto TemplateFinder::UpdateEntries(const std::shared_ptr<Node> &current_p, size_t char_pos,
                                   std::vector<size_t> &pattern_entries) -> void
{
    auto update_entries = [current_p, char_pos, &pattern_entries](size_t bias) {
        auto pattern_pos = static_cast<int64_t>(char_pos - bias - current_p->word_size_ + 1);
        if (pattern_pos >= 0 && pattern_pos < static_cast<int64_t>(pattern_entries.size())) {
            ++pattern_entries[static_cast<size_t>(pattern_pos)];
        }
    };
    std::for_each(current_p->word_bias_.begin(), current_p->word_bias_.end(), update_entries);
}

int main() {
    std::string text_template;
    std::string text;
    std::cin >> text_template >> text;

    TemplateFinder finder(text_template);
    finder.FindEntries(text, std::cout);

    std::cout << std::endl;
    return 0;
}
