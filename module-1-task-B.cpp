#include <iostream>
#include <vector>

class PalindromeHandler {
private:
    std::string text_;
    std::vector<int64_t> odd_palindrome_;
    std::vector<int64_t> even_palindrome_;
public:
    explicit PalindromeHandler(std::string text);
    void GetAmount(bool even);
    int64_t GetAmountOfPalindromes();
};

PalindromeHandler::PalindromeHandler(std::string text) : text_(std::move(text)) {}

void PalindromeHandler::GetAmount(bool even) {
    int64_t left = 0;
    int64_t right = -1;

    std::vector<int64_t> palindrome_count(text_.size());
    auto bias = static_cast<int64_t>(even);
    
    for(int64_t i = 0; i < static_cast<int64_t>(text_.size()); ++i) {
        auto answer = static_cast<int64_t>(!even); //Min palindrome length is 1 (for odd-sized) or 0 (for even-sized)
        if(i <= right) {
            int64_t potential_update = palindrome_count[static_cast<size_t>(left + (right - i) + bias)];
            answer = right - i + 1 < potential_update ? right - i + 1 : potential_update;
        }
        //Trying to improve the answer in naive way
        while(i - answer - bias > -1 && static_cast<size_t>(i + answer) < text_.size()) {
            if(text_[static_cast<size_t>(i - answer - bias)] == text_[static_cast<size_t>(i + answer)]) {
                ++answer;
            } else {
                break;
            }
        }
        palindrome_count[static_cast<size_t>(i)] = answer;
        //Updating left and right bounds
        if(right < i + answer - 1) {
            left = i - answer + (1 - bias) ;
            right = i + answer - 1;
        }
    }

    (even ? even_palindrome_ : odd_palindrome_) = palindrome_count;
}

int64_t PalindromeHandler::GetAmountOfPalindromes() {
    GetAmount(false);
    GetAmount(true);

    int64_t counter = 0;
    for(int64_t palindrome_size : odd_palindrome_) {
        counter += (palindrome_size - 1);
    }
    for(int64_t palindrome_size : even_palindrome_) {
        counter += palindrome_size;
    }
    return counter;
}

int main() {
    std::string text;
    std::cin >> text;
    PalindromeHandler p_handler(text);
    std::cout << p_handler.GetAmountOfPalindromes();
    return 0;
}
