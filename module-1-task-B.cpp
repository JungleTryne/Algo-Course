#include <iostream>
#include <vector>

class PalindromeHandler {
private:
    std::string text;
    std::vector<int64_t> odd_palindrome;
    std::vector<int64_t> even_palindrome;
public:
    explicit PalindromeHandler(std::string text);
    void GetAmount(bool even);
    int64_t GetAmountOfPalindromes();
};

PalindromeHandler::PalindromeHandler(std::string text) : text(std::move(text)) {}

void PalindromeHandler::GetAmount(bool even) {
    int64_t left = 0;
    int64_t right = -1;

    std::vector<int64_t> palindrome_count(text.size());
    auto bias = static_cast<int64_t>(even);
    
    for(int64_t i = 0; i < static_cast<int64_t>(text.size()); ++i) {
        auto answer = static_cast<int64_t>(!even); //Как минимум, палиндром имеет длину 1, если нечет, и 0, если четный палиндром
        if(i <= right) {
            int64_t potential_update = palindrome_count[left + (right - i) + bias];
            answer = right - i + 1 < potential_update ? right - i + 1 : potential_update;
        }
        while(i - answer - bias > -1 && i + answer < text.size()) { //Пытаемся наивно улучшить результат
            if(text[i - answer - bias] == text[i + answer]) {
                ++answer;
            } else {
                break;
            }
        }
        palindrome_count[i] = answer;
        if(right < i + answer - 1) { //Обновляем правую и левую границы
            right = i + answer - 1;
            left = i - answer + (1 - bias) ;
        }
    }

    (even ? even_palindrome : odd_palindrome) = palindrome_count;
}

int64_t PalindromeHandler::GetAmountOfPalindromes() {
    GetAmount(false);
    GetAmount(true);

    int64_t counter = 0;
    for(int64_t palindrome_size : odd_palindrome) {
        counter += (palindrome_size - 1);
    }
    for(int64_t palindrome_size : even_palindrome) {
        counter += palindrome_size;
    }
    return counter;
}

int main() {
    std::string text;
    std::cin >> text;
    PalindromeHandler pHandler(text);
    std::cout << pHandler.GetAmountOfPalindromes();
    return 0;
}
