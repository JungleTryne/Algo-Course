#include <algorithm>
#include <cassert>
#include <iostream>
#include <unordered_set>
#include <vector>

class GameSolution {
private:
    size_t n_;
    std::vector<size_t> positions_;
    std::vector<size_t> grandyFunction_;
public:
    explicit GameSolution(size_t n);
    void Solve();
    std::vector<size_t> GetTrivialPositions(size_t n) const;
    std::vector<size_t> GetPositions();
};

void GameSolution::Solve() {
    for(size_t i = 3; i < n_ + 1; ++i) {
        std::unordered_set<size_t> mex;
        for(size_t j = 3; j < i - 1; ++j) {
            size_t result = grandyFunction_[j - 1] ^ grandyFunction_[i - j];

            mex.insert(result);
            if(i == n_ && !result) {
                positions_.push_back(j);
            }
        }

        mex.insert(grandyFunction_[i - 1]);
        if(i == n_ && (!grandyFunction_[i - 1])) {
            positions_.insert(positions_.end(), {1, n_});
        }

        if(i == 3) {
            mex.insert(grandyFunction_[0]);
        } else {
            mex.insert(grandyFunction_[i - 2]);
            if(i == n_ && (!grandyFunction_[i - 2])) {
                positions_.insert(positions_.end(), {2, i-1});
            }
        }

        size_t mexValue = 0;
        while(mex.find(mexValue) != mex.end()) {
            ++mexValue;
        }
        grandyFunction_[i] = mexValue;
    }

    std::sort(positions_.begin(), positions_.end());
}

std::vector<size_t> GameSolution::GetPositions() {
    if(n_ < 4) {
        return GetTrivialPositions(n_);
    }
    return positions_;
}

std::vector<size_t> GameSolution::GetTrivialPositions(size_t n) const {
    if(n == 1) return {1};
    if(n == 2) return {1, 2};
    if(n == 3) return {2};
    assert(false);
}

GameSolution::GameSolution(size_t n) : n_(n) {
    grandyFunction_.resize(n_ + 1);
    grandyFunction_[0] = 0;
    grandyFunction_[1] = 1;
    grandyFunction_[2] = 1;
}

int main() {
    size_t n = 0;
    std::cin >> n;

    GameSolution solution(n);
    solution.Solve();
    auto result = solution.GetPositions();
    std::cout << (result.empty() ? "Mueller" : "Schtirlitz") << std::endl;
    for(size_t pos : result) {
        std::cout << pos << std::endl;
    }
    return 0;
}
