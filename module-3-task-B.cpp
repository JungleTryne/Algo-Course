#include <algorithm>
#include <cmath>
#include <set>
#include <vector>
#include <iostream>
#include <iomanip>

//Point struct
template<typename Type>
struct GeneralPoint {
    Type x;
    Type y;

    explicit GeneralPoint(Type x=0, Type y=0) : x(x), y(y) {};
    GeneralPoint(const GeneralPoint& other) = default;
    ~GeneralPoint() = default;
};

using PointType = int64_t;
using Point = GeneralPoint<PointType>;

int64_t CrossProduct(const Point& first, const Point& second, const Point& third) {
    return (second.x - first.x) * (third.y - first.y) - (second.y - first.y) * (third.x - first.x);
}

enum class TurnType {
    COUNTER_CLOCKWISE,
    CLOCKWISE,
    COLLINEAR
};

TurnType GetTurnType(const Point& first, const Point& second, const Point& third) {
    int64_t crossProduct = CrossProduct(first, second, third);
    if(crossProduct > 0) {
        return TurnType::COUNTER_CLOCKWISE;
    } else if(crossProduct < 0) {
        return TurnType::CLOCKWISE;
    }
    return TurnType::COLLINEAR;
}

double GetLength(const Point& first, const Point& second) {
    return std::sqrt(std::pow((second.x - first.x), 2) + std::pow((second.y - first.y), 2));
}

class ConvexHullHandler {
private:
    std::vector<Point> points_;

    Point firstPoint_;
    Point secondPoint_;

    std::vector<Point> aboveLine_;
    std::vector<Point> belowLine_;

    void HandleLoop(bool above, const Point& point);

public:
    [[maybe_unused]] explicit ConvexHullHandler(const std::vector<Point> &points) : points_(points) {}
    explicit ConvexHullHandler(std::vector<Point>&& points) : points_(std::move(points)) {}

    void BuildConvexHull();
    double GetConvexHullLength() const;
};

void ConvexHullHandler::HandleLoop(bool above, const Point& point) {
    auto ConvexHullHandler::*current_vector = above ?
                                              &ConvexHullHandler::aboveLine_ : &ConvexHullHandler::belowLine_;

    auto condition = above ?
            [](const Point& first, const Point& second, const Point& third) {
                TurnType type = GetTurnType(first, second, third);
                return type == TurnType::CLOCKWISE || type == TurnType::COLLINEAR;
            } :
            [](const Point& first, const Point& second, const Point& third) {
                TurnType type = GetTurnType(first, second, third);
                return type == TurnType::COUNTER_CLOCKWISE || type == TurnType::COLLINEAR;
            };

    if(condition(firstPoint_, point, secondPoint_)) {

        while ((this->*current_vector).size() > 1 && !condition((this->*current_vector).rbegin()[1],
                                                                (this->*current_vector).rbegin()[0], point)) {
            (this->*current_vector).pop_back();
        }
        (this->*current_vector).push_back(point);
    }
}

void ConvexHullHandler::BuildConvexHull() {
    aboveLine_.resize(0);
    belowLine_.resize(0);

    std::sort(points_.begin(), points_.end(), [](const Point& left, const Point& right) {
        return left.x < right.x || (left.x == right.x && left.y < right.y);
    });

    firstPoint_  = points_.front();
    secondPoint_ = points_.back();

    aboveLine_.push_back(firstPoint_);
    belowLine_.push_back(firstPoint_);

    for (auto point : points_) {
        /* above line */
        HandleLoop(true, point);
        /* below line */
        HandleLoop(false, point);
    }
}

double ConvexHullHandler::GetConvexHullLength() const {
    double result = 0;

    for(size_t i = 0; i < aboveLine_.size() - 1 ; ++i) {
        result += GetLength(aboveLine_[i], aboveLine_[i + 1]);
    }

    result += GetLength(aboveLine_.back(), secondPoint_);
    result += GetLength(secondPoint_, belowLine_.back());

    for(size_t i = 0; i < belowLine_.size() - 1 ; ++i) {
        result += GetLength(belowLine_[i], belowLine_[i + 1]);
    }

    return result;
}

int main() {
    size_t number_of_points = 0;
    std::cin >> number_of_points;
    std::set<std::pair<PointType, PointType>> points_filter;

    for(size_t i = 0; i < number_of_points; ++i) {
        PointType x, y;
        std::cin >> x >> y;
        points_filter.insert({x, y});
    }

    std::vector<Point> points;
    for(auto [x,y] : points_filter) {
        points.push_back(Point{x,y});
    }

    ConvexHullHandler builder(std::move(points));
    builder.BuildConvexHull();
    double result = builder.GetConvexHullLength();

    std::cout << std::setprecision(15) << result;

    return 0;
}
