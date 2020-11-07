#include <algorithm>
#include <cmath>
#include <set>
#include <vector>
#include <iostream>
#include <iomanip>

//Структура точки
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

/* If return > 0 -> OAB -- counter-clockwise turn
 * if return == 0 -> collinear
 * if return < 0 -> OAB -- clockwise turn
 */
double CrossProduct(const Point& first, const Point& second, const Point& third) {
    return (second.x - first.x) * (third.y - first.y) - (second.y - first.y) * (third.x - first.x);
}

double GetLength(const Point& first, const Point& second) {
    return std::sqrt(std::pow((second.x - first.x), 2) + std::pow((second.y - first.y), 2));
}

class ConvexHullBuilder {
private:
    std::vector<Point> points;

    Point first_point;
    Point second_point;

    std::vector<Point> above_line;
    std::vector<Point> below_line;

    void HandleLoop(bool above, const Point& point);

public:
    [[maybe_unused]] explicit ConvexHullBuilder(const std::vector<Point> &points) : points(points) {}
    explicit ConvexHullBuilder(std::vector<Point>&& points) : points(std::move(points)) {}

    double GetConvexHullLength();
};

void ConvexHullBuilder::HandleLoop(bool above, const Point& point) {
    auto ConvexHullBuilder::*current_vector = above ?
            &ConvexHullBuilder::above_line : &ConvexHullBuilder::below_line;

    auto condition = above ?
            [](const Point& first, const Point& second, const Point& third)
                {return CrossProduct(first, second, third) <= 0;} :
            [](const Point& first, const Point& second, const Point& third) {
                return CrossProduct(first, second, third) >= 0;};

    if(condition(first_point, point, second_point)) {

        while ((this->*current_vector).size() > 1 && !condition((this->*current_vector).rbegin()[1],
                                                                (this->*current_vector).rbegin()[0], point)) {
            (this->*current_vector).pop_back();
        }
        (this->*current_vector).push_back(point);
    }
}

double ConvexHullBuilder::GetConvexHullLength() {
    std::sort(points.begin(), points.end(), [](const Point& left, const Point& right) {
        return left.x < right.x || (left.x == right.x && left.y < right.y);
    });

    first_point  = points.front();
    second_point = points.back();

    above_line.push_back(first_point);
    below_line.push_back(first_point);

    for (auto point : points) {
        /* above line */
        HandleLoop(true, point);
        /* below line */
        HandleLoop(false, point);
    }

    double result = 0;

    for(size_t i = 0; i < above_line.size() - 1 ;++i) {
        result += GetLength(above_line[i], above_line[i+1]);
    }

    result += GetLength(above_line.back(), second_point);
    result += GetLength(second_point, below_line.back());

    for(size_t i = 0; i < below_line.size() - 1 ;++i) {
        result += GetLength(below_line[i], below_line[i+1]);
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

    ConvexHullBuilder builder(std::move(points));
    double result = builder.GetConvexHullLength();

    std::cout << std::setprecision(15) << result;

    return 0;
}
