#include <bits/stdc++.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <type_traits>
#include <unordered_set>
#include <queue>

long double eps() {
    return 1e-12;
}

int64_t max_coordinate() {
    return 100000;
}

long double pi() {
    return 3.14159265358979323846;
}

template<typename PointType>
struct Point {
    PointType x;
    PointType y;
    PointType z;

    explicit Point(PointType x=0, PointType y=0, PointType z=0) : x(x), y(y), z(z) {};
    Point(const Point& other) = default;
    ~Point() = default;

    bool operator==(const Point& other) const;
    bool operator!=(const Point& other) const;
    bool operator<(const Point& other) const;
    explicit operator Point<long double>() const;

    [[nodiscard]] long double GetLength() const;
};

namespace std {
    template <> struct hash<Point<int64_t>> {
        size_t operator()(const Point<int64_t>& point) const {
            size_t x_hash = static_cast<size_t>(point.x);
            size_t y_hash = 31*static_cast<size_t>(point.y);
            size_t z_hash = 127*static_cast<size_t>(point.z);
            return x_hash ^ y_hash ^ z_hash;
        }
    };
}


template<typename PointType>
bool Point<PointType>::operator==(const Point &other) const {
    if constexpr (std::is_same<PointType, long double>::value) {
        return  std::abs(this->x - other.x) <= eps() &&
                std::abs(this->y - other.y) <= eps() &&
                std::abs(this->z - other.z) <= eps();
    } else {
        return this->x == other.x && this->y == other.y && this->z == other.z;
    }

}

template<typename PointType>
bool Point<PointType>::operator!=(const Point<PointType> &other) const {
    return !(*this == other);
}

template<typename PointType>
long double Point<PointType>::GetLength() const {
    return sqrt(x*x + y*y + z*z);
}

template<typename PointType>
Point<PointType>::operator Point<long double>() const {
    return Point<long double>(static_cast<long double>(x), static_cast<long double>(y), static_cast<long double>(z));
}

template<typename PointType>
bool Point<PointType>::operator<(const Point &other) const {
    return std::tuple(x, y, z) < std::tuple(other.x, other.y, other.z);
}

template<typename PointType>
Point<PointType> operator-(const Point<PointType>& one, const Point<PointType>& two) {
    Point newPoint(one.x - two.x, one.y - two.y, one.z - two.z);
    return newPoint;
}

template<typename PointType>
Point<PointType> operator+(const Point<PointType>& one, const Point<PointType>& two) {
    Point newPoint(one.x + two.x, one.y + two.y, one.z + two.z);
    return newPoint;
}

template<typename PointType>
Point<PointType> operator*(const Point<PointType>& one, PointType scalar) {
    Point newPoint(one.x*scalar, one.y*scalar, one.z*scalar);
    return newPoint;
}

template<typename PointType>
Point<PointType> operator*(long double scalar, const Point<PointType>& one) {
    Point newPoint(one.x*scalar, one.y*scalar, one.z*scalar);
    return newPoint;
}

template<typename PointType>
PointType GetDotProduct(const Point<PointType>& one, const Point<PointType>& two) {
    return one.x*two.x + one.y*two.y + one.z*two.z;
}

template <typename PointType>
Point<PointType> GetCrossProduct(const Point<PointType>& first, const Point<PointType>& second) {
    return Point<PointType>(
            first.y*second.z - first.z*second.y,
            -(first.x*second.z - first.z*second.x),
            first.x*second.y - first.y*second.x
    );
}

template<typename PointType>
constexpr long double GetThirdDeterminant(std::array<std::array<PointType, 3>, 3> matrix) {
    return  matrix[0][0]*(matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1]) -
            matrix[0][1]*(matrix[1][0]*matrix[2][2] - matrix[1][2]*matrix[2][0]) +
            matrix[0][2]*(matrix[1][0]*matrix[2][1] - matrix[1][1]*matrix[2][0]);
}

template <typename PointType>
struct Hull {
    Point<PointType> first;
    Point<PointType> second;
    Point<PointType> third;
    explicit Hull(const Point<PointType>& first, const Point<PointType>& second, const Point<PointType>& third) :
            first(first), second(second), third(third) {}

    Point<PointType> GetOtherPoint(const Point<PointType>& edge_first, const Point<PointType>& edge_second) const;
    bool operator==(const Hull& other) const;
};

template<typename PointType>
Point<PointType> Hull<PointType>::GetOtherPoint(const Point<PointType> &edge_first,
                                                const Point<PointType> &edge_second) const {
    if(first == edge_first && second == edge_second) return third;
    if(first == edge_second && second == edge_first) return third;
    if(first == edge_first && third == edge_second) return second;
    if(first == edge_second && third == edge_first) return second;
    if(second == edge_first && third == edge_second) return first;
    if(second == edge_second && third == edge_first) return first;

    assert(false);
}

template<typename PointType>
bool Hull<PointType>::operator==(const Hull &other) const {
    std::vector<Point<PointType>> sorted = {first, second, third};
    std::vector<Point<PointType>> sorted_other = {other.first, other.second, other.third};
    std::sort(sorted.begin(), sorted.end());
    std::sort(sorted_other.begin(), sorted_other.end());
    return sorted == sorted_other;
}

namespace std {
    template<> struct hash<Hull<int64_t>> {
        size_t operator()(const Hull<int64_t>& hull) const {
            size_t first_hash = std::hash<Point<int64_t>>()(hull.first);
            size_t second_hash = std::hash<Point<int64_t>>()(hull.second);
            size_t third_hash = std::hash<Point<int64_t>>()(hull.third);
            return first_hash ^ second_hash ^ third_hash;
        }
    };
}

template <typename PointType>
struct Edge {
    Point<PointType> first;
    Point<PointType> second;
    explicit Edge(const Point<PointType>& first, const Point<PointType>& second) :
            first(first), second(second) {}
    bool operator==(const Edge<PointType>& other) const;
};

template<typename PointType>
bool Edge<PointType>::operator==(const Edge<PointType> &other) const {
    return (first == other.first && second == other.second) || (first == other.second && second == other.first);
}

namespace std {
    template<> struct hash<Edge<int64_t>> {
        size_t operator()(const Edge<int64_t>& edge) const {
            size_t hash_first = std::hash<Point<int64_t>>()(edge.first);
            size_t hash_second = std::hash<Point<int64_t>>()(edge.second);
            return hash_first ^ hash_second;
        }
    };
}

long double GetOXYAngle(const Point<int64_t>& point) {
    auto orthoPoint = point;
    orthoPoint.z = 0;
    return acosl(orthoPoint.GetLength() / point.GetLength());
}

long double GetDistanceBetweenPointAndHull(const Point<int64_t>& point, const Hull<int64_t>& hull) {
    const auto& [first, second, third] = hull;
    auto numerator = GetThirdDeterminant<int64_t>({{
                                                           {point.x - first.x, point.y - first.y, point.z - first.z},
                                                           {second.x - first.x, second.y - first.y, second.z - first.z},
                                                           {third.x - first.x, third.y - first.y, third.z - first.z}
                                                   }});
    auto denominator = GetCrossProduct(first - second, first - third).GetLength();
    numerator = std::abs(numerator);
    return static_cast<long double>(numerator) / static_cast<long double>(denominator);
}

class ConvexHullBuilder {
private:
    std::unordered_set<Hull<int64_t>> hulls{};
    std::vector<Point<int64_t>> points{};

    std::unordered_set<Edge<int64_t>> used_edges_;

    [[nodiscard]] Hull<int64_t> GetFirstHull() const;

    [[nodiscard]] long double GetHullDegree(const Point<int64_t> &commonFirst, const Point<int64_t> &commonSecond,
                                            const Point<int64_t> &first, const Point<int64_t> &second) const;
    bool IsUsed(const Edge<int64_t>& edge) const;
    void Use(const Edge<int64_t>& edge);

    long double GetHullOXYDegree(const Point<int64_t> &common_first, const Point<int64_t> &common_second,
                                 const Point<int64_t> &third) const;
public:
    explicit ConvexHullBuilder(const std::vector<Point<int64_t>>& points) : points(points) {}

    void BuildHull();
    long double GetMinDistance(const Point<int64_t>& point) const;
};

Hull<int64_t> ConvexHullBuilder::GetFirstHull() const {
    Point<int64_t> firstMinPoint  = Point(max_coordinate(), max_coordinate(), max_coordinate());

    for(auto point : points) {
        if(point.z < firstMinPoint.z) {
            firstMinPoint = point;
        }
    }

    double minAngle = 10;

    Point<int64_t> secondMinPoint = Point(max_coordinate(), max_coordinate(), max_coordinate());
    for(auto point : points) {
        if(point != firstMinPoint) {
            double angle = GetOXYAngle(point - firstMinPoint);
            if(angle < minAngle) {
                minAngle = angle;
                secondMinPoint = point;
            }
        }
    }

    minAngle = 10;

    Point<int64_t> thirdMinPoint = Point(max_coordinate(), max_coordinate(), max_coordinate());
    for(auto point : points) {
        if(point != firstMinPoint && point != secondMinPoint) {
            double angle = GetHullOXYDegree(firstMinPoint, secondMinPoint, point);
            if(angle < minAngle) {
                minAngle = angle;
                thirdMinPoint = point;
            }
        }
    }

    return Hull<int64_t>(firstMinPoint, secondMinPoint, thirdMinPoint);
}

void ConvexHullBuilder::BuildHull() {
    auto firstHull = GetFirstHull();
    hulls.insert(firstHull);

    std::queue<std::pair<Edge<int64_t>, Hull<int64_t>>> edge_queue;
    edge_queue.push(std::make_pair(Edge{firstHull.first, firstHull.second}, firstHull));

    while (!edge_queue.empty()) {
        auto [edge, hull] = edge_queue.front();
        auto thirdVertex = hull.GetOtherPoint(edge.first, edge.second);

        edge_queue.pop();

        if(IsUsed(edge)) {
            continue;
        }

        long double bestAngle = -1;
        Point bestPoint = thirdVertex;

        for(auto point : points) {
            if(point != edge.first && point != edge.second) {
                long double angle = GetHullDegree(edge.first, edge.second, thirdVertex, point);
                if(angle > bestAngle) {
                    bestAngle = angle;
                    bestPoint = point;
                }
            }
        }
        if(!(IsUsed(Edge{edge.first, bestPoint}) && IsUsed(Edge{edge.second, bestPoint}))) {
            auto bestHull = Hull{edge.first, edge.second, bestPoint};
            hulls.insert(bestHull);
            edge_queue.push(std::make_pair(Edge{bestHull.second, bestHull.third}, bestHull));
            edge_queue.push(std::make_pair(Edge{bestHull.third, bestHull.first}, bestHull));
        }
        Use(edge);
    }
}

long double ConvexHullBuilder::GetHullDegree(const Point<int64_t> &commonFirst, const Point<int64_t> &commonSecond,
                                             const Point<int64_t> &first, const Point<int64_t> &second) const {
    /* It's definite that first and second have common edge! */
    Point<int64_t> firstVector  = commonSecond - commonFirst;
    Point<int64_t> secondVector = first - commonFirst;
    Point<int64_t> thirdVector  = second - commonFirst;

    auto scalar_second = GetDotProduct(firstVector, secondVector);
    auto scalar_third  = GetDotProduct(firstVector, thirdVector);

    Point<long double> orthoSecond = static_cast<Point<long double>>(secondVector) -
                                     scalar_second * static_cast<Point<long double>>(firstVector) *
                                     (1 / (firstVector.GetLength() * firstVector.GetLength()));

    Point<long double> orthoThird = static_cast<Point<long double>>(thirdVector) -
                                    scalar_third * static_cast<Point<long double>>(firstVector) *
                                    (1 / (firstVector.GetLength() * firstVector.GetLength()));

    auto scalarAngle = GetDotProduct(orthoThird, orthoSecond);
    return acosl(scalarAngle / (orthoSecond.GetLength() * orthoThird.GetLength()));
}

bool ConvexHullBuilder::IsUsed(const Edge<int64_t> &edge) const {
    return used_edges_.find(edge) != used_edges_.end();
}

void ConvexHullBuilder::Use(const Edge<int64_t> &edge) {
    used_edges_.insert(edge);
}

long double ConvexHullBuilder::GetMinDistance(const Point<int64_t> &point) const {
    long double min_distance = GetDistanceBetweenPointAndHull(point, *hulls.begin());
    for(const auto& hull : hulls) {
        long double distance = GetDistanceBetweenPointAndHull(point, hull);
        min_distance = std::min(min_distance, distance);
    }
    return min_distance;
}

long double ConvexHullBuilder::GetHullOXYDegree(const Point<int64_t> &first, const Point<int64_t> &second,
                                                const Point<int64_t> &third) const {
    auto firstVector = static_cast<Point<long double>>(second - first);
    auto secondVector = static_cast<Point<long double>>(third - first);
    firstVector = firstVector * (1 / firstVector.GetLength());
    secondVector = secondVector * (1 / secondVector.GetLength());
    auto crossProduct = GetCrossProduct(firstVector, secondVector);
    auto orthoCross = crossProduct;
    orthoCross.z = 0;
    return pi() / 2 - acosl(orthoCross.GetLength() / crossProduct.GetLength());
}

int main() {
    size_t number_of_points = 0;
    std::cin >> number_of_points;

    std::vector<Point<int64_t>> points;
    for(size_t i = 0; i < number_of_points; ++i) {
        int64_t x = 0;
        int64_t y = 0;
        int64_t z = 0;
        std::cin >> x >> y >> z;
        points.emplace_back(x, y, z);
    }

    ConvexHullBuilder builder(points);
    builder.BuildHull();

    size_t requests_number = 0;
    std::cin >> requests_number;
    for(size_t i = 0; i < requests_number; ++i) {
        Point<int64_t> point;
        std::cin >> point.x >> point.y >> point.z;
        std::cout << std::fixed << std::setprecision(4) << builder.GetMinDistance(point) << std::endl;
    }

    return 0;
}