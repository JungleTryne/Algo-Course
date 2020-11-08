#include <bits/stdc++.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <type_traits>
#include <unordered_set>
#include <queue>

const double eps = 1e-9;
const int32_t max_coordinate = 10000;

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
    explicit operator Point<double>() const;

    [[nodiscard]] double GetLength() const;
};

namespace std {
    template <> struct hash<Point<int32_t>> {
        size_t operator()(const Point<int32_t>& point) const {
            size_t x_hash = point.x;
            size_t y_hash = 31*point.y;
            size_t z_hash = 127*point.z;
            return x_hash ^ y_hash ^ z_hash;
        }
    };
}


template<typename PointType>
bool Point<PointType>::operator==(const Point &other) const {
    if constexpr (std::is_same<PointType, double>::value) {
        return  std::abs(this->x - other.x) <= eps &&
                std::abs(this->y - other.y) <= eps &&
                std::abs(this->z - other.z) <= eps;
    } else {
        return this->x == other.x && this->y == other.y && this->z == other.z;
    }

}

template<typename PointType>
bool Point<PointType>::operator!=(const Point<PointType> &other) const {
    return !(*this == other);
}

template<typename PointType>
double Point<PointType>::GetLength() const {
    return sqrt(x*x + y*y + z*z);
}

template<typename PointType>
Point<PointType>::operator Point<double>() const {
    return Point<double>(static_cast<double>(x), static_cast<double>(y), static_cast<double>(z));
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
Point<PointType> operator*(double scalar, const Point<PointType>& one) {
    Point newPoint(one.x*scalar, one.y*scalar, one.z*scalar);
    return newPoint;
}

template<typename PointType>
PointType GetDotProduct(const Point<PointType>& one, const Point<PointType>& two) {
    return one.x*two.x + one.y*two.y + one.z*two.z;
}

template<typename PointType>
struct Segment {
    Point<PointType> first;
    Point<PointType> second;
    explicit Segment(Point<PointType> first, Point<PointType> second) : first(first), second(second) {}
    bool operator==(const Segment<PointType>& other) {
        return this->first == other.first && this->second == other.second;
    }
    [[nodiscard]] double GetLength() const;
};

template<typename PointType>
double Segment<PointType>::GetLength() const {
    return (first - second).GetLength();
}

template <typename PointType>
Point<double> GetNormal(const Point<PointType>& vector) {
    double length = vector.GetLength();
    double new_x = static_cast<double>(vector.x) / length;
    double new_y = static_cast<double>(vector.y) / length;
    double new_z = static_cast<double>(vector.z) / length;
    return Point<double>(new_x, new_y, new_z);
}

template <typename PointType>
Point<PointType> GetCrossProduct(const Point<PointType>& first, const Point<PointType>& second) {
    return Point<PointType>(
            first.y*second.z - first.z*second.y,
            -(first.x*second.z - first.z*second.x),
            first.x*second.y - first.y*second.x
    );
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
    std::unordered_set<Point<PointType>> first_points {first, second, third};
    std::unordered_set<Point<PointType>> second_points {other.first, other.second, other.third};
    return std::equal(first_points.begin(), first_points.end, second_points.begin());
}

namespace std {
    template<> struct hash<Hull<int32_t>> {
        size_t operator()(const Hull<int32_t>& hull) const {
            size_t first_hash = std::hash<Point<int32_t>>()(hull.first);
            size_t second_hash = std::hash<Point<int32_t>>()(hull.second);
            size_t third_hash = std::hash<Point<int32_t>>()(hull.third);
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
    template<> struct hash<Edge<int32_t>> {
        size_t operator()(const Edge<int32_t>& edge) const {
            size_t hash_first = std::hash<Point<int32_t>>()(edge.first);
            size_t hash_second = std::hash<Point<int32_t>>()(edge.second);
            return hash_first ^ hash_second;
        }
    };
}

class ConvexHullBuilder {
private:
    std::vector<Hull<int32_t>> hulls{};
    std::vector<Point<int32_t>> points{};

    std::unordered_set<Edge<int32_t>> used_edges_;

    [[nodiscard]] Hull<int32_t> GetFirstHull() const;

    [[nodiscard]] double GetHullDegree(const Point<int32_t> &common_first, const Point<int32_t> &common_second,
                         const Point<int32_t> &first, const Point<int32_t> &second) const;
    bool IsUsed(const Edge<int32_t>& edge) const;
    void Use(const Edge<int32_t>& edge);
public:
    explicit ConvexHullBuilder(const std::vector<Point<int32_t>>& points) : points(points) {}
    void BuildHull();
    double GetShortestEscapeDistance(const Point<int32_t>& point) const;
};

Hull<int32_t> ConvexHullBuilder::GetFirstHull() const {
    Point<int32_t> firstMinPoint  = Point(max_coordinate, max_coordinate, max_coordinate);
    Point<int32_t> secondMinPoint = Point(max_coordinate, max_coordinate, max_coordinate);
    Point<int32_t> thirdMinPoint  = Point(max_coordinate, max_coordinate, max_coordinate);

    for(auto point : points) {
        if(point.y < firstMinPoint.y || (point.y == firstMinPoint.y && point.x < firstMinPoint.x)) {
            firstMinPoint = point;
        }
    }
    for(auto point : points) {
        if (point != firstMinPoint) {
            if (point.y < secondMinPoint.y || (point.y == secondMinPoint.y && point.x < secondMinPoint.x)) {
                secondMinPoint = point;
            }
        }
    }

    for(auto point : points) {
        if (point != firstMinPoint && point != secondMinPoint) {
            if (point.y < thirdMinPoint.y || (point.y == thirdMinPoint.y && point.x < thirdMinPoint.x)) {
                thirdMinPoint = point;
            }
        }
    }
    return Hull<int32_t>(firstMinPoint, secondMinPoint, thirdMinPoint);
}

void ConvexHullBuilder::BuildHull() {
    auto firstHull = GetFirstHull();
    hulls.push_back(firstHull);

    std::queue<std::pair<Edge<int32_t>, Hull<int32_t>>> edge_queue;
    edge_queue.push(std::make_pair(Edge{firstHull.first, firstHull.second}, firstHull));

    while (!edge_queue.empty()) {
        auto [edge, hull] = edge_queue.front();
        auto thirdVertex = hull.GetOtherPoint(edge.first, edge.second);

        edge_queue.pop();

        if(IsUsed(edge)) {
            continue;
        }

        double bestAngle = -1;
        Point bestPoint = thirdVertex;

        for(auto point : points) {
            if(point != edge.first && point != edge.second) {
                double angle = GetHullDegree(edge.first, edge.second, thirdVertex, point);
                if(angle > bestAngle) {
                    bestAngle = angle;
                    bestPoint = point;
                }
            }
        }
        if(!(IsUsed(Edge{edge.first, bestPoint}) && IsUsed(Edge{edge.second, bestPoint}))) {
            auto bestHull = Hull{edge.first, edge.second, bestPoint};
            hulls.push_back(bestHull);
            edge_queue.push(std::make_pair(Edge{bestHull.second, bestHull.third}, bestHull));
            edge_queue.push(std::make_pair(Edge{bestHull.third, bestHull.first}, bestHull));
        }
        Use(edge);
    }
}

double ConvexHullBuilder::GetHullDegree(const Point<int32_t> &common_first, const Point<int32_t> &common_second,
                                        const Point<int32_t> &first, const Point<int32_t> &second) const {
    /* It's definite that first and second have common edge! */
    Point<int32_t> first_vector  = common_second - common_first;
    Point<int32_t> second_vector = first - common_first;
    Point<int32_t> third_vector  = second - common_first;

    auto scalar_second = GetDotProduct(first_vector, second_vector);
    auto scalar_third  = GetDotProduct(first_vector, third_vector);

    Point<double> orthoSecond = static_cast<Point<double>>(second_vector) -
        scalar_second*static_cast<Point<double>>(first_vector)*(1/first_vector.GetLength());

    Point<double> orthoThird = static_cast<Point<double>>(third_vector) -
        scalar_third*static_cast<Point<double>>(first_vector)*(1/first_vector.GetLength());

    auto scalar_angle = GetDotProduct(orthoThird, orthoSecond);
    return acos(scalar_angle/(orthoSecond.GetLength()*orthoThird.GetLength()));
}

bool ConvexHullBuilder::IsUsed(const Edge<int32_t> &edge) const {
    return used_edges_.find(edge) != used_edges_.end();
}

void ConvexHullBuilder::Use(const Edge<int32_t> &edge) {
    used_edges_.insert(edge);
}

double GetDistanceBetweenPlaneAndPoint(const Point<int32_t>& first, const Point<int32_t>& second,
                                       const Point<int32_t>& third, const Point<int32_t>& fourth) {

}

double ConvexHullBuilder::GetShortestEscapeDistance(const Point<int32_t> &point) const {
    for(auto hull : hulls) {

    }
}

int main() {
    size_t number_of_points = 0;
    std::cin >> number_of_points;

    std::vector<Point<int32_t>> points;
    for(size_t i = 0; i < number_of_points; ++i) {
        int32_t x = 0;
        int32_t y = 0;
        int32_t z = 0;
        std::cin >> x >> y >> z;
        points.emplace_back(x, y, z);
    }

    ConvexHullBuilder builder(points);
    builder.BuildHull();
    return 0;
}