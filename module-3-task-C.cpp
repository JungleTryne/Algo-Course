#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <type_traits>
#include <set>
#include <queue>

const double eps = 1e-9;

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

    double GetLength() const;
};


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

//Вычитаем точки как вектора
template<typename PointType>
Point<PointType> operator-(const Point<PointType>& one, const Point<PointType>& two) {
    Point newPoint(one.x - two.x, one.y - two.y, one.z - two.z);
    return newPoint;
}

//Складываем точки как вектора
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
    double GetLength() const;
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

template <typename PointType>
struct Edge {
    Point<PointType> first;
    Point<PointType> second;
    explicit Edge(const Point<PointType>& first, const Point<PointType>& second) :
            first(first), second(second) {}
};

class ConvexHullBuilder {
private:
    std::vector<Hull<int64_t>> hulls{};
    std::vector<Point<int64_t>> points{};

    Hull<int64_t> GetFirstHull() const;

    double GetHullDegree(const Point<int64_t> &common_first, const Point<int64_t> &common_second,
                         const Point<int64_t> &first, const Point<int64_t> &second) const;
    bool IsUsed(const Edge<int64_t>& edge) const;
    void Use(const Edge<int64_t>& edge);
public:
    explicit ConvexHullBuilder(const std::vector<Point<int64_t>>& points) : points(points) {}
    void BuildHull();
};

Hull<int64_t> ConvexHullBuilder::GetFirstHull() const {
    Point<int64_t> firstMinPoint  = points[0];
    Point<int64_t> secondMinPoint = points[0];
    Point<int64_t> thirdMinPoint  = points[0];

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
    return Hull<int64_t>(firstMinPoint, secondMinPoint, thirdMinPoint);
}

void ConvexHullBuilder::BuildHull() {
    auto firstHull = GetFirstHull();
    hulls.push_back(firstHull);

    std::queue<std::pair<Edge<int64_t>, Hull<int64_t>>> edge_queue;
    edge_queue.push(std::make_pair(Edge{firstHull.first, firstHull.second}, firstHull));

    while (!edge_queue.empty()) {
        auto [edge, hull] = edge_queue.front();
        auto thirdVertex = hull.GetOtherPoint(edge.first, edge.second);

        edge_queue.pop();

        if(IsUsed(edge)) {
            continue;
        }

        double bestAngle = -1;
        Hull bestHull = hull;

        for(auto point : points) {
            if(point != edge.first && point != edge.second) {
                auto newHull = Hull(edge.first, edge.second, point);

                double angle = GetHullDegree(edge.first, edge.second, thirdVertex, point);
                if(angle > bestAngle) {
                    bestAngle = angle;
                    bestHull = newHull;
                }
            }
        }

        hulls.push_back(bestHull);
        edge_queue.push(std::make_pair(Edge{bestHull.first, bestHull.second}, bestHull));
        edge_queue.push(std::make_pair(Edge{bestHull.first, bestHull.third}, bestHull));
        edge_queue.push(std::make_pair(Edge{bestHull.second, bestHull.third}, bestHull));
        Use(edge);
    }
}

double ConvexHullBuilder::GetHullDegree(const Point<int64_t> &common_first, const Point<int64_t> &common_second,
                                        const Point<int64_t> &first, const Point<int64_t> &second) const {
    /* It's definite that first and second have common edge! */
    Point<int64_t> first_vector  = common_second - common_first;
    Point<int64_t> second_vector = common_second - first;
    Point<int64_t> third_vector  = common_second - second;

    auto scalar_second = GetDotProduct(first_vector, second_vector);
    auto scalar_third  = GetDotProduct(first_vector, third_vector);

    Point<double> orthoSecond = static_cast<Point<double>>(second_vector) -
        scalar_second*static_cast<Point<double>>(first_vector)*(1/first_vector.GetLength());

    Point<double> orthoThird = static_cast<Point<double>>(third_vector) -
        scalar_third*static_cast<Point<double>>(first_vector)*(1/first_vector.GetLength());

    auto scalar_angle = GetDotProduct(orthoThird, orthoSecond);
    return acos(scalar_angle/(orthoSecond.GetLength()*orthoThird.GetLength()));
}

int main() {

}