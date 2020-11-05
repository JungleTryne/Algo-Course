#include <algorithm>
#include <cmath>
#include <iostream>
#include <type_traits>
#include <iomanip>

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

constexpr double GetThirdDeterminant(std::array<std::array<double, 3>, 3> matrix) {
    return  matrix[0][0]*(matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1]) -
            matrix[0][1]*(matrix[1][0]*matrix[2][2] - matrix[1][2]*matrix[2][0]) +
            matrix[0][2]*(matrix[1][0]*matrix[2][1] - matrix[1][1]*matrix[2][0]);
}

class DistanceFinder {
private:
    const Segment<int64_t> first_;
    const Segment<int64_t> second_;

    Point<double> normalDeltaFirst;
    Point<double> normalDeltaSecond;
    Point<double> crossProduct;
    double denominator{};

    Point<double> closestFirst;
    Point<double> closestSecond;

    void HandleNotLineFirst();
    void HandleNotLineSecond();

    double HandleParallelSegments();
    double HandleCrossingSegments();

    double GetDistanceBetweenPointAndSegment(const Point<int64_t>& point, const Segment<int64_t> segment);

    auto GetCoefficients() const -> std::pair<double, double>;

public:
    explicit DistanceFinder(const Segment<int64_t>& first, const Segment<int64_t>& second) : first_(first),
        second_(second) {}

    double FindDistance();
};

double DistanceFinder::GetDistanceBetweenPointAndSegment(const Point<int64_t> &point, const Segment<int64_t> segment) {
    auto deltaPoint = point - segment.first;
    if(std::abs((segment.second - segment.first).GetLength()) < eps) {
        return (point - segment.second).GetLength();
    }
    auto normalDeltaSegment = GetNormal(segment.second - segment.first);
    double dotProduct = GetDotProduct(static_cast<Point<double>>(deltaPoint), normalDeltaSegment);
    if(dotProduct < 0) {
        dotProduct = 0;
    } else if(dotProduct > (segment.second - segment.first).GetLength()) {
        dotProduct = (segment.second - segment.first).GetLength();
    }
    auto ResultPoint = dotProduct * normalDeltaSegment + static_cast<Point<double>>(segment.first);
    return (ResultPoint - static_cast<Point<double>>(point)).GetLength();
}

double DistanceFinder::FindDistance() {
    auto deltaFirst  = first_.second - first_.first;
    auto deltaSecond = second_.second - second_.first;

    if(deltaFirst.GetLength() == 0) {
        return GetDistanceBetweenPointAndSegment(first_.first, second_);
    }
    if(deltaSecond.GetLength() == 0) {
        return GetDistanceBetweenPointAndSegment(second_.first, first_);
    }

    normalDeltaFirst  = GetNormal(deltaFirst);
    normalDeltaSecond = GetNormal(deltaSecond);

    crossProduct = GetCrossProduct(normalDeltaFirst, normalDeltaSecond);

    denominator = crossProduct.GetLength();
    denominator = denominator*denominator;

    if(denominator == 0) {
        auto result = HandleParallelSegments();
        return result;
    }

    return HandleCrossingSegments();
}

double DistanceFinder::HandleParallelSegments() {
    if((first_.first - first_.second).GetLength() == 0) {

    }

    double firstEdgeDot  = GetDotProduct(normalDeltaFirst, static_cast<Point<double>>(second_.first - first_.first));
    double secondEdgeDot = GetDotProduct(normalDeltaFirst, static_cast<Point<double>>(second_.second - first_.first));

    if(firstEdgeDot < 0 && secondEdgeDot < 0) {
        return std::abs(firstEdgeDot) < std::abs(secondEdgeDot) ?
               (second_.first - first_.first).GetLength() : (second_.second - first_.first).GetLength();
    }

    if(firstEdgeDot > first_.GetLength() && secondEdgeDot > first_.GetLength()) {
        return std::abs(firstEdgeDot) < std::abs(secondEdgeDot) ?
               (second_.first - first_.second).GetLength() : (second_.second - first_.second).GetLength();
    }

    return ((firstEdgeDot * normalDeltaFirst) + static_cast<Point<double>>(first_.first) -
        static_cast<Point<double>>(second_.first)).GetLength();
}

auto DistanceFinder::GetCoefficients() const -> std::pair<double, double>
{
    auto edgeDelta = second_.first - first_.first;

    double DetFirst = GetThirdDeterminant({{
           {static_cast<double>(edgeDelta.x), static_cast<double>(edgeDelta.y), static_cast<double>(edgeDelta.z)},
           {normalDeltaSecond.x, normalDeltaSecond.y, normalDeltaSecond.z},
           {crossProduct.x, crossProduct.y, crossProduct.z}
   }});

    double DetSecond = GetThirdDeterminant({{
            {static_cast<double>(edgeDelta.x), static_cast<double>(edgeDelta.y), static_cast<double>(edgeDelta.z)},
            {normalDeltaFirst.x, normalDeltaFirst.y, normalDeltaFirst.z},
            {crossProduct.x, crossProduct.y, crossProduct.z}
    }});

    double coeffFirst  = DetFirst / denominator;
    double coeffSecond = DetSecond / denominator;
    return std::make_pair(coeffFirst, coeffSecond);
}

double DistanceFinder::HandleCrossingSegments()
{
    auto [coeffFirst, coeffSecond] = GetCoefficients();

    closestFirst  = static_cast<Point<double>>(first_.first) + (coeffFirst * normalDeltaFirst);
    closestSecond = static_cast<Point<double>>(second_.first) + (coeffSecond * normalDeltaSecond);

    bool notLineShortestFirst = false;
    bool notLineShortestSecond = false;

    if(coeffFirst < 0 || coeffFirst > first_.GetLength()) {
        notLineShortestFirst = true;
        closestFirst = (coeffFirst < 0) ? static_cast<Point<double>>(first_.first) :
                       static_cast<Point<double>>(first_.second);
    }

    if(coeffSecond < 0 || coeffSecond > second_.GetLength()) {
        notLineShortestSecond = true;
        closestSecond = (coeffSecond < 0) ? static_cast<Point<double>>(second_.first) :
                        static_cast<Point<double>>(second_.second);
    }

    if(notLineShortestFirst) {
        HandleNotLineFirst();
    }

    if(notLineShortestSecond) {
        HandleNotLineSecond();
    }

    return (closestFirst - closestSecond).GetLength();
}

void DistanceFinder::HandleNotLineFirst() {
    double dotProduct = GetDotProduct(normalDeltaSecond,
                                      closestFirst - static_cast<Point<double>>(second_.first));
    if(dotProduct < 0) {
        dotProduct = 0;
    } else if (dotProduct > second_.GetLength()) {
        dotProduct = second_.GetLength();
    }
    closestSecond = static_cast<Point<double>>(second_.first) + (normalDeltaSecond * dotProduct);
}

void DistanceFinder::HandleNotLineSecond() {
    double dotProduct = GetDotProduct(normalDeltaFirst,
                                      closestSecond - static_cast<Point<double>>(first_.first));
    if(dotProduct < 0) {
        dotProduct = 0;
    } else if (dotProduct > first_.GetLength()) {
        dotProduct = first_.GetLength();
    }
    closestFirst = static_cast<Point<double>>(first_.first) + (normalDeltaFirst * dotProduct);
}

int main() {
    int64_t x = 0;
    int64_t y = 0;
    int64_t z = 0;
    std::cin >> x >> y >> z;

    Point<int64_t> first_point = Point(x,y,z);
    std::cin >> x >> y >> z;
    Point<int64_t> second_point = Point(x,y,z);

    Segment<int64_t> first_segment = Segment(first_point, second_point);

    std::cin >> x >> y >> z;
    first_point = Point(x,y,z);
    std::cin >> x >> y >> z;
    second_point = Point(x,y,z);

    Segment<int64_t> second_segment = Segment(first_point, second_point);

    DistanceFinder finder(first_segment, second_segment);
    std::cout << std::fixed << std::setprecision(8) << finder.FindDistance() << std::endl;

    return 0;
}
