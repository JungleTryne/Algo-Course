#include <algorithm>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <exception>
#include <set>
#include <vector>

const long double eps = 1e-9;

//Структура точки
template<typename PointType>
struct Point {
    PointType x;
    PointType y;

    explicit Point(PointType x= 0, PointType y= 0) : x(x), y(y) {};
    Point(const Point<PointType>& other) = default;
    ~Point() = default;

    bool operator==(const Point<PointType>& other) const;
    bool operator!=(const Point<PointType>& other) const;
};


//Класс многоугольника
template<typename PointType>
class Polygon {
protected:
    std::vector<Point<PointType>> vertices;
public:
    explicit Polygon(const std::vector<Point<PointType>> &points);

    template<typename... T>
    explicit Polygon(const Point<PointType>& first, const T&... other);

    size_t GetVerticesNumber() const;
    long double area() const;

    Polygon GetMinkowskiSum(const Polygon& first_vector) const;
};

template<typename PointType>
bool Point<PointType>::operator==(const Point<PointType>& other) const {
    if constexpr (std::is_same<PointType, long double>::value) {
        return std::abs(this->x - other.x) <= eps && std::abs(this->y - other.y) <= eps;
    } else {
        return (this->x == other.x) && (this->y == other.y);
    }

}

template<typename PointType>
bool Point<PointType>::operator!=(const Point<PointType>& other) const {
    return !(*this==other);
}


//Вычитаем точки как вектора
template<typename PointType>
Point<PointType> operator-(const Point<PointType>& one, const Point<PointType>& two) {
    Point newPoint(one.x - two.x, one.y - two.y);
    return newPoint;
}

//Складываем точки как вектора
template<typename PointType>
Point<PointType> operator+(const Point<PointType>& one, const Point<PointType>& two) {
    Point newPoint(one.x + two.x, one.y + two.y);
    return newPoint;
}

//Умножение на скаляр справа
template<typename PointType>
Point<PointType> operator*(const Point<PointType>& one, PointType coefficient) {
    Point newPoint(one.x * coefficient, one.y*coefficient);
    return newPoint;
}

//Умножение на скаляр слева
template<typename PointType>
Point<PointType> operator*(PointType coefficient, Point<PointType>& one) {
    return one*coefficient;
}

template<typename PointType>
PointType GetDeterminant(PointType x11, PointType x12, PointType x21, PointType x22) {
    return x11 * x22 - x12 * x21;
}

//==================Polygon class methods========================//
//Используется формула площади Гаусса
template<typename PointType>
long double Polygon<PointType>::area() const {
    size_t pointer = 0;
    PointType area = 0;
    for (size_t i = 0; i < this->vertices.size(); ++i) {
        area += GetDeterminant(this->vertices[pointer].x,
                               this->vertices[pointer].y,
                               this->vertices[(pointer + 1) % this->vertices.size()].x,
                               this->vertices[(pointer + 1) % this->vertices.size()].y
        );
        pointer++;
    }
    return static_cast<long double>(std::abs(area))*0.5;
}

template<typename PointType>
Polygon<PointType>::Polygon(const std::vector<Point<PointType>>& points) {
    this->vertices = points;
}

template<typename PointType>
template<typename... T>
Polygon<PointType>::Polygon(const Point<PointType>& first, const T &... other) : Polygon(other...) {
    this->vertices.insert(this->vertices.begin(), first);
}

template<typename PointType>
Polygon<PointType> Polygon<PointType>::GetMinkowskiSum(const Polygon &other) const {

    std::vector<Point<PointType>> MinkowskiSum;

    std::vector<Point<PointType>> FirstSides;
    for(size_t i = 0; i < this->vertices.size(); ++i) {
        FirstSides.push_back(this->vertices[(i+1) % this->vertices.size()] - this->vertices[i]);
    }

    std::vector<Point<PointType>> SecondSides;
    for(size_t i = 0; i < other.vertices.size(); ++i) {
        SecondSides.push_back(other.vertices[(i+1) % other.vertices.size()] - other.vertices[i]);
    }

    auto Comparator = [](const Point<PointType>& first_vector, const Point<PointType>& second_vector) {
        auto first = first_vector;
        auto second = first_vector + second_vector;
        return (first.x) * (second.y) - (first.y) * (second.x) >= 0;
    };

    std::merge(FirstSides.begin(), FirstSides.end(), SecondSides.begin(), SecondSides.end(), std::back_inserter(MinkowskiSum), Comparator);

    std::vector<Point<PointType>> MinkowskiSumFiltered = {MinkowskiSum[0]};
    for(size_t i = 1; i < MinkowskiSum.size(); ++i) {
        auto first = MinkowskiSum[i-1];
        auto second = MinkowskiSum[i];
        if(GetDeterminant(first.x, first.y, second.x, second.y) == 0) {
            MinkowskiSumFiltered.back() = MinkowskiSumFiltered.back() + second;
            continue;
        }
        MinkowskiSumFiltered.push_back(second);
    }

    std::vector<Point<PointType>> MinkowskiVertices = {this->vertices[0] + other.vertices[0]};
    for(size_t i = 0; i < MinkowskiSumFiltered.size(); ++i) {
        MinkowskiVertices.push_back(MinkowskiSumFiltered[i] + MinkowskiVertices[i]);
    }

    if(MinkowskiVertices[0] == MinkowskiVertices.back()) {
        MinkowskiVertices.pop_back();
    }

    return Polygon(MinkowskiVertices);
}

template<typename PointType>
size_t Polygon<PointType>::GetVerticesNumber() const {
    return this->vertices.size();
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    std::vector<Polygon<int64_t>> polygons = {};
    for (size_t i = 0; i < 2; ++i) {
        size_t vertexNumber = 0;
        std::cin >> vertexNumber;
        std::vector<Point<int64_t>> vertices;
        for (size_t j = 0; j < vertexNumber; ++j) {
            int64_t x = 0;
            int64_t y = 0;
            std::cin >> x >> y;
            vertices.push_back(Point(x, y));
        }
        size_t min_point = 0;
        for(size_t j = 0; j < vertexNumber; ++j) {
            if(vertices[j].x < vertices[min_point].x ||
               (vertices[j].x == vertices[min_point].x && vertices[j].y < vertices[min_point].y))  {
                min_point = j;
            }
        }
        std::rotate(vertices.begin(), vertices.begin() + min_point, vertices.end());

        auto polygon = Polygon(vertices);
        polygons.push_back(polygon);
    }


    auto sum = polygons[0].GetMinkowskiSum(polygons[1]);
    std::cout << std::fixed << std::setprecision(6) <<
        static_cast<long double>((sum.area() - polygons[0].area() - polygons[1].area())) / 2 << std::endl;

    return 0;
}
