#include <algorithm>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <complex>
#include <exception>
#include <set>
#include <vector>

const double eps = 1e-9;
const double pi = 3.14159265358979323;

//Структура точки
struct Point {
    double x;
    double y;

    explicit Point(double x=0, double y=0) : x(x), y(y) {};
    Point(const Point& other) = default;
    ~Point() = default;

    bool operator==(const Point& other) const;
    bool operator!=(const Point& other) const;

    double getLength() const;
    double getScalar(const Point& other) const;
    void rotate(double angle);
};

//Класс прямой
class Line {
private:
    Point one;
    Point two;
public:
    Line(const Point& first, const Point& second);
    Line(double angle, double height);
    Line(const Point& point, double angle);

    bool operator==(const Line& other) const;
    bool operator!=(const Line& other) const;

    Point getNormalVector() const;
};

//Абстрактный класс фигуры
class Shape {
public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool containsPoint(const Point& point) const = 0;

    virtual bool operator==(const Shape& other) const = 0;
    virtual bool operator!=(const Shape& other) const = 0;
    virtual bool isCongruentTo(const Shape& other) const = 0;
    virtual bool isSimilarTo(const Shape& other) const = 0;

    virtual ~Shape() = default;

    virtual void rotate(const Point& center, double angle) = 0;
    virtual void reflex(const Line& axis) = 0;
    virtual void scale(const Point& center, double coefficient) = 0;
};

//Класс многоугольника
class Polygon : public Shape {
protected:
    std::vector<Point> vertices;
    static bool verticesVectorsEqual(std::vector<Point>& thisCopy, std::vector<Point>& otherCopy) ;
    static bool scalarAndEdgesEqual(std::vector<double>& scalarThis, std::vector<double>& scalarOther,
                                    std::vector<double>& edgesThis, std::vector<double>& edgesOther);
public:
    Polygon(const Polygon& other);
    explicit Polygon(const std::vector<Point>& points);
    explicit Polygon(const Point& first);

    template<typename... T>
    explicit Polygon(const Point& first, const T&... other);

    double perimeter() const override;
    double area() const override;
    bool operator==(const Shape& other) const override;
    bool operator!=(const Shape& other) const override;
    bool isCongruentTo(const Shape& other) const override;
    bool isSimilarTo(const Shape& other) const override;
    bool containsPoint(const Point& point) const override;
    bool isConvex() const;
    std::vector<Point> getVertices() const;

    void rotate(const Point& center, double angle) override;
    void reflex(const Line& axis) override;
    void scale(const Point& center, double coefficient) override;

    Polygon GetMinkowskiSum(const Polygon& other) const;
};


bool Point::operator==(const Point &other) const {
    return std::abs(this->x - other.x) <= eps && std::abs(this->y - other.y) <= eps;
}

bool Point::operator!=(const Point &other) const {
    return !(*this==other);
}

double Point::getLength() const {
    return sqrt(this->x * this->x + this->y * this->y);
}

void Point::rotate(double angle) {
    //angle в радианах
    double newX = this->x*cos(angle) - this->y*sin(angle);
    double newY = this->x*sin(angle) + this->y*cos(angle);
    this->x = newX;
    this->y = newY;
}

double Point::getScalar(const Point& other) const {
    return this->x*other.x + this->y*other.y;
}

//Вычитаем точки как вектора
Point operator-(const Point& one, const Point& two) {
    Point newPoint(one.x - two.x, one.y - two.y);
    return newPoint;
}

//Складываем точки как вектора
Point operator+(const Point& one, const Point& two) {
    Point newPoint(one.x + two.x, one.y + two.y);
    return newPoint;
}

//Умножение на скаляр справа
Point operator*(const Point& one, double coefficient) {
    Point newPoint(one.x * coefficient, one.y*coefficient);
    return newPoint;
}

//Умножение на скаляр слева
Point operator*(double coefficient, Point& one) {
    return one*coefficient;
}

double GetVectorsAngle(const Point& one, const Point& two) {
    return acos((one.x*two.x + one.y*two.y)/(one.getLength()*two.getLength()));
}

//==================Line class methods========================//

//Класс выброса исклучения при невозможности создания прямой
class NoLineException : public std::exception {
    const char* what() const noexcept override;
};

const char *NoLineException::what() const noexcept {
    return "Нельзя провести прямую через одну точку";
}

Line::Line(const Point &first, const Point &second) {
    this->one = first;
    this->two = second;
    if (first == second) {
        throw NoLineException();
    }
}

Line::Line(double angle, double height) {
    this->one = Point(0, height);
    this->two = Point(-(height/angle), 0);
}

Line::Line(const Point &point, double angle) {
    this->one = point;
    this->two = Point(point.x + 1, point.y + angle);
}

bool Line::operator==(const Line& other) const {
    Point pointOne = other.one;
    Point pointTwo = other.two;

    //Условие 1: Первая точка второй прямой находится на первой прямой
    bool condOne = (this->one.x - this->two.x)*(pointOne.y - this->two.y)
                   - (pointOne.x - this->two.x)*(this->one.y - this->two.y) < eps;

    //Условие 2: Вторая точка второй прямой находится на первой прямой
    bool condTwo = (this->one.x - this->two.x)*(pointTwo.y - this->two.y)
                   - (pointTwo.x - this->two.x)*(this->one.y - this->two.y) < eps;

    //Оба условия должны выполняться
    return condOne && condTwo;
}

bool Line::operator!=(const Line& other) const {
    return !(*this == other);
}

Point Line::getNormalVector() const {
    return Point(this->one.y - this->two.y, this->two.x - this->one.x);
}

Point GetReflectedPoint(const Point& point, const Line& line) {
    Point normalVector = line.getNormalVector();
    double normalVectorLength = normalVector.getLength();
    normalVector.x /= normalVectorLength;
    normalVector.y /= normalVectorLength;

    double ortoVectorLength = 2 * (normalVector.x * point.x + normalVector.y * point.y);

    Point reflected(point.x - ortoVectorLength * normalVector.x, point.y - ortoVectorLength * normalVector.y);
    return reflected;
}

Point GetScaledPoint(const Point& point, const Point& center, double coefficient) {
    Point vector = point - center;
    vector = vector * coefficient;
    Point newPoint = vector + center;
    return newPoint;
}

double GetDeterminant(double x11, double x12, double x21, double x22) {
    return x11 * x22 - x12 * x21;
}

double GetAngle(const Point& firstVector, const Point& secondVector) {
    double angle = acos((firstVector.x*secondVector.x + firstVector.y*secondVector.y)/
                        (firstVector.getLength()* secondVector.getLength()));
    return angle;
}

//==================Polygon class methods========================//

double Polygon::perimeter() const {
    double sum = 0.0;
    for (size_t i = 0; i < this->vertices.size(); ++i) {
        sum += (this->vertices[(i+1) % this->vertices.size()] - this->vertices[i]).getLength();
    }
    return sum;
}

//Используется формула площади Гаусса
double Polygon::area() const {
    size_t pointer = 0;
    double area = 0;
    for (size_t i = 0; i < this->vertices.size(); ++i) {
        area += GetDeterminant(this->vertices[pointer].x,
                               this->vertices[pointer].y,
                               this->vertices[(pointer + 1) % this->vertices.size()].x,
                               this->vertices[(pointer + 1) % this->vertices.size()].y
        );
        pointer++;
    }
    return std::abs(area)*0.5;
}

bool Polygon::verticesVectorsEqual(std::vector<Point>& thisCopy, std::vector<Point>& otherCopy) {
    for (size_t j = 0; j < thisCopy.size()+1; ++j) {
        bool equal = true;
        for (size_t i = 0; i < thisCopy.size(); ++i) {
            if (thisCopy[i] != otherCopy[i]) {
                equal = false;
                break;
            }
        }
        std::rotate(thisCopy.begin(), thisCopy.begin() + 1, thisCopy.end());
        if (equal) {
            return true;
        }
    }
    return false;
}

bool Polygon::operator==(const Shape &otherShape) const {
    //Проверяем, что является многоугольником или ребенком класса Polygon
    const auto otherPointer = dynamic_cast<const Polygon*>(&otherShape);
    if (!otherPointer) {
        return false;
    }
    Polygon other = *otherPointer;

    if (this->vertices.size() != other.vertices.size()) {
        return false;
    }

    //Путем поворота проверяем, что вершины совпадают
    std::vector<Point> thisCopy = this->vertices;
    bool equal = this->verticesVectorsEqual(thisCopy, other.vertices);
    if(equal) {
        return true;
    }

    //Обходим в другую сторону
    std::reverse(thisCopy.begin(), thisCopy.end());
    equal = this->verticesVectorsEqual(thisCopy, other.vertices);

    return equal;
}


bool Polygon::isSimilarTo(const Shape &otherShape) const {
    const auto otherPointer = dynamic_cast<const Polygon*>(&otherShape);
    if (!otherPointer) {
        return false;
    }
    Polygon other = *otherPointer;

    if (this->vertices.size() != other.vertices.size()) {
        return false;
    }

    std::vector<Point> edgesThis;
    std::vector<Point> edgesOther;

    for (size_t i = 0; i < this->vertices.size(); ++i) {
        Point firstEdge = this->vertices[i]-this->vertices[(i+1) % this->vertices.size()];
        Point firstEdgeOther = other.vertices[i]-other.vertices[(i+1) % other.vertices.size()];
        edgesThis.push_back(firstEdge);
        edgesOther.push_back(firstEdgeOther);
    }

    //Проходимся по каждому ребру данного многоугольника, и делим его длину на длину нулевого ребра
    //другого многоугольника. Потом на данный коэфициент меняем все ребра данного многоугольника
    //и выполняем проверку равенства. Если многоугольники подобны, то рано или поздно мы найдем ребро,
    //дающее походящий коэфициент подобия
    for (size_t i = 0; i < this->vertices.size(); ++i) {
        double koef = edgesThis[i].getLength() / edgesOther[0].getLength();

        std::vector<Point> otherCopy = edgesOther;
        for (size_t j = 0; j < edgesOther.size(); ++j) {
            otherCopy[j] = GetScaledPoint(otherCopy[j], Point(0,0), koef);
        }
        std::vector<Point> scaledPolygonVector{Point(0, 0)};
        for (size_t j = 0; j < edgesOther.size()-1; ++j) {
            Point nextPoint = scaledPolygonVector.back() + otherCopy[j];
            scaledPolygonVector.push_back(nextPoint);
        }
        Polygon scaledPolygon(scaledPolygonVector);
        if (this->isCongruentTo(scaledPolygon)) {
            return true;
        }
    }
    return false;
}

bool Polygon::scalarAndEdgesEqual(std::vector<double> &scalarThis, std::vector<double> &scalarOther,
                                  std::vector<double> &edgesThis, std::vector<double> &edgesOther)
{
    for (size_t i = 0; i < scalarOther.size(); ++i) {
        bool equal = true;
        for (size_t j = 0; j < scalarOther.size(); ++j) {
            if (! ((std::abs(scalarThis[j] - scalarOther[j]) < eps) &&
                   (std::abs(edgesThis[j]-edgesOther[j]) < eps)))
            {
                equal = false;
                break;
            }
        }
        if (equal) {
            return true;
        }
        std::rotate(scalarOther.begin(), scalarOther.begin() + 1, scalarOther.end());
        std::rotate(edgesOther.begin(), edgesOther.begin() + 1, edgesOther.end());
    }
    return false;
}

bool Polygon::isCongruentTo(const Shape &otherShape) const {
    const auto otherPointer = dynamic_cast<const Polygon*>(&otherShape);
    if (!otherPointer) {
        return false;
    }
    Polygon other = *otherPointer;

    if (this->vertices.size() != other.vertices.size()) {
        return false;
    }

    //Получаем списки скалярных произведений
    std::vector<double> scalarThis;
    std::vector<double> scalarOther;
    std::vector<double> edgesThis;
    std::vector<double> edgesOther;

    for (size_t i = 0; i < this->vertices.size(); ++i) {
        size_t firstPointer = i;
        size_t secondPointer = (i+1) % this->vertices.size();
        size_t thirdPointer = (i+2) % this->vertices.size();

        Point firstEdge = this->vertices[firstPointer] - this->vertices[secondPointer];
        Point secondEdge = this->vertices[secondPointer] - this->vertices[thirdPointer];

        Point firstEdgeOther = other.vertices[firstPointer]-other.vertices[secondPointer];
        Point secondEdgeOther = other.vertices[secondPointer] - other.vertices[thirdPointer];

        scalarThis.push_back((firstEdge).getScalar(secondEdge));
        scalarOther.push_back((firstEdgeOther).getScalar(secondEdgeOther));
        edgesThis.push_back(firstEdge.getLength());
        edgesOther.push_back(firstEdgeOther.getLength());
    }


    bool equal = this->scalarAndEdgesEqual(scalarThis, scalarOther, edgesThis, edgesOther);
    if (equal) {
        return true;
    }

    std::reverse(edgesOther.begin(), edgesOther.end());
    std::reverse(scalarOther.begin(), scalarOther.end());
    std::rotate(scalarOther.begin(), scalarOther.begin()+1, scalarOther.end());

    equal = this->scalarAndEdgesEqual(scalarThis, scalarOther, edgesThis, edgesOther);

    return equal;

}

//Испускаем луч, проверяем, сколько раз пересечет грани многоугольника
bool Polygon::containsPoint(const Point& point) const {
    bool inside = false;
    for (size_t i = 0, j = this->vertices.size()-1; i < this->vertices.size(); j = i++) {
        if ((this->vertices[i].y > point.y) != (this->vertices[j].y > point.y) &&
            (point.x < (this->vertices[j].x-this->vertices[i].x) * (point.y-this->vertices[i].y) /
                       (this->vertices[j].y-this->vertices[i].y) +this->vertices[i].x)
                )
        {
            inside = !inside;
        }
    }
    return inside;
}

void Polygon::rotate(const Point &center, double angle) {
    for (Point& vertex : this->vertices) {
        Point vector = vertex - center;
        vector.rotate((angle/360)*2*pi);
        vertex = vector + center;
    }
}

void Polygon::reflex(const Line &axis) {
    for (Point& vertex : this->vertices) {
        vertex = GetReflectedPoint(vertex, axis);
    }
}

void Polygon::scale(const Point& center, double coefficient) {
    for (Point& vertex : this->vertices) {
        vertex = GetScaledPoint(vertex, center, coefficient);
    }
}

Polygon::Polygon(const Polygon& other) {
    this->vertices = other.vertices;
}

Polygon::Polygon(const std::vector<Point>& points) {
    this->vertices = points;
}

bool Polygon::operator!=(const Shape& otherShape) const {
    const auto otherPointer = dynamic_cast<const Polygon*>(&otherShape);
    if (!otherPointer) {
        return false;
    }
    Polygon other = *otherPointer;
    return !(*this == other);
}

//Узнаем, совершен ли был поворот вектора при обходе многоугольника по часовой стрелке и против часовой
//Если оба поворота были совершены, то значит многоугольник невыпуклый
bool Polygon::isConvex() const {
    bool positive = false;
    bool negative = false;
    for (size_t i = 0; i < this->vertices.size(); ++i) {
        size_t firstIndex = i;
        size_t secondIndex = (i + 1) % this->vertices.size();
        size_t thirdIndex = (i + 2) % this->vertices.size();

        double deltaxOne = this->vertices[secondIndex].x - this->vertices[firstIndex].x;
        double deltaxTwo = this->vertices[thirdIndex].x - this->vertices[secondIndex].x;
        double deltayOne = this->vertices[secondIndex].y - this->vertices[firstIndex].y;
        double deltayTwo = this->vertices[thirdIndex].y - this->vertices[secondIndex].y;
        double zComponentProduct = deltaxOne*deltayTwo - deltayOne*deltaxTwo;
        if (zComponentProduct > 0) {
            positive = true;
        }
        if (zComponentProduct < 0) {
            negative = true;
        }
    }
    return positive ^ negative;
}

std::vector<Point> Polygon::getVertices() const {
    return this->vertices;
}

Polygon::Polygon(const Point& first) {
    this->vertices.insert(this->vertices.begin(), first);
}

template<typename... T>
Polygon::Polygon(const Point &first, const T &... other) : Polygon(other...) {
    this->vertices.insert(this->vertices.begin(), first);
}

Polygon Polygon::GetMinkowskiSum(const Polygon &other) const {
    size_t firstPos  = 0;
    size_t secondPos = 0;
    std::vector<Point> MinkowskiVertices = {};

    while(firstPos < this->vertices.size() || secondPos < other.vertices.size()) {
        if(firstPos + secondPos > this->vertices.size() + other.vertices.size()) {
            break;
        }
        MinkowskiVertices.push_back(this->vertices[firstPos % this->vertices.size()] +
                                    other.vertices[secondPos % other.vertices.size()]);

        size_t nextFirstPos  = (firstPos + 1) % this->vertices.size();
        size_t nextSecondPos = (secondPos + 1) % other.vertices.size();

        auto fFirst = this->vertices[firstPos % this->vertices.size()];
        auto fSecond = this->vertices[nextFirstPos];

        auto sFirst = other.vertices[secondPos % other.vertices.size()];
        auto sSecond = other.vertices[nextSecondPos];

        auto fSub = fSecond - fFirst;
        auto sSub = sSecond - sFirst;

        std::complex<double> fcSub(fSub.x, fSub.y);
        std::complex<double> scSub(sSub.x, sSub.y);

        double fAngle = std::arg(fcSub);
        double sAngle = std::arg(scSub);

        if(std::abs(fAngle - sAngle) < eps) {
            ++firstPos;
            ++secondPos;
        }
        else if(fAngle < sAngle) {
            ++firstPos;
        } else if(fAngle > sAngle) {
            ++secondPos;
        }
    }

    return Polygon(MinkowskiVertices);
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    size_t firstPolygonSize = 0;
    std::cin >> firstPolygonSize;

    std::vector<Point> verticesFirst = {};

    for(size_t i = 0; i < firstPolygonSize; ++i) {
        double x = 0;
        double y = 0;
        std::cin >> x >> y;
        verticesFirst.push_back(Point{x,y});
    }

    size_t min_point = 0;
    for(size_t i = 0; i < verticesFirst.size(); ++i) {
        if(verticesFirst[i].x < verticesFirst[min_point].x ||
           (std::abs(verticesFirst[i].x - verticesFirst[min_point].x) &&
            verticesFirst[i].y < verticesFirst[min_point].y)) {
            min_point = i;
        }
    }

    std::rotate(verticesFirst.begin(), verticesFirst.begin() + min_point, verticesFirst.end());

    auto first = Polygon(verticesFirst);

    size_t secondPolygonSize = 0;
    std::cin >> secondPolygonSize;

    std::vector<Point> verticesSecond = {};

    for(size_t i = 0; i < secondPolygonSize; ++i) {
        double x = 0;
        double y = 0;
        std::cin >> x >> y;
        verticesSecond.push_back(Point{x,y});
    }

    min_point = 0;
    for(size_t i = 0; i < verticesSecond.size(); ++i) {
        if(verticesSecond[i].x < verticesSecond[min_point].x ||
           (std::abs(verticesSecond[i].x - verticesSecond[min_point].x) &&
            verticesSecond[i].y < verticesSecond[min_point].y)) {
            min_point = i;
        }
    }

    std::rotate(verticesSecond.begin(), verticesSecond.begin() + min_point, verticesSecond.end());

    auto second = Polygon(verticesSecond);

    auto result = first.GetMinkowskiSum(second);

    std::cout << std::fixed << std::setprecision(6) << (Polygon(result).area() - first.area() - second.area()) / 2 << std::endl;

    return 0;
}