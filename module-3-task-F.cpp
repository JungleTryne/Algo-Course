#include <algorithm>
#include <cmath>
#include <list>
#include <memory>
#include <vector>
#include <set>

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

struct Face;

struct SemiEdge;

struct Site {
    size_t index;
    Point<double> point;
    Face* face;
};

struct Face {
    Point<double>* vertex;
    SemiEdge* semiEdge;
};

struct SemiEdge {
    Point<double>* from;
    Point<double>* to;

    SemiEdge* previous;
    SemiEdge* next;

    SemiEdge* twin;
    Face* face;
};

enum class EventType {
    SITE,
    CIRCLE
};

class Event {
private:
    Point<double> point;
    EventType type;
public:
    Event(Point<double> point, EventType type) : point(point), type(type) {}

    Point<double> GetPoint() const;
    double GetY() const;
    EventType GetType() const;
};

double Event::GetY() const {
    return point.y;
}

EventType Event::GetType() const {
    return type;
}

Point<double> Event::GetPoint() const {
    return this->point;
}

bool operator<(const Event& first, const Event& second) {
    return first.GetY() < second.GetY();
}

struct Arch {
    Point<double> site;
    explicit Arch(Point<double> site) : site(site) {}
};

class BeachLine {
private:

public:
    void AddArch(Point<double> point);
    bool IsEmpty() const;
    std::optional<Arch> GetArchAbove(Point<double> point, double currentY);
    void SetRoot(Arch arch);
};

class VoronoiDiagram {
private:
    std::vector<Site> sites_;
    std::vector<Face> faces_;

    std::set<Event> priority_queue;

    //DCEL structure
    std::list<Point<double>> vertices_;
    std::list<SemiEdge> semiEdges_;

    BeachLine beach_line_;
    double currentY;

    void HandleSiteEvent(const Event& event);
    void HandleCircleEvent(const Event& event);
public:
    explicit VoronoiDiagram(const std::vector<Point<double>>& vertices);
};

void VoronoiDiagram::HandleSiteEvent(const Event& event) {
    if(beach_line_.IsEmpty()) {
        beach_line_.SetRoot(Arch(event.GetPoint()));
    }
    std::optional<Arch> ArchToDestroy = beach_line_.GetArchAbove(event.GetPoint(), currentY);

    if(ArchToDestroy.has_value()) {
    }

}

void VoronoiDiagram::HandleCircleEvent(const Event &event) {

}

VoronoiDiagram::VoronoiDiagram(const std::vector<Point<double>>& vertices) {
    for (auto vertex : vertices) {
        priority_queue.insert(Event(vertex, EventType::SITE));
    }
    while (!priority_queue.empty()) {
        auto event = *priority_queue.begin();
        priority_queue.erase(priority_queue.begin());
        if(event.GetType() == EventType::SITE) {
            HandleSiteEvent(event);
        } else {
            HandleCircleEvent(event);
        }
    }
}

int main() {

    return 0;
}