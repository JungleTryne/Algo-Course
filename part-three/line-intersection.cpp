#include <algorithm>
#include <iostream>
#include <set>
#include <type_traits>
#include <vector>

const double eps = 1e-9;

template<typename PointType>
struct Point {
    PointType x;
    PointType y;

    explicit Point(PointType x=0, PointType y=0) : x(x), y(y) {};
    Point(const Point& other) = default;
    ~Point() = default;

    bool operator==(const Point& other) const;
    bool operator!=(const Point& other) const;
};

template<typename PointType>
bool Point<PointType>::operator==(const Point &other) const {
    if constexpr (std::is_same<PointType, double>::value) {
        return std::abs(this->x - other.x) <= eps && std::abs(this->y - other.y) <= eps;
    } else {
        return this->x == other.x && this->y == other.y;
    }

}

template<typename PointType>
bool Point<PointType>::operator!=(const Point<PointType> &other) const {
    return !(*this == other);
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

template<typename PointType>
PointType GetArea(const Point<PointType>& first, const Point<PointType>& second, const Point<PointType>& third) {
    return (second.x - first.x) * (third.y - first.y) - (second.y - first.y) * (third.x - first.x);
}

template<typename PointType>
struct Segment {
        Point<PointType> first;
        Point<PointType> second;
        explicit Segment(Point<PointType> first, Point<PointType> second) : first(first), second(second) {}
        bool operator==(const Segment<PointType>& other) {
                return this->first == other.first && this->second == other.second;
        }
};

template<typename PointType>
bool operator<(const Segment<PointType>& first, const Segment<PointType>& second) {
    PointType left_x = std::max(std::min(first.first.x, first.second.x), std::min(second.first.x, second.second.x));

    auto relevant_y = [](const Segment<PointType>& segment, PointType x) -> double {
        if constexpr (std::is_same<PointType, double>::value) {
            if(std::abs(segment.first.x - segment.second.x) < eps) return segment.first.y;
        } else {
            if(segment.first.x == segment.second.x) return segment.first.y;
        }

        double tangent = static_cast<double>(segment.first.y - segment.second.y) /
                static_cast<double>(segment.first.x - segment.second.x);

        return static_cast<double>(segment.first.y) + tangent * static_cast<double>(x - segment.first.x);
    };

    return relevant_y(first, left_x) < relevant_y(second, left_x) - eps;
}

template<typename PointType>
[[nodiscard]] bool DoIntersect(const Segment<PointType>& first, const Segment<PointType>& second) {
    auto planar_intersection = [](PointType firstA, PointType firstB, PointType secondA, PointType secondB){
        return std::max(std::min(firstA, firstB), std::min(secondA, secondB)) <=
            std::min(std::max(firstA, firstB), std::max(secondA, secondB));
    };

    return GetArea(first.first, first.second, second.first) * GetArea(first.first, first.second, second.second) <= 0 &&
           GetArea(second.first, second.second, first.first) * GetArea(second.first, second.second, first.second) <= 0 &&
           planar_intersection(first.first.x, first.second.x, second.first.x, second.second.x) &&
           planar_intersection(first.first.y, first.second.y, second.first.y, second.second.y);
}

enum class EventType {
    NEW_SEG,
    END_SEG
};

template <typename PointType>
struct Event {
    PointType x_coord;
    EventType type;
    Segment<PointType> segment;

    Event(PointType x_coord, EventType type, Segment<PointType> segment) : x_coord(x_coord), type(type), segment(segment) {}
    bool operator<(const Event& other) const;
};

template<typename PointType>
bool Event<PointType>::operator<(const Event &other) const {
    if constexpr (std::is_same<PointType, double>::value) {
        return std::abs(this->x_coord - other.x_coord) < eps ? this->type < other.type : this->x_coord <= other.x_coord - eps;
    } else {
        return this->x_coord == other.x_coord ? this->type < other.type : this->x_coord < other.x_coord;
    }
}

template<typename PointType>
class IntersectionFinder {
private:
    const std::vector<Segment<PointType>> segments_;
    std::set<Segment<PointType>> priority_queue_;
public:
    explicit IntersectionFinder(const std::vector<Segment<PointType>>& segments) : segments_(segments) {
        priority_queue_.clear();
    }
    std::optional<std::pair<Segment<PointType>, Segment<PointType>>> IsThereIntersection();
};

template<typename PointType>
std::optional<std::pair<Segment<PointType>, Segment<PointType>>> IntersectionFinder<PointType>::IsThereIntersection() {
    std::vector<Event<PointType>> events;

    for(auto segment : segments_) {
        events.push_back(Event(std::min(segment.first.x, segment.second.x), EventType::NEW_SEG, segment));
        events.push_back(Event(std::max(segment.first.x, segment.second.x), EventType::END_SEG, segment));
    }
    std::stable_sort(events.begin(), events.end());

    for(auto event : events) {
        auto currentSegment = event.segment;
        if(event.type == EventType::NEW_SEG) {
                auto next_segment = priority_queue_.lower_bound(currentSegment);
                auto prev_segment = (next_segment == priority_queue_.begin()) ? priority_queue_.end() : std::prev(next_segment);

                if (next_segment != priority_queue_.end() && DoIntersect(*next_segment, currentSegment)) {
                    return std::make_pair(*next_segment, currentSegment);
                }

                else if(prev_segment != priority_queue_.end() && DoIntersect(*prev_segment, currentSegment))
                {
                    return std::make_pair(*prev_segment, currentSegment);
                }
                priority_queue_.insert(currentSegment);

        } else if (event.type == EventType::END_SEG) {
            auto segment_iterator = priority_queue_.find(event.segment);
            auto next_segment = std::next(segment_iterator);
            auto prev_segment = (next_segment == priority_queue_.begin()) ? priority_queue_.end() : std::prev(segment_iterator);
            if(next_segment != priority_queue_.end() && prev_segment != priority_queue_.end() && DoIntersect(*next_segment, *prev_segment)) {
                return std::make_pair(*prev_segment, *next_segment);
            }
            priority_queue_.erase(segment_iterator);
        }
    }

    return std::nullopt;
}

int main() {
    size_t segment_number = 0;
    std::cin >> segment_number;
    std::vector<Segment<int64_t>> segments = {};

    for(size_t i = 0; i < segment_number; ++i) {
        int64_t x1, y1, x2, y2;
        std::cin >> x1 >> y1 >> x2 >> y2;
        segments.emplace_back(Point(x1, y1), Point(x2, y2));
    }

    IntersectionFinder<int64_t> finder(segments);

    auto answer = finder.IsThereIntersection();
    if(answer.has_value()) {
        std::cout << "YES" << std::endl;
        const auto [first, second] = *answer;
        auto first_pos  = std::find(segments.begin(), segments.end(), first);

        auto second_pos = segments.end();
        for(auto pos = segments.begin(); pos != segments.end(); ++pos) {
            if(*pos == second && pos != first_pos) {
                second_pos = pos;
                break;
            }
        }

        if(std::distance(first_pos, second_pos) < 0) {
            std::swap(first_pos, second_pos);
        }
        std::cout << std::distance(segments.begin(), first_pos) + 1 << ' ' << std::distance(segments.begin(), second_pos) + 1 << std::endl;
    } else {
        std::cout << "NO" << std::endl;
    }

    return 0;
}
