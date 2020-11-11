#include <algorithm>
#include <cassert>
#include <cmath>
#include <list>
#include <memory>
#include <vector>
#include <utility>
#include <queue>
#include <set>

inline static constexpr long double eps() {
    return 1e-9;
};

inline static constexpr long double infty() {
    return std::numeric_limits<double>::infinity();
}

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
        return std::abs(this->x - other.x) <= eps() && std::abs(this->y - other.y) <= eps();
    } else {
        return (this->x == other.x) && (this->y == other.y);
    }

}

template<typename PointType>
bool Point<PointType>::operator!=(const Point<PointType>& other) const {
    return !(*this==other);
}

template<typename PointType>
Point<PointType> operator-(const Point<PointType>& one, const Point<PointType>& two) {
    Point newPoint(one.x - two.x, one.y - two.y);
    return newPoint;
}

template<typename PointType>
Point<PointType> operator+(const Point<PointType>& one, const Point<PointType>& two) {
    Point newPoint(one.x + two.x, one.y + two.y);
    return newPoint;
}

template<typename PointType>
Point<PointType> operator*(const Point<PointType>& one, PointType coefficient) {
    Point newPoint(one.x * coefficient, one.y*coefficient);
    return newPoint;
}

template<typename PointType>
Point<PointType> operator*(PointType coefficient, Point<PointType>& one) {
    return one*coefficient;
}

template<typename PointType>
PointType GetDeterminant(PointType x11, PointType x12, PointType x21, PointType x22) {
    return x11 * x22 - x12 * x21;
}

std::vector<Point<double>> GetArchIntersection(const Point<double>& first, const Point<double>& second, double dir) {
    std::vector<Point<double>> result;

    if(std::abs(first.y - second.y) < eps()) {
        double y = (first.y + second.y) / 2;
        double delta = sqrt(dir * dir - dir * (first.y + second.y) + first.y * second.y);
        result.push_back(Point(first.x - delta, y));
        result.push_back(Point(first.x + delta, y));
        return result;
    }

    double length = pow(first.x - second.x, 2) + pow(first.y - second.y, 2);

    double delta = 2 * sqrt(pow(first.x - second.x, 2) * (dir - first.y) * (dir - second.y) * length);
    double fixed = -2 * dir * pow(first.x - second.x, 2) + (first.y + second.y) * length;

    double denominator = 2 * pow(first.y - second.y, 2);
    double delta_x = first.x - second.x;

    double y1 = (fixed - delta) / denominator;
    double y2 = (fixed + delta) / denominator;

    double x1 = (first.x * first.x - second.x * second.x +
            (2 * y1 - second.y - first.y) * (second.y - first.y)) / delta_x;

    double x2 = (first.x * first.x - second.x * second.x +
            (2 * y2 - second.y - first.y) * (second.y - first.y)) / delta_x;

    x1 /= 2;
    x2 /= 2;

    if (x1 > x2) {
        std::swap(x1, x2);
        std::swap(y1, y2);
    }

    result.push_back(Point<double>(x1, y1));
    result.push_back(Point<double>(x2, y2));
    return result;
}

Point<double> GetCenterOfCircle(const Point<double>& first, const Point<double>& second, const Point<double>& third) {

}


namespace DCEL {

    class Face;
    class HalfEdge;
    class Site;
    class Vertex;

    class Face {
        HalfEdge* edge;
        Site* site;
    };

    class HalfEdge {
        Face* face;

        HalfEdge* next;
        HalfEdge* prev;
        HalfEdge* twin;

        Vertex* start;
        Vertex* finish;
    };

    class Site {
        Face* face;
        Point<double> point;
        size_t index;
    };

    class Vertex {
        Point<double> point;
    };
}

struct VoronoiDiagram {
    std::vector<std::unique_ptr<DCEL::Face>> faces;
    std::vector<std::unique_ptr<DCEL::HalfEdge>> halfEdges;
    std::vector<std::unique_ptr<DCEL::Site>> sites;
    std::vector<std::unique_ptr<DCEL::Vertex>> vertices;
};

class VoronoiDiagramBuilder {
private:
    enum class EventType {SITE, CIRCLE};
    struct Event {
        double sweep_line_coord;
        EventType type;
        bool valid;

        std::optional<Point<double>> point;

        explicit Event(double sweep_line, EventType type, bool valid);
        bool operator<(const Event& other) const;
    };

    class BeachLine;

    class BeachLineNode {
    private:
        BeachLineNode* parent_;
        std::unique_ptr<BeachLineNode> left_ = nullptr;
        std::unique_ptr<BeachLineNode> right_ = nullptr;

        int32_t height_;
        double* sweepLineCoord_;

        std::pair<Point<double>, Point<double>> archs_;

        Event* circleEvent_ = nullptr;
        DCEL::HalfEdge* halfEdge_ = nullptr;

        BeachLineNode* next_ = nullptr;
        BeachLineNode* prev_ = nullptr;

        bool IsLeaf() const { return archs_.first == archs_.second; }

        bool IsRoot() const { return !parent_; }
        int32_t GetHeight() const { return height_; }

        int32_t GetBalance() const {return left_->height_ - right_->height_; }
        Point<double> GetLeaf() const { return archs_.first; }

        double GetValue() const;
        friend BeachLine;

    public:
        BeachLineNode(BeachLineNode* parent, double* sweepLineCoord, std::pair<Point<double>, Point<double>> archs);
    };

    class BeachLine {
    private:
        double* sweepLineCoord_ = nullptr;
        std::unique_ptr<BeachLineNode> root = nullptr;

        BeachLineNode* GetArchToInsertTo(const Point<double>& newFocus) const;
        BeachLineNode* InsertToTheLeft(BeachLineNode* node, const Point<double>& newFocus);
        BeachLineNode* InsertToTheRight(BeachLineNode* node, const Point<double>& newFocus);

        std::optional<Event> CheckCircleEvent(BeachLineNode* first, BeachLineNode* second, BeachLineNode* third) const;
    public:
        BeachLineNode* InsertNewArch(const Point<double>& newFocus);
        BeachLine(double* sweepLine) : sweepLineCoord_(sweepLine) {}
    };

    std::priority_queue<std::shared_ptr<Event>> queue_;
    std::vector<Point<double>> points_;
    VoronoiDiagram diagram_;
    double sweepLineCoord_ = -infty();

    void HandleSiteEvent(const Event& event);
    void HandleCircleEvent(const Event& event);

    BeachLine line{&sweepLineCoord_};

public:
    [[maybe_unused]] explicit VoronoiDiagramBuilder(const std::vector<Point<double>>& points) : points_(points) {}
    explicit VoronoiDiagramBuilder(std::vector<Point<double>>&& points) : points_(std::move(points)) {}
    void Build();
};

VoronoiDiagramBuilder::BeachLineNode::BeachLineNode(VoronoiDiagramBuilder::BeachLineNode *parent,
                                                    double *sweepLineCoord,
                                                    std::pair<Point<double>, Point<double>> archs) :
                                                    parent_(parent), sweepLineCoord_(sweepLineCoord),
                                                    archs_(archs) {}

double VoronoiDiagramBuilder::BeachLineNode::GetValue() const {
    if(this->IsLeaf()) {
        return this->archs_.first.x;
    } else {
        auto first_point = this->archs_.first;
        auto second_point = this->archs_.second;
        auto intersections = GetArchIntersection(first_point, second_point, *this->sweepLineCoord_);
        if(intersections.size() == 1) {
            return intersections[0].x;
        }

        if(first_point.y < second_point.y) {
            return intersections[0].x;
        }
        return intersections[1].x;
    }
}


VoronoiDiagramBuilder::BeachLineNode *
VoronoiDiagramBuilder::BeachLine::GetArchToInsertTo(const Point<double> &newFocus) const {
    BeachLineNode* currentNode = root.get();
    while(!currentNode->IsLeaf()) {
        double value = currentNode->GetValue();
        if(value <= newFocus.x - eps()) {
            currentNode = currentNode->right_.get();
        } else {
            currentNode = currentNode->left_.get();
        }
    }
    return currentNode;
}

VoronoiDiagramBuilder::BeachLineNode *VoronoiDiagramBuilder::BeachLine::InsertNewArch(const Point<double> &newFocus) {
    if(!root) {
        root = std::make_unique<BeachLineNode> (
                nullptr,
                sweepLineCoord_,
                std::make_pair(newFocus, newFocus)
        );
        root->height_ = 0;
        return root.get();
    }
    BeachLineNode* ArchToBreak = GetArchToInsertTo(newFocus);

    ArchToBreak->circleEvent_->valid = false;

    double value = ArchToBreak->GetValue();
    if(value <= newFocus.x - eps()) {
        return InsertToTheRight(ArchToBreak, newFocus);
    }
    return InsertToTheLeft(ArchToBreak, newFocus);
 }

VoronoiDiagramBuilder::BeachLineNode *
VoronoiDiagramBuilder::BeachLine::InsertToTheLeft(VoronoiDiagramBuilder::BeachLineNode *node,
                                                  const Point<double> &newFocus) {
    node->right_ = std::make_unique<BeachLineNode> (
            node,
            sweepLineCoord_,
            node->archs_
    );

    node->archs_ = std::make_pair(newFocus, node->archs_.second);

    node->left_ = std::make_unique<BeachLineNode>(
            node,
            sweepLineCoord_,
            std::make_pair(node->archs_.second, node->archs_.first)
    );

    auto leftNode = node->left_.get();

    leftNode->left_ = std::make_unique<BeachLineNode>(
            leftNode,
            sweepLineCoord_,
            std::make_pair(leftNode->archs_.first, leftNode->archs_.first)
    );

    leftNode->right_ = std::make_unique<BeachLineNode>(
            leftNode,
            sweepLineCoord_,
            std::make_pair(leftNode->archs_.second, leftNode->archs_.second)
    );

    return leftNode->right_.get();
}

VoronoiDiagramBuilder::BeachLineNode *
VoronoiDiagramBuilder::BeachLine::InsertToTheRight(VoronoiDiagramBuilder::BeachLineNode *node,
                                                   const Point<double> &newFocus) {
    node->left_ = std::make_unique<BeachLineNode> (
            node,
            sweepLineCoord_,
            node->archs_
    );

    node->archs_ = std::make_pair(node->archs_.second, newFocus);

    node->right_ = std::make_unique<BeachLineNode>(
            node,
            sweepLineCoord_,
            std::make_pair(node->archs_.second, node->archs_.first)
    );

    auto rightNode = node->right_.get();

    rightNode->right_ = std::make_unique<BeachLineNode>(
            rightNode,
            sweepLineCoord_,
            std::make_pair(rightNode->archs_.first, rightNode->archs_.first)
    );

    rightNode->left_ = std::make_unique<BeachLineNode>(
            rightNode,
            sweepLineCoord_,
            std::make_pair(rightNode->archs_.second, rightNode->archs_.second)
    );

    return rightNode->right_.get();
}

std::optional<VoronoiDiagramBuilder::Event> VoronoiDiagramBuilder::BeachLine::CheckCircleEvent(
        VoronoiDiagramBuilder::BeachLineNode *first,
        VoronoiDiagramBuilder::BeachLineNode *second,
        VoronoiDiagramBuilder::BeachLineNode *third) const
{
    if(!(first && second && third)) {
        return std::nullopt;
    }

    auto first_point = first->archs_.first;
    auto second_point = second->archs_.first;
    auto third_point = third->archs_.first;

    if(second_point.y > first_point.y && second_point.y > third_point.y) {
        return std::nullopt;
    }



}

void VoronoiDiagramBuilder::Build() {
    for(const auto& point : points_) {
        auto new_event = std::make_shared<Event>(point.y, EventType::SITE, true);
        queue_.push(std::move(new_event));
    }
    while(!queue_.empty()) {
        auto event = *queue_.top();
        queue_.pop();

        sweepLineCoord_ = event.sweep_line_coord;

        if(!event.valid) {
            continue;
        }
        if(event.type == EventType::SITE) {
            HandleSiteEvent(event);
        }
        else {
            HandleCircleEvent(event);
        }
    }
}

void VoronoiDiagramBuilder::HandleSiteEvent(const VoronoiDiagramBuilder::Event &event) {
    BeachLineNode* new_leaf = this->line.InsertNewArch(*event.point);

}

void VoronoiDiagramBuilder::HandleCircleEvent(const VoronoiDiagramBuilder::Event &event) {
    assert(false);
}

bool VoronoiDiagramBuilder::Event::operator<(const VoronoiDiagramBuilder::Event &other) const {
    return sweep_line_coord < other.sweep_line_coord;
}

VoronoiDiagramBuilder::Event::Event(double sweep_line, VoronoiDiagramBuilder::EventType type, bool valid) :
    sweep_line_coord(sweep_line), type(type), valid(valid) {}

int main() {

    return 0;
}