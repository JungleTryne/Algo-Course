#define DEBUG

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <list>
#include <memory>
#include <vector>
#include <utility>
#include <queue>
#include <set>

inline static constexpr long double eps() {
    return 1e-9;
};

inline static constexpr double subeps() {
    return 1e-5;
}

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

    double GetLength() const;
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
double Point<PointType>::GetLength() const {
    return sqrt(static_cast<double>(x) * static_cast<double>(x) + static_cast<double>(y) * static_cast<double>(y));
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
    //Getting circle's center, simplifying variables' names
    double x1 = first.x;
    double x2 = second.x;
    double x3 = third.x;

    double y1 = first.y;
    double y2 = second.y;
    double y3 = third.y;

    /*
     * This formula can be gotten if one expand matrix determinant formula
     * Actually circleCenterX = -B/2A, on the numerator it is -B and 2A for denominator
     */
    double circleCenterX =
            ( (x1 * x1 + y1 * y1) * (y2 - y3) + (x2 * x2 + y2 * y2) * (y3 - y1) + (x3 * x3 + y3 * y3) * (y1 - y2) ) /
            //---------------------------------------------------------------------------------------------------
                            ( 2 * (x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2) );

    /*
     * This formula can be gotten if one expand matrix determinant formula
     * Actually circleCenterY = -C/2A, on the numerator it is -C and 2A for denominator
     */
    double circleCenterY =
            ( (x1 * x1 + y1 * y1) * (x3 - x2) + (x2 * x2 + y2 * y2) * (x1 - x3) + (x3 * x3 + y3 * y3) * (x2 - x1) ) /
            //---------------------------------------------------------------------------------------------------
                            ( 2 * (x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2) );

    return Point(circleCenterX, circleCenterY);
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
    class BeachLineNode;

    enum class EventType {SITE, CIRCLE};
    struct Event {
        double sweep_line_coord;
        EventType type;
        bool valid;

        std::optional<Point<double>> point;
        BeachLineNode* archToDelete = nullptr;

        explicit Event(double sweep_line, EventType type, bool valid);
        bool operator<(const Event& other) const;
    };

    struct SharedPtrComparator {
        bool operator()(const std::shared_ptr<Event>& first, const std::shared_ptr<Event>& second) const {
            return first->sweep_line_coord < second->sweep_line_coord;
        }
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

        std::shared_ptr<Event> circleEvent_ = nullptr;
        DCEL::HalfEdge* halfEdge_ = nullptr;

        BeachLineNode* next_ = nullptr;
        BeachLineNode* prev_ = nullptr;

        bool IsLeaf() const { return archs_.first == archs_.second; }

        bool IsRoot() const { return !parent_; }
        int32_t GetHeight() const { return height_; }

        int32_t GetBalance() const {return left_->height_ - right_->height_; }

        double GetValue() const;
        friend BeachLine;


        void line_print() const;

    public:
        BeachLineNode* GetPrev() const { return this->prev_; }
        BeachLineNode* GetNext() const { return this->next_; }
        Point<double> GetLeaf() const { return archs_.first; }
        BeachLineNode(BeachLineNode* parent, double* sweepLineCoord, std::pair<Point<double>, Point<double>> archs);

        void SetCircleEvent(std::shared_ptr<Event> circleEvent) {
            circleEvent_ = circleEvent;
        }
    };

    class BeachLine {
    private:
        double* sweepLineCoord_ = nullptr;
        std::unique_ptr<BeachLineNode> root = nullptr;

        BeachLineNode* FindIntersectingArch(const Point<double>& newFocus) const;
        BeachLineNode* InsertToTheLeft(BeachLineNode* node, const Point<double>& newFocus);
        BeachLineNode* InsertToTheRight(BeachLineNode* node, const Point<double>& newFocus);

        BeachLineNode* GetNextNodeNeighbor(BeachLineNode* node) const;
        BeachLineNode* GetPrevNodeNeighbor(BeachLineNode* node) const;

        void print_line() const;

    public:
        BeachLineNode* InsertNewArch(const Point<double>& newFocus);
        std::pair<BeachLineNode*, BeachLineNode*> DeleteArch(BeachLineNode* node);
        BeachLine(double* sweepLine) : sweepLineCoord_(sweepLine) {}

        void printTree() const;

        friend VoronoiDiagramBuilder;
    };

    std::priority_queue<std::shared_ptr<Event>, std::vector<std::shared_ptr<Event>>, SharedPtrComparator> queue_;
    std::vector<Point<double>> points_;
    VoronoiDiagram diagram_;
    double sweepLineCoord_ = infty();

    void HandleSiteEvent(const Event& event);
    void HandleCircleEvent(const Event& event);

    std::shared_ptr<Event> CheckCircleEvent(BeachLineNode* first, BeachLineNode* second, BeachLineNode* third) const;

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
            return intersections[1].x;
        }
        return intersections[0].x;

    }
}

void VoronoiDiagramBuilder::BeachLineNode::line_print() const {
    std::cout << "(" << this->GetLeaf().x << ", " << this->GetLeaf().y << ") ";
    if(this->next_) {
        this->next_->line_print();
    }
}


VoronoiDiagramBuilder::BeachLineNode*
VoronoiDiagramBuilder::BeachLine::FindIntersectingArch(const Point<double> &newFocus) const {
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

VoronoiDiagramBuilder::BeachLineNode* VoronoiDiagramBuilder::BeachLine::InsertNewArch(const Point<double> &newFocus) {
    if(!root) {
        root = std::make_unique<BeachLineNode> (
                nullptr,
                sweepLineCoord_,
                std::make_pair(newFocus, newFocus)
        );
        root->height_ = 0;
        return root.get();
    }
    BeachLineNode* ArchToBreak = FindIntersectingArch(newFocus);

    if(ArchToBreak->circleEvent_) {
        ArchToBreak->circleEvent_->valid = false;
    }

    double value = ArchToBreak->GetValue();
    if(value <= newFocus.x - eps()) {
        return InsertToTheRight(ArchToBreak, newFocus);
    }
    return InsertToTheLeft(ArchToBreak, newFocus);
 }

VoronoiDiagramBuilder::BeachLineNode*
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

    node->right_->prev_ = leftNode->right_.get();
    leftNode->right_->prev_ = leftNode->left_.get();
    leftNode->left_->prev_ = node->prev_;

    leftNode->left_->next_ = leftNode->right_.get();
    leftNode->right_->next_ = node->right_.get();
    node->right_->next_ = node->next_;

    if(node->next_) {
        node->next_->prev_ = node->right_.get();
    }
    if(node->prev_) {
        node->prev_->next_ = leftNode->left_.get();
    }

    node->next_ = nullptr;
    node->prev_ = nullptr;

    return leftNode->right_.get();
}

VoronoiDiagramBuilder::BeachLineNode*
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
            std::make_pair(rightNode->archs_.second, rightNode->archs_.second)
    );

    rightNode->left_ = std::make_unique<BeachLineNode>(
            rightNode,
            sweepLineCoord_,
            std::make_pair(rightNode->archs_.first, rightNode->archs_.first)
    );

    node->left_->next_ = rightNode->left_.get();
    rightNode->left_->next_ = rightNode->right_.get();
    rightNode->right_->next_ = node->next_;

    rightNode->right_->prev_ = rightNode->left_.get();
    rightNode->left_->prev_ = node->left_.get();
    node->left_->prev_ = node->prev_;

    if(node->next_) {
        node->next_->prev_ = rightNode->right_.get();
    }
    if(node->prev_) {
        node->prev_->next_ = node->left_.get();
    }

    node->next_ = nullptr;
    node->prev_ = nullptr;

    return rightNode->left_.get();
}

std::shared_ptr<VoronoiDiagramBuilder::Event> VoronoiDiagramBuilder::CheckCircleEvent(
        VoronoiDiagramBuilder::BeachLineNode *first,
        VoronoiDiagramBuilder::BeachLineNode *second,
        VoronoiDiagramBuilder::BeachLineNode *third) const
{
    if(!(first && second && third)) {
        return nullptr;
    }

    auto first_point = first->GetLeaf();
    auto second_point = second->GetLeaf();
    auto third_point = third->GetLeaf();

    //checking the right turn
    if((second_point.x - first_point.x) * (third_point.y - first_point.y) -
            (third_point.x - first_point.x) * (second_point.y - first_point.y) > 0) {
        return nullptr;
    }

    auto circleCenter = GetCenterOfCircle(first_point, second_point, third_point);
    if(circleCenter.y - (circleCenter - first_point).GetLength() > sweepLineCoord_) {
        return nullptr;
    }

    double circleRadius = (circleCenter - first_point).GetLength();
    std::shared_ptr<Event> circleEvent = std::make_shared<Event>(circleCenter.y - circleRadius, EventType::CIRCLE, true);
    circleEvent->point = circleCenter;
    circleEvent->archToDelete = second;

    second->SetCircleEvent(circleEvent);
    return circleEvent;
}

VoronoiDiagramBuilder::BeachLineNode *
VoronoiDiagramBuilder::BeachLine::GetNextNodeNeighbor(VoronoiDiagramBuilder::BeachLineNode *node) const {
    /*
     * Using visitor pattern to get next neighbour for an arch
     * node must be leaf initially
     */

    auto currentNode = node;
    if(!node->parent_) {
        return node;
    }

    currentNode = currentNode->parent_;
    while(currentNode->parent_->right_.get() == currentNode) {
        currentNode = currentNode->parent_;
    }

    currentNode = currentNode->right_.get();
    while(!currentNode->IsLeaf()) {
        currentNode = currentNode->left_.get();
    }

    return currentNode;
}

VoronoiDiagramBuilder::BeachLineNode *
VoronoiDiagramBuilder::BeachLine::GetPrevNodeNeighbor(VoronoiDiagramBuilder::BeachLineNode *node) const {
    /*
 * Using visitor pattern to get next neighbour for an arch
 * node must be leaf initially
 */

    auto currentNode = node;
    if(!node->parent_) {
        return node;
    }

    currentNode = currentNode->parent_;
    while(currentNode->parent_->left_.get() == currentNode) {
        currentNode = currentNode->parent_;
    }

    currentNode = currentNode->left_.get();
    while(!currentNode->IsLeaf()) {
        currentNode = currentNode->right_.get();
    }

    return currentNode;
}

auto VoronoiDiagramBuilder::BeachLine::DeleteArch(VoronoiDiagramBuilder::BeachLineNode *node)
    -> std::pair<BeachLineNode*, BeachLineNode*>
{
    BeachLineNode* parent = node->parent_;
    assert(parent->parent_);
    auto superParent = parent->parent_;

    auto focusToDelete = node->GetLeaf();
    auto deletedPrev = node->prev_->GetLeaf();
    auto deletedNext = node->next_->GetLeaf();

    if(parent->right_.get() == node) {
        if(superParent->right_.get() == parent) {
            auto next = superParent->right_->right_->next_;

            superParent->right_->left_->parent_ = superParent;

            std::unique_ptr<BeachLineNode> tmp = nullptr;
            tmp.swap(superParent->right_->left_);
            tmp.swap(superParent->right_);

            auto mostRight = superParent->right_.get();
            while(mostRight->right_) {
                mostRight = mostRight->right_.get();
            }

            mostRight->next_ = next;
            if(next) {
                next->prev_ = mostRight;
            }

            auto currentNode = superParent;
            while(currentNode && !(currentNode->archs_.first == focusToDelete && currentNode->archs_.second == deletedNext)) {
                currentNode = currentNode->parent_;
            }
            if(currentNode) {
                currentNode->archs_ = std::make_pair(deletedPrev, deletedNext);
            }

            return std::make_pair(superParent->right_.get(), superParent->right_->next_);
        } else {
            auto next = superParent->left_->right_->next_;

            superParent->left_->left_->parent_ = superParent;

            std::unique_ptr<BeachLineNode> tmp = nullptr;
            tmp.swap(superParent->left_->left_);
            tmp.swap(superParent->left_);

            auto mostRight = superParent->left_.get();
            while(mostRight->right_) {
                mostRight = mostRight->right_.get();
            }

            mostRight->next_ = next;
            if(next) {
                next->prev_ = mostRight;
            }

            auto currentNode = superParent;
            while(currentNode && !(currentNode->archs_.first == deletedPrev && currentNode->archs_.second == focusToDelete)) {
                currentNode = currentNode->parent_;
            }
            if(currentNode) {
                currentNode->archs_ = std::make_pair(deletedPrev, deletedNext);
            }

            return std::make_pair(superParent->left_.get(), superParent->left_->next_);
        }

    } else {
        if(superParent->right_.get() == parent) {
            auto prev = superParent->right_->left_->prev_;

            superParent->right_->right_->parent_ = superParent;

            std::unique_ptr<BeachLineNode> tmp = nullptr;
            tmp.swap(superParent->right_->right_);
            tmp.swap(superParent->right_);

            auto mostLeft = superParent->right_.get();
            while(mostLeft->left_) {
                mostLeft = mostLeft->left_.get();
            }

            mostLeft->prev_ = prev;
            if(prev) {
                prev->next_ =mostLeft;
            }

            auto currentNode = superParent;
            while(currentNode && !(currentNode->archs_.first == focusToDelete && currentNode->archs_.second == deletedNext)) {
                currentNode = currentNode->parent_;
            }
            if(currentNode) {
                currentNode->archs_ = std::make_pair(deletedPrev, deletedNext);
            }

            return std::make_pair(superParent->right_->prev_, superParent->right_.get());
        } else {
            auto prev = superParent->left_->left_->prev_;

            superParent->left_->right_->parent_ = superParent;

            std::unique_ptr<BeachLineNode> tmp = nullptr;
            tmp.swap(superParent->left_->right_);
            tmp.swap(superParent->left_);

            auto mostLeft = superParent->left_.get();
            while(mostLeft->left_) {
                mostLeft = mostLeft->left_.get();
            }

            mostLeft->prev_ = prev;
            if(prev) {
                prev->next_ =mostLeft;
            }

            auto currentNode = superParent;
            while(currentNode && !(currentNode->archs_.first == deletedPrev && currentNode->archs_.second == focusToDelete)) {
                currentNode = currentNode->parent_;
            }
            if(currentNode) {
                currentNode->archs_ = std::make_pair(deletedPrev, deletedNext);
            }

            return std::make_pair(superParent->left_->prev_, superParent->left_.get());
        }
    }
}

void VoronoiDiagramBuilder::BeachLine::print_line() const {
    if(root) {
        auto currentNode = root.get();
        while(currentNode->left_) {
            currentNode = currentNode->left_.get();
        }
        currentNode->line_print();
        std::cout << '\n';
    }
}

void VoronoiDiagramBuilder::BeachLine::printTree() const {
    std::queue<std::pair<BeachLineNode*, size_t>> queue;
    if(!this->root) {
        return;
    }

    size_t currentLevel = 0;
    queue.emplace(this->root.get(), 0);
    while(!queue.empty()) {
        auto currentNode = queue.front().first;
        size_t poppedLevel = queue.front().second;
        queue.pop();

        if(currentLevel != poppedLevel) {
            currentLevel = poppedLevel;
            std::cout << std::endl;
        }

        std::cout << '[' << "(" << currentNode->archs_.first.x << ", " <<  currentNode->archs_.first.y << "), "
            << "(" << currentNode->archs_.second.x << ", " <<  currentNode->archs_.second.y<< ")] ";
        if(currentNode->left_.get()) {
            queue.emplace(currentNode->left_.get(), poppedLevel + 1);
            queue.emplace(currentNode->right_.get(), poppedLevel + 1);
        }
    }
    std::cout << std::endl;
    std::cout << std::endl;
}

void VoronoiDiagramBuilder::Build() {
    for(const auto& point : points_) {
        auto new_event = std::make_shared<Event>(point.y, EventType::SITE, true);
        new_event->point = point;
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

#ifdef DEBUG
        line.print_line();
        //std::cout << "====================" << std::endl;
        //line.printTree();
#endif
    }
}

void VoronoiDiagramBuilder::HandleSiteEvent(const VoronoiDiagramBuilder::Event &event) {
    assert(event.point.has_value());

    BeachLineNode* newLeaf = this->line.InsertNewArch(*event.point);
    BeachLineNode* nextNewLeaf = newLeaf->GetNext();
    BeachLineNode* nextNextNewLeaf = nullptr;
    if(nextNewLeaf) {
        nextNextNewLeaf = nextNewLeaf->GetNext();
    }

    std::shared_ptr<Event> firstCircleEvent = nullptr;

    if(newLeaf && nextNewLeaf && nextNextNewLeaf) {
        firstCircleEvent = CheckCircleEvent(newLeaf, nextNewLeaf, nextNextNewLeaf);
    }

    BeachLineNode* prevNewLeaf = newLeaf->GetPrev();
    BeachLineNode* prevPrevNewLeaf = nullptr;
    if(prevNewLeaf) {
        prevPrevNewLeaf = prevNewLeaf->GetPrev();
    }

    std::shared_ptr<Event> secondCircleEvent = nullptr;

    if(newLeaf && prevNewLeaf && prevPrevNewLeaf) {
        secondCircleEvent = CheckCircleEvent(prevPrevNewLeaf, prevNewLeaf, newLeaf);
    }

    if(firstCircleEvent) {
        queue_.push(firstCircleEvent);
    }

    if(secondCircleEvent) {
        queue_.push(secondCircleEvent);
    }

    //TODO: insert lines to Voronoi diagram

}

void VoronoiDiagramBuilder::HandleCircleEvent(const VoronoiDiagramBuilder::Event &event) {
    auto circleCenter = event.point;

    //to be able to find arch in the tree

    BeachLineNode* nodeToDelete = event.archToDelete;
    assert(nodeToDelete);

    auto [prev, next] = line.DeleteArch(nodeToDelete);

    sweepLineCoord_ -= subeps();

    std::shared_ptr<Event> firstCircleEvent = nullptr;
    if(prev && prev->GetPrev() && prev->GetNext()) {
        firstCircleEvent = CheckCircleEvent(prev->GetPrev(), prev, prev->GetNext());
    }

    std::shared_ptr<Event> secondCircleEvent = nullptr;
    if(next && next->GetPrev() && next->GetNext()) {
        secondCircleEvent = CheckCircleEvent(next->GetPrev(), next, next->GetNext());
    }

    sweepLineCoord_ += subeps();

    if(firstCircleEvent) {
        queue_.push(firstCircleEvent);
    }

    if(secondCircleEvent) {
        queue_.push(secondCircleEvent);
    }

    //TODO: insert vertex to Voronoi diagram
}

bool VoronoiDiagramBuilder::Event::operator<(const VoronoiDiagramBuilder::Event &other) const {
    return sweep_line_coord < other.sweep_line_coord;
}

VoronoiDiagramBuilder::Event::Event(double sweep_line, VoronoiDiagramBuilder::EventType type, bool valid) :
    sweep_line_coord(sweep_line), type(type), valid(valid) {}

int main() {
    double x = 0.0;
    double y = 0.0;

    std::vector<Point<double>> points;

    while(std::cin >> x) {
        std::cin >> y;
        points.emplace_back(x, y);
    }

    VoronoiDiagramBuilder builder(std::move(points));
    builder.Build();
    return 0;
}