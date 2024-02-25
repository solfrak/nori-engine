#include <nori/common.h>
#include <set>
NORI_NAMESPACE_BEGIN

class Node
{
    public:
        Node *build(TBoundingBox<Point3f> box, std::vector<int> *listTriangleIndice, Mesh *mesh, int depth);
        Node();
        TBoundingBox<Point3f> *nodeBox = nullptr;
        void rayIntersection(Ray3f ray, std::vector<int> *hitTriangle);
        std::vector<int> *triangles = nullptr;
        std::vector<Node *> child;
};

NORI_NAMESPACE_END