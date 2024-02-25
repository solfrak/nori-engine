#include <nori/common.h>
#include <nori/mesh.h>
#include <nori/node.h>
NORI_NAMESPACE_BEGIN

#define MAX_DEPTH_OCTREE 8

Node::Node()
{
    this->child = std::vector<Node *>(8);
}

Node* Node::build(TBoundingBox<Point3f> box, std::vector<int> *listTriangleIndice, Mesh *mesh, int depth)
{
     
    int triangleCount = listTriangleIndice->size();
    if(triangleCount == 0)
    {
        return nullptr;
    }

    //only a few triangles
    if(triangleCount <= 10 || depth >= 8)
    {
        Node *leaf = new Node();
        leaf->nodeBox = new TBoundingBox<Point3f>(box);
        leaf->triangles = new std::vector(*listTriangleIndice);
        // delete listTriangleIndice;
        return leaf;
    }


    std::vector<std::vector<int>> list = std::vector<std::vector<int>>(8);

    // std::vector<TBoundingBox<Point3f>> boxlist = std::vector<TBoundingBox<Point3f>>(8);
    TBoundingBox<Point3f> boxlist [8];

    //create 8 sub-box from the parent box
    Vector3f offsetX = Vector3f((box.max.x() - box.min.x()) / 2, 0, 0);
    Vector3f offsetY = Vector3f(0, (box.max.y() - box.min.y()) / 2, 0);
    Vector3f offsetZ = Vector3f(0, 0, (box.max.z() - box.min.z()) / 2);

    boxlist[0] = TBoundingBox<Point3f>(box.min, box.getCenter());
    boxlist[1] = BoundingBox3f(box.min + offsetX, box.getCenter() + offsetX);
    boxlist[2] = BoundingBox3f(box.min + offsetY, box.getCenter() + offsetY);
    boxlist[3] = BoundingBox3f(box.min + offsetX + offsetY, box.getCenter() + offsetX + offsetY);
    boxlist[4] = BoundingBox3f(box.min + offsetZ, box.getCenter() + offsetZ);
    boxlist[5] = BoundingBox3f(box.min + offsetZ + offsetX, box.getCenter() + offsetZ + offsetX);
    boxlist[6] = BoundingBox3f(box.min + offsetZ + offsetY, box.getCenter() + offsetZ + offsetY);
    boxlist[7] = BoundingBox3f(box.getCenter(), box.max);

    for(int i = 0; i < triangleCount; i++)
    {
        for(int j = 0; j < 8; j++)
        {
            if(boxlist[j].overlaps(mesh->getBoundingBox(listTriangleIndice->at(i))))
            {
                list[j].push_back(listTriangleIndice->at(i));
            }
        }
    }

    Node *node = new Node();
    node->nodeBox = new BoundingBox3f(box);
    for(int i = 0; i < 8; i++)
    {
        node->child[i] = build(boxlist[i], &list[i], mesh, depth + 1);
    }

    return node;

}

void Node::rayIntersection(Ray3f ray, std::vector<int> *hitTriangle)
{
    if(triangles != nullptr)
    {
        for(int i = 0; i < (int)triangles->size(); i++)
        {
            hitTriangle->push_back(triangles->at(i));
            // std::cout << triangles->at(i) << std::endl;
        }
    }
    //Si notre nodeBox est touche par le ray
    else if(nodeBox->rayIntersect(ray))
    {
        //On regarde s'il s'agit d'une leaf et on mets les triangles dans la liste si c'est la cas
        //On doit checker pour chaque enfant de maniere recursive
        for(int i = 0; i < 8; i++)
        {
            //On s'assure que l'enfant est existe
            if(child[i] != nullptr)
            {
                child[i]->rayIntersection(ray, hitTriangle);
            }
        }
    }
}
NORI_NAMESPACE_END