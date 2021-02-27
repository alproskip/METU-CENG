#ifndef __HW1__PARSER__
#define __HW1__PARSER__

#include <string>
#include <vector>
#include <algorithm>

namespace parser
{
    //Notice that all the structures are as simple as possible
    //so that you are not enforced to adopt any style or design.
    struct Vec3f
    {
        double x, y, z;
        Vec3f(){}
        Vec3f(const Vec3f& rhs){
            x=rhs.x;
            y=rhs.y;
            z=rhs.z;
        }
        Vec3f(double a, double b, double c){
            x=a;
            y=b;
            z=c;
        }
        Vec3f operator+(Vec3f b){  
            return Vec3f(x+b.x,y+b.y,z+b.z);
        }
        Vec3f operator-(Vec3f b){  
            return Vec3f(x-b.x,y-b.y,z-b.z);
        }
        Vec3f operator-(){  
            return Vec3f(x*-1,y*-1,z*-1);
        }
        Vec3f operator*(double d){  
            return Vec3f(x*d,y*d,z*d);
        }
        Vec3f operator/(double d){
            return Vec3f(x/d,y/d,z/d);
        }

    };

    struct Vec3i
    {
        int x, y, z;
    };


    struct Ray
    {
        Vec3f a;
        Vec3f b;
    };
    
    struct Vec4f
    {
        double x, y, z, w;
    };

    struct Camera
    {
        Vec3f position;
        Vec3f gaze;
        Vec3f up;
        Vec4f near_plane;
        double near_distance;
        int image_width, image_height;
        std::string image_name;
    };

    struct PointLight
    {
        Vec3f position;
        Vec3f intensity;
    };

    struct Material
    {
        Vec3f ambient;
        Vec3f diffuse;
        Vec3f specular;
        Vec3f mirror;
        double phong_exponent;
    };
    
    struct Face
    {
        int v0_id;
        int v1_id;
        int v2_id;
    };


    struct Triangle
    {
        int v0_id;
        int v1_id;
        int v2_id;
        int id;
        int material_id;
        Vec3f center;
        Face indices;
        Triangle(){}
        Triangle(int a,int b, int c, int mat){
            v0_id = a;
            v1_id = b;
            v2_id = c;
            indices.v0_id = a;
            indices.v1_id = b;
            indices.v2_id = c;
            material_id=mat;
        }
        bool operator<(const Triangle& rhs){
            return center.z < rhs.center.z;
        }
    };

    struct Mesh
    {
        int material_id;
        std::vector<Face> faces;
        std::vector<Triangle> tris;
        std::vector<Triangle> left;
        std::vector<Triangle> right;
        double max_x;
        double max_y;
        double max_z;
        double min_x;
        double min_y;
        double min_z;
        double max_xLEFT;
        double max_yLEFT;
        double max_zLEFT;
        double min_xLEFT;
        double min_yLEFT;
        double min_zLEFT;
        double max_xRIGHT;
        double max_yRIGHT;
        double max_zRIGHT;
        double min_xRIGHT;
        double min_yRIGHT;
        double min_zRIGHT;
    };

    struct Sphere
    {
        int material_id;
        int center_vertex_id;
        double radius;
    };

    struct Scene
    {
        //Data
        Vec3f refS;
        Vec3i background_color;
        double shadow_ray_epsilon;
        int max_recursion_depth;
        std::vector<Camera> cameras;
        Vec3f ambient_light;
        std::vector<PointLight> point_lights;
        std::vector<Material> materials;
        std::vector<Vec3f> vertex_data;
        std::vector<Mesh> meshes;
        std::vector<Triangle> triangles;
        std::vector<Sphere> spheres;
        //Functions
        void loadFromXml(const std::string& filepath);
        
        void AllMeshToTri(){
            for(int i=0; i<meshes.size(); i++){
                int len = meshes[i].faces.size();
                for(int j=0; j<meshes[i].faces.size(); j++)
                    meshes[i].tris.push_back(Triangle(meshes[i].faces[j].v0_id, meshes[i].faces[j].v1_id, meshes[i].faces[j].v2_id, meshes[i].material_id));
                if(len < 50)
                    continue;
                for(int j=0; j<meshes[i].tris.size(); j++){
                    meshes[i].tris[j].center = (vertex_data[meshes[i].tris[j].v0_id-1] + vertex_data[meshes[i].tris[j].v1_id-1] + vertex_data[meshes[i].tris[j].v2_id-1])/3.0;
                }
                std::sort(meshes[i].tris.begin(),meshes[i].tris.end());
                for(int j=0; j<meshes[i].tris.size(); j++){
                    meshes[i].tris[j].id = j;
                }
                meshes[i].left = {meshes[i].tris.begin(),meshes[i].tris.begin()+len/2};
                meshes[i].right = {meshes[i].tris.begin()+len/2+1,meshes[i].tris.end()};
            }
        }
        void calculateEdges(){
            double a, b, c;
            double d, e ,f;
            for(int i=0; i<meshes.size(); i++){
                a = b = c = 50000;
                d = e = f = -50000;
                if(meshes[i].tris.size()<50)
                    continue;
                int len = meshes[i].tris.size();
                for(int j=0; j<len; j++){
                    if(vertex_data[meshes[i].tris[j].v0_id-1].x < a)
                        a = vertex_data[meshes[i].tris[j].v0_id-1].x;
                    if(vertex_data[meshes[i].tris[j].v1_id-1].x < a)
                        a = vertex_data[meshes[i].tris[j].v1_id-1].x;
                    if(vertex_data[meshes[i].tris[j].v2_id-1].x < a)
                        a = vertex_data[meshes[i].tris[j].v2_id-1].x;
                    if(vertex_data[meshes[i].tris[j].v0_id-1].y < b)
                        b = vertex_data[meshes[i].tris[j].v0_id-1].y;
                    if(vertex_data[meshes[i].tris[j].v1_id-1].y < b)
                        b = vertex_data[meshes[i].tris[j].v1_id-1].y;
                    if(vertex_data[meshes[i].tris[j].v2_id-1].y < b)
                        b = vertex_data[meshes[i].tris[j].v2_id-1].y;
                    if(vertex_data[meshes[i].tris[j].v0_id-1].z < c)
                        c = vertex_data[meshes[i].tris[j].v0_id-1].z;
                    if(vertex_data[meshes[i].tris[j].v1_id-1].z < c)
                        c = vertex_data[meshes[i].tris[j].v1_id-1].z;
                    if(vertex_data[meshes[i].tris[j].v2_id-1].z < c)
                        c = vertex_data[meshes[i].tris[j].v2_id-1].z;
                    if(vertex_data[meshes[i].tris[j].v0_id-1].x > d)
                        d = vertex_data[meshes[i].tris[j].v0_id-1].x;
                    if(vertex_data[meshes[i].tris[j].v1_id-1].x > d)
                        d = vertex_data[meshes[i].tris[j].v1_id-1].x;
                    if(vertex_data[meshes[i].tris[j].v2_id-1].x > d)
                        d = vertex_data[meshes[i].tris[j].v2_id-1].x;
                    if(vertex_data[meshes[i].tris[j].v0_id-1].y > e)
                        e = vertex_data[meshes[i].tris[j].v0_id-1].y;
                    if(vertex_data[meshes[i].tris[j].v1_id-1].y > e)
                        e = vertex_data[meshes[i].tris[j].v1_id-1].y;
                    if(vertex_data[meshes[i].tris[j].v2_id-1].y > e)
                        e = vertex_data[meshes[i].tris[j].v2_id-1].y;                 
                    if(vertex_data[meshes[i].tris[j].v0_id-1].z > f)
                        f = vertex_data[meshes[i].tris[j].v0_id-1].z;
                    if(vertex_data[meshes[i].tris[j].v1_id-1].z > f)
                        f = vertex_data[meshes[i].tris[j].v1_id-1].z;
                    if(vertex_data[meshes[i].tris[j].v2_id-1].z > f)
                        f = vertex_data[meshes[i].tris[j].v2_id-1].z;
                }
                meshes[i].min_x = a - 0.01;
                meshes[i].min_y = b - 0.01;
                meshes[i].min_z = c - 0.01;
                meshes[i].max_x = d + 0.01;
                meshes[i].max_y = e + 0.01;
                meshes[i].max_z = f + 0.01;
                if(len < 50)
                    continue;
                int lenL = meshes[i].left.size();
                int lenR = meshes[i].right.size();
                a = b = c = 50000;
                d = e = f = -50000;
                for(int j=0; j<lenL; j++){
                    if(vertex_data[meshes[i].left[j].v0_id-1].x < a)
                        a = vertex_data[meshes[i].left[j].v0_id-1].x;
                    if(vertex_data[meshes[i].left[j].v1_id-1].x < a)
                        a = vertex_data[meshes[i].left[j].v1_id-1].x;
                    if(vertex_data[meshes[i].left[j].v2_id-1].x < a)
                        a = vertex_data[meshes[i].left[j].v2_id-1].x;
                    if(vertex_data[meshes[i].left[j].v0_id-1].y < b)
                        b = vertex_data[meshes[i].left[j].v0_id-1].y;
                    if(vertex_data[meshes[i].left[j].v1_id-1].y < b)
                        b = vertex_data[meshes[i].left[j].v1_id-1].y;
                    if(vertex_data[meshes[i].left[j].v2_id-1].y < b)
                        b = vertex_data[meshes[i].left[j].v2_id-1].y;
                    if(vertex_data[meshes[i].left[j].v0_id-1].z < c)
                        c = vertex_data[meshes[i].left[j].v0_id-1].z;
                    if(vertex_data[meshes[i].left[j].v1_id-1].z < c)
                        c = vertex_data[meshes[i].left[j].v1_id-1].z;
                    if(vertex_data[meshes[i].left[j].v2_id-1].z < c)
                        c = vertex_data[meshes[i].left[j].v2_id-1].z;
                    if(vertex_data[meshes[i].left[j].v0_id-1].x > d)
                        d = vertex_data[meshes[i].left[j].v0_id-1].x;
                    if(vertex_data[meshes[i].left[j].v1_id-1].x > d)
                        d = vertex_data[meshes[i].left[j].v1_id-1].x;
                    if(vertex_data[meshes[i].left[j].v2_id-1].x > d)
                        d = vertex_data[meshes[i].left[j].v2_id-1].x;
                    if(vertex_data[meshes[i].left[j].v0_id-1].y > e)
                        e = vertex_data[meshes[i].left[j].v0_id-1].y;
                    if(vertex_data[meshes[i].left[j].v1_id-1].y > e)
                        e = vertex_data[meshes[i].left[j].v1_id-1].y;
                    if(vertex_data[meshes[i].left[j].v2_id-1].y > e)
                        e = vertex_data[meshes[i].left[j].v2_id-1].y;                 
                    if(vertex_data[meshes[i].left[j].v0_id-1].z > f)
                        f = vertex_data[meshes[i].left[j].v0_id-1].z;
                    if(vertex_data[meshes[i].left[j].v1_id-1].z > f)
                        f = vertex_data[meshes[i].left[j].v1_id-1].z;
                    if(vertex_data[meshes[i].left[j].v2_id-1].z > f)
                        f = vertex_data[meshes[i].left[j].v2_id-1].z;
                }
                meshes[i].min_xLEFT = a - 0.01;
                meshes[i].min_yLEFT = b - 0.01;
                meshes[i].min_zLEFT = c - 0.01;
                meshes[i].max_xLEFT = d + 0.01;
                meshes[i].max_yLEFT = e + 0.01;
                meshes[i].max_zLEFT = f + 0.01;
                a = b = c = 50000;
                d = e = f = -50000;
                for(int j=0; j<lenR; j++){
                    if(vertex_data[meshes[i].right[j].v0_id-1].x < a)
                        a = vertex_data[meshes[i].right[j].v0_id-1].x;
                    if(vertex_data[meshes[i].right[j].v1_id-1].x < a)
                        a = vertex_data[meshes[i].right[j].v1_id-1].x;
                    if(vertex_data[meshes[i].right[j].v2_id-1].x < a)
                        a = vertex_data[meshes[i].right[j].v2_id-1].x;
                    if(vertex_data[meshes[i].right[j].v0_id-1].y < b)
                        b = vertex_data[meshes[i].right[j].v0_id-1].y;
                    if(vertex_data[meshes[i].right[j].v1_id-1].y < b)
                        b = vertex_data[meshes[i].right[j].v1_id-1].y;
                    if(vertex_data[meshes[i].right[j].v2_id-1].y < b)
                        b = vertex_data[meshes[i].right[j].v2_id-1].y;
                    if(vertex_data[meshes[i].right[j].v0_id-1].z < c)
                        c = vertex_data[meshes[i].right[j].v0_id-1].z;
                    if(vertex_data[meshes[i].right[j].v1_id-1].z < c)
                        c = vertex_data[meshes[i].right[j].v1_id-1].z;
                    if(vertex_data[meshes[i].right[j].v2_id-1].z < c)
                        c = vertex_data[meshes[i].right[j].v2_id-1].z;
                    if(vertex_data[meshes[i].right[j].v0_id-1].x > d)
                        d = vertex_data[meshes[i].right[j].v0_id-1].x;
                    if(vertex_data[meshes[i].right[j].v1_id-1].x > d)
                        d = vertex_data[meshes[i].right[j].v1_id-1].x;
                    if(vertex_data[meshes[i].right[j].v2_id-1].x > d)
                        d = vertex_data[meshes[i].right[j].v2_id-1].x;
                    if(vertex_data[meshes[i].right[j].v0_id-1].y > e)
                        e = vertex_data[meshes[i].right[j].v0_id-1].y;
                    if(vertex_data[meshes[i].right[j].v1_id-1].y > e)
                        e = vertex_data[meshes[i].right[j].v1_id-1].y;
                    if(vertex_data[meshes[i].right[j].v2_id-1].y > e)
                        e = vertex_data[meshes[i].right[j].v2_id-1].y;                 
                    if(vertex_data[meshes[i].right[j].v0_id-1].z > f)
                        f = vertex_data[meshes[i].right[j].v0_id-1].z;
                    if(vertex_data[meshes[i].right[j].v1_id-1].z > f)
                        f = vertex_data[meshes[i].right[j].v1_id-1].z;
                    if(vertex_data[meshes[i].right[j].v2_id-1].z > f)
                        f = vertex_data[meshes[i].right[j].v2_id-1].z;
                }
                meshes[i].min_xRIGHT = a - 0.01;
                meshes[i].min_yRIGHT = b - 0.01;
                meshes[i].min_zRIGHT = c - 0.01;
                meshes[i].max_xRIGHT = d + 0.01;
                meshes[i].max_yRIGHT = e + 0.01;
                meshes[i].max_zRIGHT = f + 0.01;
            }
        }
    };
}

#endif