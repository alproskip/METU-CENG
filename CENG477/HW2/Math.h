#include<cmath>
#include<cfloat>
#include<iostream>
#include"parser.h"
using namespace parser;

int intersectMesh(Ray& r, Mesh& mesh, bool forpixel, Scene& sc);

double DOT(Vec3f& A, Vec3f& B){
    return A.x*B.x+A.y*B.y+A.z*B.z;
}

double distance(Vec3f& A, Vec3f& B)
{
    return sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
}

double length(Vec3f v)
{
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
}

Vec3f cross(Vec3f a, Vec3f b)
{
    return Vec3f(a.y*b.z-b.y*a.z,b.x*a.z-a.x*b.z,a.x*b.y-b.x*a.y);
}

Vec3f normalize(Vec3f v)
{
    double d = length(v);
    return Vec3f(v.x/d,v.y/d,v.z/d);
}

void translate(Vec3f& v, const Vec3f& t){
    v.x += t.x;
    v.y += t.y;
    v.z += t.z;
}

void scale(Vec3f& v, Vec3f& t){
    v.x *= t.x;
    v.y *= t.y;
    v.z *= t.z;
}

void x_rotate(Vec3f& v, double angle){
    double temp = 0;
    angle *= 0.0174533;
    double c = cos(angle);
    double s = sin(angle);

    if (c<0.000001 && c>-0.000001) c = 0;
    if (s<0.000001 && s>-0.000001) s = 0;

    v.x = v.x;
    temp = v.y*c-v.z*s;
    v.z = v.y*s+v.z*c; 
    v.y = temp;
}

void y_rotate(Vec3f& v, double angle){
    double temp = 0;
    angle *= 0.0174533;
    double c = cos(angle);
    double s = sin(angle);

    if (c<0.000001 && c>-0.000001) c = 0;
    if (s<0.000001 && s>-0.000001) s = 0;

    temp = v.x*c + v.z*s;
    v.y = v.y;
    v.z = v.z*c - v.x*s;
    v.x = temp;
}

void z_rotate(Vec3f& v, double angle){
    double temp = 0;
    angle *= 0.0174533;
    double c = cos(angle);
    double s = sin(angle);

    if (c<0.000001 && c>-0.000001) c = 0;
    if (s<0.000001 && s>-0.000001) s = 0;
    
    temp = v.x*c - v.y*s;
    v.y = v.x*s + v.y*c;
    v.z = v.z;
    v.x = temp;
}

void rotate(Vec3f& vertex, Rotation& r){
    double min, mid, max;
    Vec3f u, v, w;
    if (r.x == 1 and r.y == 0 and r.z == 0) // Rotate to X
        x_rotate(vertex, r.angle);

    else if (r.x == 0 and r.y == 1 and r.z == 0) // Rotate to Y
        y_rotate(vertex, r.angle);
    
    else if (r.x == 0 and r.y == 0 and r.z == 1) // Rotate to Z
        z_rotate(vertex, r.angle);
    
    else {
        // Rotate the Axis
        // create orthonormal basis
        u = Vec3f(r.x, r.y, r.z);
        min = std::min(u.x, std::min(u.y,u.z));
        if (min == u.x) {
            v = Vec3f(0, u.z, -1*u.y);
            w = normalize(cross(u,v));
            v = normalize(v);
        }
        else if (min == u.y) {
            v = Vec3f(v.z, 0, -1*v.x);
            w = normalize(cross(u,v));
            v = normalize(v);
        }
        else if (min == u.z) {
            v = Vec3f(v.y, -1*v.x, 0);
            w = normalize(cross(u,v));
            v = normalize(v);
        }

        // Find center in orthonormal basis
        Vec3f temp = vertex;
        temp.x = u.x*temp.x + u.y*temp.y + u.z*temp.z;
        temp.y = v.x*temp.x + v.y*temp.y + v.z*temp.z;
        temp.z = w.x*temp.x + w.y*temp.y + w.z*temp.z;
        // Do the rotation
        x_rotate(temp, r.angle);
        // Turn it back to xyz coordinates
        temp.x = u.x*temp.x + v.x*temp.y + w.x*temp.z;
        temp.y = u.y*temp.x + v.y*temp.y + w.y*temp.z;
        temp.z = u.z*temp.x + v.z*temp.y + w.z*temp.z;
        // Set new center of the sphere
        vertex = temp;
    }
}

void split(std::string& trans, char& t_flag, int& ind){
    t_flag = trans[0];
    trans.erase(trans.begin());
    std::string toint = "";
    while(!trans.empty() && trans[0]!=' '){
        toint.push_back(trans[0]);
        trans.erase(trans.begin());
    }
    if(!trans.empty())
        trans.erase(trans.begin());
    ind = std::stoi(toint);
}

void split(std::string& trans, std::vector<int>& rotind){
    std::string temp = trans;
    std::string toint = "";
    while(!temp.empty()){
        if(temp[0] == 's' || temp[0] == 't'){
            while(!temp.empty() && temp[0]!=' ')
                temp.erase(temp.begin());
            if(!temp.empty())
                temp.erase(temp.begin());
        }
        else if(temp[0] == 'r'){
            temp.erase(temp.begin());
            while(!temp.empty() && temp[0]!=' '){
                toint.push_back(temp[0]);
                temp.erase(temp.begin());
            }
            if(!temp.empty())
                temp.erase(temp.begin());
            rotind.push_back(std::stoi(toint));
            toint = "";
        }
        else
            temp.erase(temp.begin());
    }
}

void transform(Scene& sc){
    int tlen = sc.triangles.size();
    
    for(int i = 0 ; i<tlen; i++){

        std::string trans = sc.triangles[i].transformations; 
        
        while(!trans.empty()){
            char t_flag; //What type of transformation?
            int ind; // get transformation index
            split(trans,t_flag,ind);
            if(t_flag == 't'){
                translate(sc.vertex_data[sc.triangles[i].indices.v0_id-1], sc.translations[ind-1]);
                translate(sc.vertex_data[sc.triangles[i].indices.v1_id-1], sc.translations[ind-1]);
                translate(sc.vertex_data[sc.triangles[i].indices.v2_id-1], sc.translations[ind-1]);
            }
            else if(t_flag == 'r'){
                rotate(sc.vertex_data[sc.triangles[i].indices.v0_id-1], sc.rotations[ind-1]);
                rotate(sc.vertex_data[sc.triangles[i].indices.v1_id-1], sc.rotations[ind-1]);
                rotate(sc.vertex_data[sc.triangles[i].indices.v2_id-1], sc.rotations[ind-1]);
            }
            else if(t_flag == 's'){
                scale(sc.vertex_data[sc.triangles[i].indices.v0_id-1], sc.scalings[ind-1]);
                scale(sc.vertex_data[sc.triangles[i].indices.v1_id-1], sc.scalings[ind-1]);
                scale(sc.vertex_data[sc.triangles[i].indices.v2_id-1], sc.scalings[ind-1]);
            }
        }
        
    }
    
    int slen = sc.spheres.size();
    
    for(int i = 0 ; i<slen; i++){

        std::string trans = sc.spheres[i].transformations;

        while(!trans.empty()){
            char t_flag; //What type of transformation?
            int ind; // get transformation index
            split(trans,t_flag,ind);

            if(t_flag == 't'){
                translate(sc.spheres[i].center, sc.translations[ind-1]);
            }

            else if(t_flag == 'r'){
                rotate(sc.spheres[i].center, sc.rotations[ind-1]);
            }

            else if(t_flag == 's'){
                scale(sc.spheres[i].center, sc.scalings[ind-1]);
                sc.spheres[i].radius *= sc.scalings[ind-1].x;
            }

        }
    }
    for(int i = 0 ; i<sc.meshes.size(); i++){

        std::string trans = sc.meshes[i].transformations;

        while(!trans.empty()){
            char t_flag; // What type of transformation?
            int ind; // get transformation index
            split(trans,t_flag,ind);
            if(t_flag == 't'){
                bool vertexToTranslate[sc.vertex_data.size()] = {0};
                for(int j = 0; j<sc.meshes[i].faces.size(); j++){
                    if(!vertexToTranslate[sc.meshes[i].faces[j].v0_id-1]){
                        vertexToTranslate[sc.meshes[i].faces[j].v0_id-1] = true;
                        translate(sc.vertex_data[sc.meshes[i].faces[j].v0_id-1], sc.translations[ind-1]);
                    }
                    if(!vertexToTranslate[sc.meshes[i].faces[j].v1_id-1]){
                        vertexToTranslate[sc.meshes[i].faces[j].v1_id-1] = true;
                        translate(sc.vertex_data[sc.meshes[i].faces[j].v1_id-1], sc.translations[ind-1]);
                    }
                    if(!vertexToTranslate[sc.meshes[i].faces[j].v2_id-1]){
                        vertexToTranslate[sc.meshes[i].faces[j].v2_id-1] = true;
                        translate(sc.vertex_data[sc.meshes[i].faces[j].v2_id-1], sc.translations[ind-1]);
                    }
                }
            }
            
            else if(t_flag == 'r'){
                bool vertexToRotate[sc.vertex_data.size()] = {0};
                for(int j = 0; j<sc.meshes[i].faces.size(); j++){
                    if(!vertexToRotate[sc.meshes[i].faces[j].v0_id-1]){
                        vertexToRotate[sc.meshes[i].faces[j].v0_id-1] = true;
                        rotate(sc.vertex_data[sc.meshes[i].faces[j].v0_id-1], sc.rotations[ind-1]);
                    }
                    if(!vertexToRotate[sc.meshes[i].faces[j].v1_id-1]){
                        vertexToRotate[sc.meshes[i].faces[j].v1_id-1] = true;
                        rotate(sc.vertex_data[sc.meshes[i].faces[j].v1_id-1], sc.rotations[ind-1]);
                    }
                    if(!vertexToRotate[sc.meshes[i].faces[j].v2_id-1]){
                        vertexToRotate[sc.meshes[i].faces[j].v2_id-1] = true;
                        rotate(sc.vertex_data[sc.meshes[i].faces[j].v2_id-1], sc.rotations[ind-1]);
                    }
                }
            }
            else if(t_flag == 's'){
                bool vertexToScale[sc.vertex_data.size()] = {0};
                for(int j = 0; j<sc.meshes[i].faces.size(); j++){
                    if(!vertexToScale[sc.meshes[i].faces[j].v0_id-1]){
                        vertexToScale[sc.meshes[i].faces[j].v0_id-1] = true;
                        scale(sc.vertex_data[sc.meshes[i].faces[j].v0_id-1], sc.scalings[ind-1]);
                    }
                    if(!vertexToScale[sc.meshes[i].faces[j].v1_id-1]){
                        vertexToScale[sc.meshes[i].faces[j].v1_id-1] = true;
                        scale(sc.vertex_data[sc.meshes[i].faces[j].v1_id-1], sc.scalings[ind-1]);
                    }
                    if(!vertexToScale[sc.meshes[i].faces[j].v2_id-1]){
                        vertexToScale[sc.meshes[i].faces[j].v2_id-1] = true;
                        scale(sc.vertex_data[sc.meshes[i].faces[j].v2_id-1], sc.scalings[ind-1]);
                    }
                }
            }
        }
    } 
}

double intersectSphere(Ray& r, Sphere& s, bool forpixel, Scene& sc)
{
    double A,B,C;    
    double delta;    
    Vec3f scenter;
    double sradius;    
    Vec3f p;        
    double t,t1,t2;
    int i;
    if(forpixel)
        r.b = r.b - r.a;
    scenter = s.center;
    sradius = s.radius;    
    C = (r.a.x-scenter.x)*(r.a.x-scenter.x)+(r.a.y-scenter.y)*(r.a.y-scenter.y)+(r.a.z-scenter.z)*(r.a.z-scenter.z)-sradius*sradius;
    B = 2*r.b.x*(r.a.x-scenter.x)+2*r.b.y*(r.a.y-scenter.y)+2*r.b.z*(r.a.z-scenter.z); 
    A = r.b.x*r.b.x+r.b.y*r.b.y+r.b.z*r.b.z;    
    delta = B*B-4*A*C;  
    if(forpixel)
        r.b = r.b + r.a;  
    if (delta<0) return -1;
    else if (delta==0)
    {
        t = -B / (2*A);
    }
    else
    {
        double tmp;
        delta = sqrt(delta);
        A = 2*A;
        t1 = (-B + delta) / A;
        t2 = (-B - delta) / A;    
        if (t2<t1){
            tmp = t2;
            t2 = t1;
            t1 = tmp;
        }
        t = t1;
    }
    return t;
}

double intersectTriangle(Ray& r, Triangle& tr, bool forpixel, Scene& sc)
{
    double  a,b,c,d,e,f,g,h,i,j,k,l;
    double beta,gamma,t;
    double eimhf,gfmdi,dhmeg,akmjb,jcmal,blmkc;
    double M;
    double dd;
    if(forpixel)
        r.b = r.b - r.a;
    Vec3f ma,mb,mc;
    ma = sc.vertex_data[tr.indices.v1_id-1];
    mb = sc.vertex_data[tr.indices.v2_id-1];
    mc = sc.vertex_data[tr.indices.v0_id-1];
    a = ma.x-mb.x;
    b = ma.y-mb.y;
    c = ma.z-mb.z;
    d = ma.x-mc.x;
    e = ma.y-mc.y;
    f = ma.z-mc.z;
    g = r.b.x;
    h = r.b.y;
    i = r.b.z;
    j = ma.x-r.a.x;
    k = ma.y-r.a.y;
    l = ma.z-r.a.z; 
    eimhf = e*i-h*f;
    gfmdi = g*f-d*i;
    dhmeg = d*h-e*g;
    akmjb = a*k-j*b;
    jcmal = j*c-a*l;
    blmkc = b*l-k*c;
    M = a*eimhf+b*gfmdi+c*dhmeg;
    if(forpixel)
        r.b = r.b + r.a;

    if (M==0) return -1;
    
    t = -(f*akmjb+e*jcmal+d*blmkc)/M;
        
    gamma = (i*akmjb+h*jcmal+g*blmkc)/M;
    
    if (gamma<0 || gamma>1) return -1;
    
    beta = (j*eimhf+k*gfmdi+l*dhmeg)/M;
    
    if (beta<0 || beta>(1.0000000001-gamma)) return -1; 

    return t;
}

double intersectTriangle(Ray& r, Triangle& tr, bool forpixel, double& BETA, double& GAMMA, Scene& sc)
{
    double  a,b,c,d,e,f,g,h,i,j,k,l;
    double beta,gamma,t;
    double eimhf,gfmdi,dhmeg,akmjb,jcmal,blmkc;
    double M;
    double dd;
    if(forpixel)
        r.b = r.b - r.a;
    Vec3f ma,mb,mc;
    ma = sc.vertex_data[tr.indices.v1_id-1];
    mb = sc.vertex_data[tr.indices.v2_id-1];
    mc = sc.vertex_data[tr.indices.v0_id-1];
    a = ma.x-mb.x;
    b = ma.y-mb.y;
    c = ma.z-mb.z;
    d = ma.x-mc.x;
    e = ma.y-mc.y;
    f = ma.z-mc.z;
    g = r.b.x;
    h = r.b.y;
    i = r.b.z;
    j = ma.x-r.a.x;
    k = ma.y-r.a.y;
    l = ma.z-r.a.z; 
    eimhf = e*i-h*f;
    gfmdi = g*f-d*i;
    dhmeg = d*h-e*g;
    akmjb = a*k-j*b;
    jcmal = j*c-a*l;
    blmkc = b*l-k*c;
    M = a*eimhf+b*gfmdi+c*dhmeg;
    if(forpixel)
        r.b = r.b + r.a;

    if (M==0) return -1;
    
    t = -(f*akmjb+e*jcmal+d*blmkc)/M;
        
    gamma = (i*akmjb+h*jcmal+g*blmkc)/M;
    
    if (gamma<0 || gamma>1) return -1;
    
    beta = (j*eimhf+k*gfmdi+l*dhmeg)/M;
    
    if (beta<0 || beta>(1.0000000001-gamma)) return -1; 

    BETA = beta;
    GAMMA = gamma;
    return t;
}

Vec3f computeTriangleNormal(Triangle& tri, Scene& sc){
    return normalize(cross(sc.vertex_data[tri.indices.v1_id-1]-sc.vertex_data[tri.indices.v0_id-1], sc.vertex_data[tri.indices.v2_id-1]-sc.vertex_data[tri.indices.v1_id-1]));
}

Vec3f computeSphereNormal(Sphere& sp, Vec3f& point, Scene& sc){
    return normalize(point - sp.center); // sc.vertex_data[sp.center_vertex_id-1]);
}

void rayHitsWhat(Ray& ray, int& triorsphere, int& itsindex, int& moreindex, double& t, bool pixel , double& BETA, double& GAMMA, Scene& scene){
    double temp_t, t_min=t;
    int numTriangles = scene.triangles.size();
    int numSphere = scene.spheres.size();
    int numMesh = scene.meshes.size();
    for (int tri_index = 0; tri_index < numTriangles; tri_index++){
        temp_t = intersectTriangle(ray, scene.triangles[tri_index], pixel, BETA, GAMMA, scene);
        if (temp_t>0 && !pixel || pixel && temp_t>=1){
            if (temp_t < t_min){
                t_min = temp_t;
                triorsphere = 1;
                itsindex = tri_index;
            }
        }
    }
    for(int s_index=0; s_index <numSphere ; s_index++){
        temp_t = intersectSphere(ray,scene.spheres[s_index], pixel, scene);
        if (temp_t>0 && !pixel || pixel && temp_t>=1){
            if (temp_t < t_min){
                t_min = temp_t;
                triorsphere = 2;
                itsindex = s_index;
            }
        }
    }
    for(int m=0; m<numMesh; m++){
        if(scene.meshes[m].tris.size()<50){
            temp_t = 3;
        }
        else
            temp_t = intersectMesh(ray,scene.meshes[m], pixel, scene);
        if(temp_t==1){
            for (int mt = 0; mt < scene.meshes[m].left.size(); mt++){
                temp_t = intersectTriangle(ray, scene.meshes[m].left[mt], pixel, BETA, GAMMA, scene);
                if (temp_t>0 && !pixel || pixel && temp_t>=1){
                    if (temp_t < t_min){
                        t_min = temp_t;
                        triorsphere = 3;
                        itsindex = m;
                        moreindex = scene.meshes[m].left[mt].id;
                    }
                }
            } 
        }
        else if(temp_t==2){
            for (int mt = 0; mt < scene.meshes[m].right.size(); mt++){
                temp_t = intersectTriangle(ray, scene.meshes[m].right[mt], pixel, BETA, GAMMA, scene);
                if (temp_t>0 && !pixel || pixel && temp_t>=1){
                    if (temp_t < t_min){
                        t_min = temp_t;
                        triorsphere = 3;
                        itsindex = m;
                        moreindex = scene.meshes[m].right[mt].id;
                    }
                }
            } 
        }
        else if(temp_t==3){
            for (int mt = 0; mt < scene.meshes[m].tris.size(); mt++){
                temp_t = intersectTriangle(ray, scene.meshes[m].tris[mt], pixel,BETA, GAMMA, scene);
                if (temp_t>0 && !pixel || pixel && temp_t>=1){
                    if (temp_t < t_min){
                        t_min = temp_t;
                        triorsphere = 3;
                        itsindex = m;
                        moreindex = mt;
                    }
                }
            } 
        }
    }
    t = t_min;
}
//u(β, γ) = ua + β(ub – ua) + γ(uc – ua)
//v(β, γ) = va + β(vb – va) + γ(vc – va)
Vec3f getKD(double& BETA, double& GAMMA, Triangle& t ,Scene& sc){
    double KR, KG, KB;
    double u,v;
    u = sc.tex_coord_data[t.indices.v0_id-1].x + \
    BETA*(sc.tex_coord_data[t.indices.v1_id-1].x - sc.tex_coord_data[t.indices.v0_id-1].x) + \
    GAMMA*(sc.tex_coord_data[t.indices.v2_id-1].x - sc.tex_coord_data[t.indices.v0_id-1].x);
    v = sc.tex_coord_data[t.indices.v0_id-1].y + \
    BETA*(sc.tex_coord_data[t.indices.v1_id-1].y - sc.tex_coord_data[t.indices.v0_id-1].y) + \
    GAMMA*(sc.tex_coord_data[t.indices.v2_id-1].y - sc.tex_coord_data[t.indices.v0_id-1].y);
    //u = ubarycentric – floor(ubarycentric)
    //v = vbarycentric – floor(vbarycentric)
    u *= abs(sc.tex_coord_data[t.indices.v1_id-1].x - sc.tex_coord_data[t.indices.v0_id-1].x );
    v *= abs(sc.tex_coord_data[t.indices.v2_id-1].y - sc.tex_coord_data[t.indices.v0_id-1].y );
    if(sc.textures[t.texture_id-1].interpolation=="nearest"){
        int x = (u-floor(u)) * sc.textures[t.texture_id-1].width;
        int y = (v-floor(v)) * sc.textures[t.texture_id-1].height;
        KR = sc.textures[t.texture_id-1].imgarray[3*(x+sc.textures[t.texture_id-1].width*y)];
        KG = sc.textures[t.texture_id-1].imgarray[3*(x+sc.textures[t.texture_id-1].width*y)+1];
        KB = sc.textures[t.texture_id-1].imgarray[3*(x+sc.textures[t.texture_id-1].width*y)+2];
    }
    else{
        int p = floor(u);
        int q = floor(v);
        double dx = u-p;
        double dy = v-q;
        KR = sc.textures[t.texture_id-1].imgarray[3*(p+sc.textures[t.texture_id-1].width*q)]*(1-dx)*(1-dy) + \
             sc.textures[t.texture_id-1].imgarray[3*(p+1+sc.textures[t.texture_id-1].width*q)]*(dx)*(1-dy) + \
             sc.textures[t.texture_id-1].imgarray[3*(p+sc.textures[t.texture_id-1].width*(q+1))]*(1-dx)*(dy) + \
             sc.textures[t.texture_id-1].imgarray[3*(p+1+sc.textures[t.texture_id-1].width*(q+1))]*(dx)*(dy);

        KG = sc.textures[t.texture_id-1].imgarray[3*(p+sc.textures[t.texture_id-1].width*q)+1]*(1-dx)*(1-dy) + \
             sc.textures[t.texture_id-1].imgarray[3*(p+1+sc.textures[t.texture_id-1].width*q)+1]*(dx)*(1-dy) + \
             sc.textures[t.texture_id-1].imgarray[3*(p+sc.textures[t.texture_id-1].width*(q+1))+1]*(1-dx)*(dy) + \
             sc.textures[t.texture_id-1].imgarray[3*(p+1+sc.textures[t.texture_id-1].width*(q+1))+1]*(dx)*(dy);

        KB = sc.textures[t.texture_id-1].imgarray[3*(p+sc.textures[t.texture_id-1].width*q)+2]*(1-dx)*(1-dy) + \
             sc.textures[t.texture_id-1].imgarray[3*(p+1+sc.textures[t.texture_id-1].width*q)+2]*(dx)*(1-dy) + \
             sc.textures[t.texture_id-1].imgarray[3*(p+sc.textures[t.texture_id-1].width*(q+1))+2]*(1-dx)*(dy) + \
             sc.textures[t.texture_id-1].imgarray[3*(p+1+sc.textures[t.texture_id-1].width*(q+1))+2]*(dx)*(dy);       
    }

    return Vec3f(KR,KG,KB);
}

Vec3f getKD(Vec3f& p, Sphere& s, Scene& sc){
    std::cout<<"DAYI \n"<<std::flush;
    double KR, KG, KB;
    Vec3f temppoint = p;
    temppoint = temppoint - s.center;

    std::vector<int> rotatesVec;
    split(s.transformations,rotatesVec);
    for(int i = 0; i< rotatesVec.size(); i++){
        Rotation tmp;
        tmp = sc.rotations[rotatesVec[i]-1];
        tmp.angle = sc.rotations[rotatesVec[i]-1].angle * -1 ;
        rotate(temppoint,tmp);
    }
    temppoint = temppoint + s.center;

    double theta =  acos( (temppoint.y-s.center.y)/s.radius );
    double phi = atan2( (temppoint.z-s.center.z) , (temppoint.x-s.center.x) );
    //double theta =  acos( (temppoint.z-s.center.z)/s.radius );
    //double phi = atan2( (temppoint.x-s.center.x) , (temppoint.y-s.center.y) );
    double u = (-phi+3.141592653589793)/(2*3.141592653589793);
    double v = theta/3.141592653589793;
    //std::cout << s.texture_id <<std::flush;
    //std::cout << sc.textures[s.texture_id-1].decalMode <<std::flush;
    if(sc.textures[s.texture_id-1].appearance=="clamp"){
        u = std::fmax(0,u);
        u = std::fmin(1,u);
        v = std::fmax(0,v);   
        v = std::fmin(1,v);        
    }

    //std::cout<<"YEGEN \n"<<std::flush; 

    if(sc.textures[s.texture_id-1].interpolation=="nearest"){
        int x = (u-floor(u)) * sc.textures[s.texture_id-1].width;
        int y = (v-floor(v)) * sc.textures[s.texture_id-1].height;
        KR = sc.textures[s.texture_id-1].imgarray[ 3* (x+ sc.textures[s.texture_id-1].width*y)  ];
        KG = sc.textures[s.texture_id-1].imgarray[ 3* (x+ sc.textures[s.texture_id-1].width*y)+1];
        KB = sc.textures[s.texture_id-1].imgarray[ 3* (x+ sc.textures[s.texture_id-1].width*y)+2];
    }
    else{
        double AA = (u-floor(u)) * sc.textures[s.texture_id-1].width;
        double BB = (v-floor(v)) * sc.textures[s.texture_id-1].height;
        int x = AA;
        int y = BB;
        double dx = AA-x;
        double dy = BB-y;
        KR = sc.textures[s.texture_id-1].imgarray[3 * (x   +sc.textures[s.texture_id-1].width * y)    ] * (1-dx) * (1-dy) + \
             sc.textures[s.texture_id-1].imgarray[3 * (x+1 +sc.textures[s.texture_id-1].width * y)    ] * (dx)   * (1-dy) + \
             sc.textures[s.texture_id-1].imgarray[3 * (x   +sc.textures[s.texture_id-1].width * (y+1))] * (1-dx) * (dy)   + \
             sc.textures[s.texture_id-1].imgarray[3 * (x+1 +sc.textures[s.texture_id-1].width * (y+1))] * (dx)   * (dy);

        KG = sc.textures[s.texture_id-1].imgarray[3 * (x   +sc.textures[s.texture_id-1].width * y)+1]     * (1-dx) * (1-dy) + \
             sc.textures[s.texture_id-1].imgarray[3 * (x+1 +sc.textures[s.texture_id-1].width * y)+1]     * (dx)   * (1-dy) + \
             sc.textures[s.texture_id-1].imgarray[3 * (x   +sc.textures[s.texture_id-1].width * (y+1))+1] * (1-dx) * (dy)   + \
             sc.textures[s.texture_id-1].imgarray[3 * (x+1 +sc.textures[s.texture_id-1].width * (y+1))+1] * (dx)   * (dy);
             
        KB = sc.textures[s.texture_id-1].imgarray[3 * (x   +sc.textures[s.texture_id-1].width * y)+2]     * (1-dx) * (1-dy) + \
             sc.textures[s.texture_id-1].imgarray[3 * (x+1 +sc.textures[s.texture_id-1].width * y)+2]     * (dx)   * (1-dy) + \
             sc.textures[s.texture_id-1].imgarray[3 * (x   +sc.textures[s.texture_id-1].width * (y+1))+2] * (1-dx) * (dy)   + \
             sc.textures[s.texture_id-1].imgarray[3 * (x+1 +sc.textures[s.texture_id-1].width * (y+1))+2] * (dx)   * (dy); 

    }

    return Vec3f(KR,KG,KB);
}

Vec3f getKD(double& BETA, double& GAMMA, Mesh& m, int moreindex, Scene& sc){
    double KR, KG, KB;
    double u,v;
    u = sc.tex_coord_data[m.tris[moreindex].indices.v1_id-1].x + \
    BETA*(sc.tex_coord_data[m.tris[moreindex].indices.v2_id-1].x - sc.tex_coord_data[m.tris[moreindex].indices.v1_id-1].x) + \
    GAMMA*(sc.tex_coord_data[m.tris[moreindex].indices.v0_id-1].x - sc.tex_coord_data[m.tris[moreindex].indices.v1_id-1].x);
    v = sc.tex_coord_data[m.tris[moreindex].indices.v1_id-1].y + \
    BETA*(sc.tex_coord_data[m.tris[moreindex].indices.v2_id-1].y - sc.tex_coord_data[m.tris[moreindex].indices.v1_id-1].y) + \
    GAMMA*(sc.tex_coord_data[m.tris[moreindex].indices.v0_id-1].y - sc.tex_coord_data[m.tris[moreindex].indices.v1_id-1].y);

    if(sc.textures[m.texture_id-1].appearance=="Clamp"){
        u = std::fmax(0,u);
        u = std::fmin(1,u);
        v = std::fmax(0,v);   
        v = std::fmin(1,v);        
    }
    else{
        u *= sqrt( (sc.tex_coord_data[m.tris[moreindex].indices.v2_id-1].x - sc.tex_coord_data[m.tris[moreindex].indices.v1_id-1].x) * \
                   (sc.tex_coord_data[m.tris[moreindex].indices.v2_id-1].x - sc.tex_coord_data[m.tris[moreindex].indices.v1_id-1].x) + \
                   (sc.tex_coord_data[m.tris[moreindex].indices.v2_id-1].y - sc.tex_coord_data[m.tris[moreindex].indices.v1_id-1].y) * \
                   (sc.tex_coord_data[m.tris[moreindex].indices.v2_id-1].y - sc.tex_coord_data[m.tris[moreindex].indices.v1_id-1].y) );
        v *= sqrt( (sc.tex_coord_data[m.tris[moreindex].indices.v0_id-1].x - sc.tex_coord_data[m.tris[moreindex].indices.v1_id-1].x) * \
                   (sc.tex_coord_data[m.tris[moreindex].indices.v0_id-1].x - sc.tex_coord_data[m.tris[moreindex].indices.v1_id-1].x) + \
                   (sc.tex_coord_data[m.tris[moreindex].indices.v0_id-1].y - sc.tex_coord_data[m.tris[moreindex].indices.v1_id-1].y) * \
                   (sc.tex_coord_data[m.tris[moreindex].indices.v0_id-1].y - sc.tex_coord_data[m.tris[moreindex].indices.v1_id-1].y) );
    }


    if(sc.textures[m.texture_id-1].interpolation=="nearest"){
        int x = (u-floor(u)) * sc.textures[m.texture_id-1].width;
        int y = (v-floor(v)) * sc.textures[m.texture_id-1].height;
        KR = sc.textures[m.texture_id-1].imgarray[ 3* (x+ sc.textures[m.texture_id-1].width*y)  ];
        KG = sc.textures[m.texture_id-1].imgarray[ 3* (x+ sc.textures[m.texture_id-1].width*y)+1];
        KB = sc.textures[m.texture_id-1].imgarray[ 3* (x+ sc.textures[m.texture_id-1].width*y)+2];
    }
    else{
        double AA = (u-floor(u)) * sc.textures[m.texture_id-1].width;
        double BB = (v-floor(v)) * sc.textures[m.texture_id-1].height;
        int x = AA;
        int y = BB;
        double dx = AA-x;
        double dy = BB-y;
        KR = sc.textures[m.texture_id-1].imgarray[3 * (x   +sc.textures[m.texture_id-1].width * y)    ] * (1-dx) * (1-dy) + \
             sc.textures[m.texture_id-1].imgarray[3 * (x+1 +sc.textures[m.texture_id-1].width * y)    ] * (dx)   * (1-dy) + \
             sc.textures[m.texture_id-1].imgarray[3 * (x   +sc.textures[m.texture_id-1].width * (y+1))] * (1-dx) * (dy)   + \
             sc.textures[m.texture_id-1].imgarray[3 * (x+1 +sc.textures[m.texture_id-1].width * (y+1))] * (dx)   * (dy);

        KG = sc.textures[m.texture_id-1].imgarray[3 * (x   +sc.textures[m.texture_id-1].width * y)+1]     * (1-dx) * (1-dy) + \
             sc.textures[m.texture_id-1].imgarray[3 * (x+1 +sc.textures[m.texture_id-1].width * y)+1]     * (dx)   * (1-dy) + \
             sc.textures[m.texture_id-1].imgarray[3 * (x   +sc.textures[m.texture_id-1].width * (y+1))+1] * (1-dx) * (dy)   + \
             sc.textures[m.texture_id-1].imgarray[3 * (x+1 +sc.textures[m.texture_id-1].width * (y+1))+1] * (dx)   * (dy);
             
        KB = sc.textures[m.texture_id-1].imgarray[3 * (x   +sc.textures[m.texture_id-1].width * y)+2]     * (1-dx) * (1-dy) + \
             sc.textures[m.texture_id-1].imgarray[3 * (x+1 +sc.textures[m.texture_id-1].width * y)+2]     * (dx)   * (1-dy) + \
             sc.textures[m.texture_id-1].imgarray[3 * (x   +sc.textures[m.texture_id-1].width * (y+1))+2] * (1-dx) * (dy)   + \
             sc.textures[m.texture_id-1].imgarray[3 * (x+1 +sc.textures[m.texture_id-1].width * (y+1))+2] * (dx)   * (dy);       
    }
    return Vec3f(KR,KG,KB);
}

void calculateColorOfPoint(Vec3f& from, Vec3f point, int* color, int which, int itsindex, int moreindex, double& BETA, double& GAMMA, Scene& scene){
    bool mapping = false;
    Vec3f THEKD;
    if(which==1){
        if(scene.triangles[itsindex].texture_id != -1){
            THEKD = getKD(BETA,GAMMA,scene.triangles[itsindex],scene);
            mapping = true;
        }
    }
    else if(which==2){
        if(scene.spheres[itsindex].texture_id != -1){
            THEKD = getKD(point,scene.spheres[itsindex],scene);
            mapping = true;
        } 
    }
    else if(which==3){
        if(scene.meshes[itsindex].texture_id != -1){
            THEKD = getKD(BETA,GAMMA,scene.meshes[itsindex],moreindex,scene);
            mapping = true;
        }
    }

    double tL, t_minL;
    int intersect;
    double am1,am2,am3;
    double dff1,dff2,dff3;
    double bp1,bp2,bp3;
    int dummy,dummy2;
    Ray tolight; // direction
    if(which==1){
        color[0] = scene.ambient_light.x * scene.materials[scene.triangles[itsindex].material_id-1].ambient.x;
        color[1] = scene.ambient_light.y * scene.materials[scene.triangles[itsindex].material_id-1].ambient.y;
        color[2] = scene.ambient_light.z * scene.materials[scene.triangles[itsindex].material_id-1].ambient.z;         
    }
    else if(which==2){
        color[0] = scene.ambient_light.x * scene.materials[scene.spheres[itsindex].material_id-1].ambient.x;
        color[1] = scene.ambient_light.y * scene.materials[scene.spheres[itsindex].material_id-1].ambient.y;
        color[2] = scene.ambient_light.z * scene.materials[scene.spheres[itsindex].material_id-1].ambient.z;
    }
    else if(which==3){
        color[0] = scene.ambient_light.x * scene.materials[scene.meshes[itsindex].tris[moreindex].material_id-1].ambient.x;
        color[1] = scene.ambient_light.y * scene.materials[scene.meshes[itsindex].tris[moreindex].material_id-1].ambient.y;
        color[2] = scene.ambient_light.z * scene.materials[scene.meshes[itsindex].tris[moreindex].material_id-1].ambient.z;
    }

    Vec3f NORMAL;
    if(which==1){
        NORMAL = computeTriangleNormal(scene.triangles[itsindex],scene);
    }
    else if(which==2){
        NORMAL = computeSphereNormal(scene.spheres[itsindex], point, scene);
    }
    else if(which==3){
        NORMAL = computeTriangleNormal(scene.meshes[itsindex].tris[moreindex],scene);
    }
    tolight.a = point + NORMAL*scene.shadow_ray_epsilon;
    for(int l=0; l<scene.point_lights.size(); l++){
        tolight.b = scene.point_lights[l].position - point;
        tolight.a = point + NORMAL*scene.shadow_ray_epsilon;
        t_minL = length(tolight.b);
        tolight.b = normalize(tolight.b);
        intersect = 0;
        // WILL BE CHANGED...
        rayHitsWhat(tolight, intersect, dummy, dummy2, t_minL, false, BETA, GAMMA, scene);
    
        if(intersect==0){
            if(which==1){
                Vec3f norm = computeTriangleNormal(scene.triangles[itsindex],scene);
                double dist = distance(tolight.a,scene.point_lights[l].position);
                if(mapping){
                    if(scene.textures[scene.triangles[itsindex].texture_id-1].decalMode=="replace_kd"){
                        dff1 = (THEKD.x/255) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.x / dist / dist;
                        dff2 = (THEKD.y/255) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.y / dist / dist;
                        dff3 = (THEKD.z/255) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.z / dist / dist;
                    }
                    else if (scene.textures[scene.triangles[itsindex].texture_id-1].decalMode=="blend_kd"){
                        dff1 = 0.5 * (scene.materials[scene.triangles[itsindex].material_id-1].diffuse.x + (THEKD.x/255)) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.x / dist / dist;
                        dff2 = 0.5 * (scene.materials[scene.triangles[itsindex].material_id-1].diffuse.y + (THEKD.y/255)) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.y / dist / dist;
                        dff3 = 0.5 * (scene.materials[scene.triangles[itsindex].material_id-1].diffuse.z + (THEKD.z/255)) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.z / dist / dist;
                    }
                    else if (scene.textures[scene.triangles[itsindex].texture_id-1].decalMode=="replace_all"){
                        dff1 = THEKD.x;
                        dff2 = THEKD.y;
                        dff3 = THEKD.z;
                    }
                }
                else{
                    dff1 = scene.materials[scene.triangles[itsindex].material_id-1].diffuse.x * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.x / dist / dist;
                    dff2 = scene.materials[scene.triangles[itsindex].material_id-1].diffuse.y * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.y / dist / dist;
                    dff3 = scene.materials[scene.triangles[itsindex].material_id-1].diffuse.z * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.z / dist / dist;
                }
                Vec3f toCam = normalize(from - tolight.a);
                Vec3f H = (tolight.b + toCam)/length(tolight.b + toCam);
                double exp = scene.materials[scene.triangles[itsindex].material_id-1].phong_exponent;
                bp1 = pow(fmax(0,DOT(norm,H)),exp) * scene.materials[scene.triangles[itsindex].material_id-1].specular.x * scene.point_lights[l].intensity.x / dist / dist;
                bp2 = pow(fmax(0,DOT(norm,H)),exp) * scene.materials[scene.triangles[itsindex].material_id-1].specular.y * scene.point_lights[l].intensity.y / dist / dist;
                bp3 = pow(fmax(0,DOT(norm,H)),exp) * scene.materials[scene.triangles[itsindex].material_id-1].specular.z * scene.point_lights[l].intensity.z / dist / dist;
                color[0] = (color[0] + dff1 + bp1) > 255 ? 255 : (color[0] + dff1 + bp1);
                color[1] = (color[1] + dff2 + bp2) > 255 ? 255 : (color[1] + dff2 + bp2);
                color[2] = (color[2] + dff3 + bp3) > 255 ? 255 : (color[2] + dff3 + bp3);
            }
            else if(which==2){
                Vec3f norm = computeSphereNormal(scene.spheres[itsindex],tolight.a,scene);
                double dist = distance(tolight.a,scene.point_lights[l].position);
                if(mapping){
                    if(scene.textures[scene.spheres[itsindex].texture_id-1].decalMode=="replace_kd"){
                        dff1 = (THEKD.x/255) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.x / dist / dist;
                        dff2 = (THEKD.y/255) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.y / dist / dist;
                        dff3 = (THEKD.z/255) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.z / dist / dist;
                    }
                    else if (scene.textures[scene.spheres[itsindex].texture_id-1].decalMode=="blend_kd"){
                        dff1 = 0.5 * (scene.materials[scene.spheres[itsindex].material_id-1].diffuse.x + (THEKD.x/255)) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.x / dist / dist;
                        dff2 = 0.5 * (scene.materials[scene.spheres[itsindex].material_id-1].diffuse.y + (THEKD.y/255)) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.y / dist / dist;
                        dff3 = 0.5 * (scene.materials[scene.spheres[itsindex].material_id-1].diffuse.z + (THEKD.z/255)) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.z / dist / dist;
                    }
                    else if (scene.textures[scene.spheres[itsindex].texture_id-1].decalMode=="replace_all"){
                        dff1 = THEKD.x;
                        dff2 = THEKD.y;
                        dff3 = THEKD.z;
                    }
                }
                else{
                    dff1 = scene.materials[scene.spheres[itsindex].material_id-1].diffuse.x * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.x / dist / dist;
                    dff2 = scene.materials[scene.spheres[itsindex].material_id-1].diffuse.y * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.y / dist / dist;
                    dff3 = scene.materials[scene.spheres[itsindex].material_id-1].diffuse.z * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.z / dist / dist;
                }
                Vec3f toCam = normalize(from - tolight.a);
                Vec3f H = (tolight.b + toCam)/length(tolight.b + toCam);
                double exp = scene.materials[scene.spheres[itsindex].material_id-1].phong_exponent;
                bp1 = pow(fmax(0,DOT(norm,H)),exp) * scene.materials[scene.spheres[itsindex].material_id-1].specular.x * scene.point_lights[l].intensity.x / dist / dist;
                bp2 = pow(fmax(0,DOT(norm,H)),exp) * scene.materials[scene.spheres[itsindex].material_id-1].specular.y * scene.point_lights[l].intensity.y / dist / dist;
                bp3 = pow(fmax(0,DOT(norm,H)),exp) * scene.materials[scene.spheres[itsindex].material_id-1].specular.z * scene.point_lights[l].intensity.z / dist / dist;
                color[0] = (color[0] + dff1 + bp1) > 255 ? 255 : (color[0] + dff1 + bp1);
                color[1] = (color[1] + dff2 + bp2) > 255 ? 255 : (color[1] + dff2 + bp2);
                color[2] = (color[2] + dff3 + bp3) > 255 ? 255 : (color[2] + dff3 + bp3);
            }
            else if(which==3){
                Vec3f norm = computeTriangleNormal(scene.meshes[itsindex].tris[moreindex],scene);
                double dist = distance(tolight.a,scene.point_lights[l].position);
                if(mapping){
                    if(scene.textures[scene.meshes[itsindex].texture_id-1].decalMode=="replace_kd"){
                        dff1 = (THEKD.x/255) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.x / dist / dist;
                        dff2 = (THEKD.y/255) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.y / dist / dist;
                        dff3 = (THEKD.z/255) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.z / dist / dist;
                    }
                    else if (scene.textures[scene.meshes[itsindex].texture_id-1].decalMode=="blend_kd"){
                        dff1 = 0.5 * (scene.materials[scene.meshes[itsindex].material_id-1].diffuse.x + (THEKD.x/255)) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.x / dist / dist;
                        dff2 = 0.5 * (scene.materials[scene.meshes[itsindex].material_id-1].diffuse.y + (THEKD.y/255)) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.y / dist / dist;
                        dff3 = 0.5 * (scene.materials[scene.meshes[itsindex].material_id-1].diffuse.z + (THEKD.z/255)) * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.z / dist / dist;
                    }
                    else if (scene.textures[scene.meshes[itsindex].texture_id-1].decalMode=="replace_all"){
                        dff1 = THEKD.x;
                        dff2 = THEKD.y;
                        dff3 = THEKD.z;
                    }
                }
                else{
                    dff1 = scene.materials[scene.meshes[itsindex].tris[moreindex].material_id-1].diffuse.x * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.x / dist / dist;
                    dff2 = scene.materials[scene.meshes[itsindex].tris[moreindex].material_id-1].diffuse.y * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.y / dist / dist;
                    dff3 = scene.materials[scene.meshes[itsindex].tris[moreindex].material_id-1].diffuse.z * fmax(0,DOT(tolight.b,norm)) * scene.point_lights[l].intensity.z / dist / dist;
                }
                Vec3f toCam = normalize(from - tolight.a);
                Vec3f H = (tolight.b + toCam)/length(tolight.b + toCam);
                double exp = scene.materials[scene.meshes[itsindex].tris[moreindex].material_id-1].phong_exponent;
                bp1 = pow(fmax(0,DOT(norm,H)),exp) * scene.materials[scene.meshes[itsindex].tris[moreindex].material_id-1].specular.x * scene.point_lights[l].intensity.x / dist / dist;
                bp2 = pow(fmax(0,DOT(norm,H)),exp) * scene.materials[scene.meshes[itsindex].tris[moreindex].material_id-1].specular.y * scene.point_lights[l].intensity.y / dist / dist;
                bp3 = pow(fmax(0,DOT(norm,H)),exp) * scene.materials[scene.meshes[itsindex].tris[moreindex].material_id-1].specular.z * scene.point_lights[l].intensity.z / dist / dist;
                color[0] = (color[0] + dff1 + bp1) > 255 ? 255 : (color[0] + dff1 + bp1);
                color[1] = (color[1] + dff2 + bp2) > 255 ? 255 : (color[1] + dff2 + bp2);
                color[2] = (color[2] + dff3 + bp3) > 255 ? 255 : (color[2] + dff3 + bp3);
            }
        }
    }
    return;
}

void reflect(Ray r, int*& color, int depth, int initial_depth, double& BETA, double& GAMMA,Scene& sc){
    double tmin=DBL_MAX;
    int closest = -1;
    int moreclosest = -1;
    int which = 0;
    if(depth==initial_depth)
        rayHitsWhat(r,which,closest,moreclosest,tmin,true,BETA, GAMMA,sc);
    else
        rayHitsWhat(r,which,closest,moreclosest,tmin,false,BETA, GAMMA,sc);
    if(which==0){
        color[0] = sc.background_color.x;
        color[1] = sc.background_color.y;
        color[2] = sc.background_color.z;
        return;
    }
    Vec3f intersection_point;
    if(depth==initial_depth)
        intersection_point = r.a + (r.b-r.a)*tmin;
    else
        intersection_point = r.a + r.b*tmin;
    if(depth==0){
        calculateColorOfPoint(r.a,intersection_point, color, which, closest, moreclosest, BETA, GAMMA,sc);
        return;
    }
    Vec3f WO = normalize(r.a - intersection_point);
    Vec3f N;
    if(which==1)
        N = computeTriangleNormal(sc.triangles[closest], sc);
    else if(which==2)
        N = computeSphereNormal(sc.spheres[closest], intersection_point, sc);
    else if(which==3)
        N = computeTriangleNormal(sc.meshes[closest].tris[moreclosest], sc);
    Vec3f WR = normalize(-WO + N*(DOT(N,WO))*2);
    Ray next;
    next.a = intersection_point + N*sc.shadow_ray_epsilon;
    next.b = WR;
    //RECURSIVE PART
    reflect(next,color,depth-1,initial_depth,BETA, GAMMA,sc);
    //RECURSIVE PART
    int TEMP[3];
    int* tempcolor = TEMP;
    tempcolor[0] = tempcolor[1] = tempcolor[2] = 0;
    if(depth!=initial_depth){
        calculateColorOfPoint(r.a,intersection_point, tempcolor, which, closest, moreclosest,BETA, GAMMA, sc);
        if(which==1){
            color[0] = color[0] * sc.materials[sc.triangles[closest].material_id-1].mirror.x + tempcolor[0] > 255 ? 255 : color[0] * sc.materials[sc.triangles[closest].material_id-1].mirror.x + tempcolor[0];
            color[1] = color[1] * sc.materials[sc.triangles[closest].material_id-1].mirror.y + tempcolor[1] > 255 ? 255 : color[1] * sc.materials[sc.triangles[closest].material_id-1].mirror.y + tempcolor[1];
            color[2] = color[2] * sc.materials[sc.triangles[closest].material_id-1].mirror.z + tempcolor[2] > 255 ? 255 : color[2] * sc.materials[sc.triangles[closest].material_id-1].mirror.z + tempcolor[2]; 
        }                                 
        else if(which==2){
            color[0] = color[0] * sc.materials[sc.spheres[closest].material_id-1].mirror.x + tempcolor[0] > 255 ? 255 :  color[0] * sc.materials[sc.spheres[closest].material_id-1].mirror.x + tempcolor[0];
            color[1] = color[1] * sc.materials[sc.spheres[closest].material_id-1].mirror.y + tempcolor[1] > 255 ? 255 :  color[1] * sc.materials[sc.spheres[closest].material_id-1].mirror.y + tempcolor[1];
            color[2] = color[2] * sc.materials[sc.spheres[closest].material_id-1].mirror.z + tempcolor[2] > 255 ? 255 :  color[2] * sc.materials[sc.spheres[closest].material_id-1].mirror.z + tempcolor[2]; 
        }
        else{
            color[0] = color[0] * sc.materials[sc.meshes[closest].tris[moreclosest].material_id-1].mirror.x + tempcolor[0] > 255 ? 255 : color[0] * sc.materials[sc.meshes[closest].tris[moreclosest].material_id-1].mirror.x + tempcolor[0];
            color[1] = color[1] * sc.materials[sc.meshes[closest].tris[moreclosest].material_id-1].mirror.y + tempcolor[1] > 255 ? 255 : color[1] * sc.materials[sc.meshes[closest].tris[moreclosest].material_id-1].mirror.y + tempcolor[1];
            color[2] = color[2] * sc.materials[sc.meshes[closest].tris[moreclosest].material_id-1].mirror.z + tempcolor[2] > 255 ? 255 : color[2] * sc.materials[sc.meshes[closest].tris[moreclosest].material_id-1].mirror.z + tempcolor[2]; 
        }   
    }
    else{
        if(which==1){
            color[0] = color[0] * sc.materials[sc.triangles[closest].material_id-1].mirror.x > 255 ? 255 : color[0] * sc.materials[sc.triangles[closest].material_id-1].mirror.x;
            color[1] = color[1] * sc.materials[sc.triangles[closest].material_id-1].mirror.y > 255 ? 255 : color[1] * sc.materials[sc.triangles[closest].material_id-1].mirror.y;
            color[2] = color[2] * sc.materials[sc.triangles[closest].material_id-1].mirror.z > 255 ? 255 : color[2] * sc.materials[sc.triangles[closest].material_id-1].mirror.z; 
        }
        else if(which==2){
            color[0] = color[0] * sc.materials[sc.spheres[closest].material_id-1].mirror.x > 255 ? 255 :  color[0] * sc.materials[sc.spheres[closest].material_id-1].mirror.x;
            color[1] = color[1] * sc.materials[sc.spheres[closest].material_id-1].mirror.y > 255 ? 255 :  color[1] * sc.materials[sc.spheres[closest].material_id-1].mirror.y;
            color[2] = color[2] * sc.materials[sc.spheres[closest].material_id-1].mirror.z > 255 ? 255 :  color[2] * sc.materials[sc.spheres[closest].material_id-1].mirror.z; 
        }
        else{
            color[0] = color[0] * sc.materials[sc.meshes[closest].tris[moreclosest].material_id-1].mirror.x > 255 ? 255 : color[0] * sc.materials[sc.meshes[closest].tris[moreclosest].material_id-1].mirror.x;
            color[1] = color[1] * sc.materials[sc.meshes[closest].tris[moreclosest].material_id-1].mirror.y > 255 ? 255 : color[1] * sc.materials[sc.meshes[closest].tris[moreclosest].material_id-1].mirror.y;
            color[2] = color[2] * sc.materials[sc.meshes[closest].tris[moreclosest].material_id-1].mirror.z > 255 ? 255 : color[2] * sc.materials[sc.meshes[closest].tris[moreclosest].material_id-1].mirror.z; 
        }   
    }
}

double intersectQuad(Ray& r, Vec3f& m1, Vec3f& m2, Vec3f& m3, bool forpixel, Mesh& mesh, Scene& sc)
{
    double  a,b,c,d,e,f,g,h,i,j,k,l;
    double beta,gamma,t;
    double eimhf,gfmdi,dhmeg,akmjb,jcmal,blmkc;
    double M;
    double dd;
    if (forpixel){
        r.b = r.b - r.a;
    }
    Vec3f ma,mb,mc;
    ma = m1;
    mb = m2;
    mc = m3;
    a = ma.x-mb.x;
    b = ma.y-mb.y;
    c = ma.z-mb.z;
    d = ma.x-mc.x;
    e = ma.y-mc.y;
    f = ma.z-mc.z;
    g = r.b.x;
    h = r.b.y;
    i = r.b.z;
    j = ma.x-r.a.x;
    k = ma.y-r.a.y;
    l = ma.z-r.a.z; 
    eimhf = e*i-h*f;
    gfmdi = g*f-d*i;
    dhmeg = d*h-e*g;
    akmjb = a*k-j*b;
    jcmal = j*c-a*l;
    blmkc = b*l-k*c;
    M = a*eimhf+b*gfmdi+c*dhmeg;
    if (forpixel){
        r.b = r.b + r.a;
    }
    if (M==0) return -1;
    t = -(f*akmjb+e*jcmal+d*blmkc)/M;
    gamma = (i*akmjb+h*jcmal+g*blmkc)/M;
    if (gamma<-1 || gamma>1) return -1;
    beta = (j*eimhf+k*gfmdi+l*dhmeg)/M;
    if ( 0 > (1-beta-gamma) || (1-beta-gamma) > 1 || 0 > beta || beta > 1) return -1; 
    return t;
}
int intersectMesh(Ray& r, Mesh& mesh, bool forpixel, Scene& sc){
    Vec3f e0(mesh.max_x, mesh.min_y, mesh.min_z);
    Vec3f e1(mesh.max_x, mesh.max_y, mesh.min_z);
    Vec3f e2(mesh.min_x, mesh.max_y, mesh.min_z);
    Vec3f e3(mesh.min_x, mesh.min_y, mesh.min_z);
    Vec3f e4(mesh.max_x, mesh.min_y, mesh.max_z);
    Vec3f e5(mesh.max_x, mesh.max_y, mesh.max_z);
    Vec3f e6(mesh.min_x, mesh.max_y, mesh.max_z);
    Vec3f e7(mesh.min_x, mesh.min_y, mesh.max_z);
    double t;
    t = intersectQuad(r, e4, e1, e5, forpixel, mesh, sc);
    if(t>0)
        goto success;
    t = intersectQuad(r, e6, e4, e5, forpixel, mesh, sc);
    if(t>0)
        goto success;
    t = intersectQuad(r, e1, e6, e5, forpixel, mesh, sc);
    if(t>0)
        goto success;
    t = intersectQuad(r, e7, e2, e3, forpixel, mesh, sc);
    if(t>0)
        goto success;
    t = intersectQuad(r, e2, e0, e3, forpixel, mesh, sc);
    if(t>0)
        goto success;
    t = intersectQuad(r, e0, e7, e3, forpixel, mesh, sc);
    if(t>0)
        goto success;
    return -1;
    success:
    if(mesh.tris.size()<50)
        return 3;
    bool signalL = false;
    bool signalR = false;
    Vec3f eLEFT0(mesh.max_xLEFT, mesh.min_yLEFT, mesh.min_zLEFT);
    Vec3f eLEFT1(mesh.max_xLEFT, mesh.max_yLEFT, mesh.min_zLEFT);
    Vec3f eLEFT2(mesh.min_xLEFT, mesh.max_yLEFT, mesh.min_zLEFT);
    Vec3f eLEFT3(mesh.min_xLEFT, mesh.min_yLEFT, mesh.min_zLEFT);
    Vec3f eLEFT4(mesh.max_xLEFT, mesh.min_yLEFT, mesh.max_zLEFT);
    Vec3f eLEFT5(mesh.max_xLEFT, mesh.max_yLEFT, mesh.max_zLEFT);
    Vec3f eLEFT6(mesh.min_xLEFT, mesh.max_yLEFT, mesh.max_zLEFT);
    Vec3f eLEFT7(mesh.min_xLEFT, mesh.min_yLEFT, mesh.max_zLEFT);
    t = intersectQuad(r, eLEFT4, eLEFT1, eLEFT5, forpixel, mesh, sc);
    if(t>0){
        signalL=true;
        goto successL;
    }
    t = intersectQuad(r, eLEFT6, eLEFT4, eLEFT5, forpixel, mesh, sc);
    if(t>0){
        signalL=true;
        goto successL;
    }
    t = intersectQuad(r, eLEFT1, eLEFT6, eLEFT5, forpixel, mesh, sc);
    if(t>0){
        signalL=true;
        goto successL;
    }
    t = intersectQuad(r, eLEFT7, eLEFT2, eLEFT3, forpixel, mesh, sc);
    if(t>0){
        signalL=true;
        goto successL;
    }
    t = intersectQuad(r, eLEFT2, eLEFT0, eLEFT3, forpixel, mesh, sc);
    if(t>0){
        signalL=true;
        goto successL;
    }
    t = intersectQuad(r, eLEFT0, eLEFT7, eLEFT3, forpixel, mesh, sc);
    if(t>0){
        signalL=true;
        goto successL;
    }
    successL:
    Vec3f eRIGHT0(mesh.max_xRIGHT, mesh.min_yRIGHT, mesh.min_zRIGHT);
    Vec3f eRIGHT1(mesh.max_xRIGHT, mesh.max_yRIGHT, mesh.min_zRIGHT);
    Vec3f eRIGHT2(mesh.min_xRIGHT, mesh.max_yRIGHT, mesh.min_zRIGHT);
    Vec3f eRIGHT3(mesh.min_xRIGHT, mesh.min_yRIGHT, mesh.min_zRIGHT);
    Vec3f eRIGHT4(mesh.max_xRIGHT, mesh.min_yRIGHT, mesh.max_zRIGHT);
    Vec3f eRIGHT5(mesh.max_xRIGHT, mesh.max_yRIGHT, mesh.max_zRIGHT);
    Vec3f eRIGHT6(mesh.min_xRIGHT, mesh.max_yRIGHT, mesh.max_zRIGHT);
    Vec3f eRIGHT7(mesh.min_xRIGHT, mesh.min_yRIGHT, mesh.max_zRIGHT);
    t = intersectQuad(r, eRIGHT4, eRIGHT1, eRIGHT5, forpixel, mesh, sc);
    if(t>0){
        signalR=true;
        goto successR;
    }
    t = intersectQuad(r, eRIGHT6, eRIGHT4, eRIGHT5, forpixel, mesh, sc);
    if(t>0){
        signalR=true;
        goto successR;
    }
    t = intersectQuad(r, eRIGHT1, eRIGHT6, eRIGHT5, forpixel, mesh, sc);
    if(t>0){
        signalR=true;
        goto successR;
    }
    t = intersectQuad(r, eRIGHT7, eRIGHT2, eRIGHT3, forpixel, mesh, sc);
    if(t>0){
        signalR=true;
        goto successR;
    }
    t = intersectQuad(r, eRIGHT2, eRIGHT0, eRIGHT3, forpixel, mesh, sc);
    if(t>0){
        signalR=true;
        goto successR;
    }
    t = intersectQuad(r, eRIGHT0, eRIGHT7, eRIGHT3, forpixel, mesh, sc);
    if(t>0){
        signalR=true;
        goto successR;
    }
    successR:
    if(signalR && signalL)
        return 3;
    else if(signalR)
        return 2;
    else if(signalL)
        return 1;
    else
        return -1;
}