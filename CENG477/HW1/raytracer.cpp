#include "parser.h"
#include "ppm.h"
#include "Math.h"
using namespace parser;


typedef unsigned char RGB[3];

int main(int argc, char* argv[])
{
    Scene scene;
    scene.loadFromXml(argv[1]);
    scene.AllMeshToTri();
    scene.calculateEdges();
    for(int c=0; c<scene.cameras.size(); c++){
        Camera camera = scene.cameras[c];
        Vec3f su,sv,s,u;
        u = normalize(cross(camera.gaze,camera.up));
        int width = camera.image_width, height = camera.image_height;
        Vec4f dim = camera.near_plane;
        double pixelW = (dim.y - dim.x)/(double)width;
        double pixelH = (dim.w - dim.z)/(double)height;
        Vec3f ref = camera.position + camera.gaze * camera.near_distance;
        ref = ref + u * dim.x + camera.up * dim.w;
        scene.refS = ref;
        unsigned char* image = new unsigned char [width * height * 3];
        unsigned char* mirror = new unsigned char [width * height * 3];
        for (int y = 0; y < height; ++y) 
        {
            for (int x = 0; x < width; ++x)
            {
                Ray ray;
                ray.a = camera.position;
                su = u * ((x+0.5)*pixelW);
                sv = camera.up * (-((y+0.5)*pixelH));
                ray.b = ref + su + sv;

                double t_min= DBL_MAX;
                int which = 0;
                int itsindex = -1;
                int moreindex = -1;
                rayHitsWhat(ray, which, itsindex, moreindex, t_min, true, scene);
        
                if(which==0){
                    image[3*(x+width*y)]=scene.background_color.x;
                    image[3*(x+width*y)+1]=scene.background_color.y;
                    image[3*(x+width*y)+2]=scene.background_color.z;
                    continue;
                }

                int TT[3];
                int* clr = TT;
                calculateColorOfPoint(scene.cameras[c].position,(ray.a + (ray.b - ray.a)*(t_min)), clr, which, itsindex, moreindex, scene);
                image[3*(x+width*y)] = clr[0];
                image[3*(x+width*y)+1] = clr[1];
                image[3*(x+width*y)+2] = clr[2];
            }
        }
        
        // MIRROR THINGS
        int mirrortemp[3];
        int* mrr = mirrortemp;
        int rec_dep = scene.max_recursion_depth;

        if(rec_dep){
            for (int y = 0; y < height; ++y)
            {
                for (int x = 0; x < width; ++x)
                {
                    Ray ray;
                    ray.a = camera.position;
                    su = u * ((x+0.5)*pixelW);
                    sv = camera.up * (-((y+0.5)*pixelH));
                    ray.b = ref + su + sv;
                    mrr[0] = 0; mrr[1]=0; mrr[2]=0;
                    reflect(ray, mrr, rec_dep, rec_dep, scene);
                    mirror[3*(x+width*y)] = mrr[0];
                    mirror[3*(x+width*y)+1] = mrr[1];
                    mirror[3*(x+width*y)+2] = mrr[2];
                }
            }
            for (int i = 0; i < height*width*3; i++)
            {
                image[i] = image[i] + mirror[i] > 255 ? 255 : image[i] + mirror[i];
            }
        }
        const char* fNAME = scene.cameras[c].image_name.c_str();
        write_ppm(fNAME, image, width, height);
    }
}

