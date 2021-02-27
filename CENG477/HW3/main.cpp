#include <iostream>
#include <chrono>
#include <cassert>
#include <cmath>
#include <sstream>
#include <cstring>
#include "parser.h"
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GLFW/glfw3.h>

using parser::Vec3f;
//////-------- Global Variables -------/////////

GLuint gpuVertexBuffer;
GLuint gpuNormalBuffer;
GLuint gpuIndexBuffer;

char gRendererInfo[512] = { 0 };
char gWindowTitle[512] = { 0 };

// Sample usage for reading an XML scene file
parser::Scene scene;
static GLFWwindow* win = NULL;

static void errorCallback(int error, const char* description) {
    fprintf(stderr, "Error: %s\n", description);
}

static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

void setCamera(){
    glViewport(0,0,scene.camera.image_width,scene.camera.image_height);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    float X = scene.camera.gaze.x / sqrt(scene.camera.gaze.x * scene.camera.gaze.x \
                                        +scene.camera.gaze.y * scene.camera.gaze.y \
                                        +scene.camera.gaze.z * scene.camera.gaze.z );
    float Y = scene.camera.gaze.y / sqrt(scene.camera.gaze.x * scene.camera.gaze.x \
                                        +scene.camera.gaze.y * scene.camera.gaze.y \
                                        +scene.camera.gaze.z * scene.camera.gaze.z );
    float Z = scene.camera.gaze.z / sqrt(scene.camera.gaze.x * scene.camera.gaze.x \
                                        +scene.camera.gaze.y * scene.camera.gaze.y \
                                        +scene.camera.gaze.z * scene.camera.gaze.z );

    gluLookAt(scene.camera.position.x,scene.camera.position.y,scene.camera.position.z \
             ,scene.camera.position.x + scene.camera.near_distance*X \
             ,scene.camera.position.y + scene.camera.near_distance*Y \
             ,scene.camera.position.z + scene.camera.near_distance*Z \
             ,scene.camera.up.x,scene.camera.up.y,scene.camera.up.z);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(scene.camera.near_plane.x,scene.camera.near_plane.y,scene.camera.near_plane.z,scene.camera.near_plane.w,scene.camera.near_distance,scene.camera.far_distance);
}

void lightsON(){
    glEnable(GL_LIGHTING);
    int len = scene.point_lights.size();
    GLfloat ambient[] = {scene.ambient_light.x,scene.ambient_light.y,scene.ambient_light.z,1.0f};
    for(int i=0; i<len; i++){
        GLfloat color[] = {scene.point_lights[i].intensity.x, scene.point_lights[i].intensity.y, scene.point_lights[i].intensity.z, 1.0f};
        GLfloat position[] = {scene.point_lights[i].position.x, scene.point_lights[i].position.y, scene.point_lights[i].position.z, 1.0f};
        glLightfv(GL_LIGHT0+i, GL_POSITION, position);
        glLightfv(GL_LIGHT0+i, GL_AMBIENT, ambient);
        glLightfv(GL_LIGHT0+i, GL_DIFFUSE, color);
        glLightfv(GL_LIGHT0+i, GL_SPECULAR, color);
        glEnable(GL_LIGHT0+i);   
    }
}

void lightsOFF(){
    glDisable(GL_LIGHTING);
    int len = scene.point_lights.size();
    for(int i=0; i<len; i++)
        glDisable(GL_LIGHT0 + i);
}

void drawMeshes(){
    static bool firstTime = true;
    static unsigned mesh_size = 0U;
    static unsigned vertex_size = 0U;
    static unsigned face_size = 0U;

    if (firstTime) {
        firstTime = false;
        mesh_size = scene.meshes.size();
        vertex_size = scene.vertex_data.size();
        face_size = 0;
        for(int i=0; i<mesh_size; i++)
            face_size += scene.meshes[i].faces.size();

        // Get vertex positions
        GLfloat vertexData[vertex_size*3];
        for(int i=0; i<vertex_size; i++){
            vertexData[3*i  ] = scene.vertex_data[i].x;
            vertexData[3*i+1] = scene.vertex_data[i].y;
            vertexData[3*i+2] = scene.vertex_data[i].z;
        }

        // Get all indices information
        GLuint indexData[face_size*3];
        unsigned ind = 0U;
        for(int m=0; m<mesh_size; m++){
            unsigned len = scene.meshes[m].faces.size();
            for(int i=0; i<len; i++){
                indexData[3*ind  ] = scene.meshes[m].faces[i].v0_id-1;
                indexData[3*ind+1] = scene.meshes[m].faces[i].v1_id-1;
                indexData[3*ind+2] = scene.meshes[m].faces[i].v2_id-1;
                ind++;
            }
        }

        // Calculate all normals
        GLfloat normalData[vertex_size*3] = {0};
        GLuint normalIncidents[vertex_size] = {0};
        for(int i=0; i<face_size; i++){
            
            parser::Vec3f a,b;
            b.x = vertexData[3*indexData[3*i  ]  ] \
                - vertexData[3*indexData[3*i+1]  ];
            a.x = vertexData[3*indexData[3*i+2]  ] \
                - vertexData[3*indexData[3*i+1]  ];

            b.y = vertexData[3*indexData[3*i  ]+1] \
                - vertexData[3*indexData[3*i+1]+1];
            a.y = vertexData[3*indexData[3*i+2]+1] \
                - vertexData[3*indexData[3*i+1]+1];

            b.z = vertexData[3*indexData[3*i  ]+2] \
                - vertexData[3*indexData[3*i+1]+2];
            a.z = vertexData[3*indexData[3*i+2]+2] \
                - vertexData[3*indexData[3*i+1]+2];

            parser::Vec3f c; // CROSS PRODUCT (V2-V1) x (V0-V1)
            c.x = a.y * b.z - a.z * b.y; c.y = a.z * b.x - a.x * b.z; c.z = a.x * b.y - a.y * b.x;
            
            
            float divider = sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
            c.x /= divider;  // Normalize the normal.
            c.y /= divider;
            c.z /= divider;
        
            normalData[3*indexData[3*i  ]  ] += c.x;
            normalData[3*indexData[3*i  ]+1] += c.y;
            normalData[3*indexData[3*i  ]+2] += c.z;
            normalIncidents[indexData[3*i  ]]++;

            normalData[3*indexData[3*i+1]  ] += c.x;
            normalData[3*indexData[3*i+1]+1] += c.y;
            normalData[3*indexData[3*i+1]+2] += c.z;
            normalIncidents[indexData[3*i+1]]++;

            normalData[3*indexData[3*i+2]  ] += c.x;
            normalData[3*indexData[3*i+2]+1] += c.y;
            normalData[3*indexData[3*i+2]+2] += c.z;
            normalIncidents[indexData[3*i+2]]++;
        }

        for(int i=0; i<vertex_size; i++){
            // normal len / how many face indicent to it.
            normalData[3*i + 0] /= normalIncidents[i];
            normalData[3*i + 1] /= normalIncidents[i];
            normalData[3*i + 2] /= normalIncidents[i];
        }

        GLuint vertexBuffer, indexBuffer;
        glGenBuffers(1, &vertexBuffer);
        glGenBuffers(1, &indexBuffer);
        glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
        glBufferData(GL_ARRAY_BUFFER, 6 * vertex_size * sizeof(GLfloat), 0, GL_STATIC_DRAW);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * face_size * sizeof(GLuint), indexData, GL_STATIC_DRAW);
        glBufferSubData(GL_ARRAY_BUFFER, 3 * vertex_size * sizeof(GLfloat), 3 * vertex_size * sizeof(GLfloat), normalData);
        glBufferSubData(GL_ARRAY_BUFFER, 0, 3 * vertex_size * sizeof(GLfloat), vertexData);
    }

    unsigned current_mesh_index = 0U;
    glVertexPointer(3, GL_FLOAT, 0, 0);
    glNormalPointer(GL_FLOAT, 0, (const void*)(vertex_size * 3 * sizeof(GLfloat)));
    for(int i=0; i<mesh_size; i++){
        glPushMatrix();

        if(scene.meshes[i].mesh_type=="Solid")
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        else
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        for(int j=scene.meshes[i].transformations.size()-1; j>=0; j--){

            std::string type = scene.meshes[i].transformations[j].transformation_type;
            int id = scene.meshes[i].transformations[j].id;
            
            if(type == "Scaling")
                glScalef(scene.scalings[id-1].x, scene.scalings[id-1].y, scene.scalings[id-1].z);
            
            else if(type == "Translation")
                glTranslatef(scene.translations[id-1].x, scene.translations[id-1].y, scene.translations[id-1].z);
            
            else if(type == "Rotation")
                glRotatef(scene.rotations[id-1].x, scene.rotations[id-1].y, scene.rotations[id-1].z, scene.rotations[id-1].w);
            
        }
        
        GLfloat ambient[] = {scene.materials[scene.meshes[i].material_id-1].ambient.x, scene.materials[scene.meshes[i].material_id-1].ambient.y, scene.materials[scene.meshes[i].material_id-1].ambient.z, 1.0f};
        GLfloat diffuse[] = {scene.materials[scene.meshes[i].material_id-1].diffuse.x, scene.materials[scene.meshes[i].material_id-1].diffuse.y, scene.materials[scene.meshes[i].material_id-1].diffuse.z, 1.0f};
        GLfloat specular[] = {scene.materials[scene.meshes[i].material_id-1].specular.x, scene.materials[scene.meshes[i].material_id-1].specular.y, scene.materials[scene.meshes[i].material_id-1].specular.z, 1.0f};
        GLfloat phong_exponent[] = {scene.materials[scene.meshes[i].material_id-1].phong_exponent};

        glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
        glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
        glMaterialfv(GL_FRONT, GL_SHININESS, phong_exponent);
        
        glDrawElements(GL_TRIANGLES, 3 * scene.meshes[i].faces.size(), GL_UNSIGNED_INT, (const void*)(3 * current_mesh_index * sizeof(GLuint)));
        glPopMatrix();
        current_mesh_index += scene.meshes[i].faces.size();
    }
}

int main(int argc, char* argv[]) {
    scene.loadFromXml(argv[1]);

    glfwSetErrorCallback(errorCallback);

    if (!glfwInit()) {
        exit(EXIT_FAILURE);
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);

    win = glfwCreateWindow(scene.camera.image_width, scene.camera.image_height, "CENG477 - HW3", NULL, NULL);
    if (!win) {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }
    glfwMakeContextCurrent(win);

    GLenum err = glewInit();
    if (err != GLEW_OK) {
        fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    glfwSetKeyCallback(win, keyCallback);

    glEnable(GL_DEPTH_TEST);

    if(scene.culling_enabled){
        glEnable(GL_CULL_FACE);
        if(scene.culling_face)
            glCullFace(GL_FRONT);
        else
            glCullFace(GL_BACK);   
    }

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_NORMALIZE);
    lightsON();
    setCamera();
    
    while(!glfwWindowShouldClose(win)){
        static int framesRendered = 0;
        static std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

        glClearColor(scene.background_color.x, scene.background_color.y, scene.background_color.z, 1);
        glClearDepth(1.0f);
        glClearStencil(0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
        glMatrixMode(GL_MODELVIEW);
        
        drawMeshes();
        
        glfwSwapBuffers(win);

        ++framesRendered;

        std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();

        std::chrono::duration<double> elapsedTime = end - start;
        if (elapsedTime.count() > 1.)
        {
            start = std::chrono::system_clock::now();

            std::stringstream stream;

            stream << framesRendered/elapsedTime.count();
            framesRendered = 0;

            strcpy(gWindowTitle, gRendererInfo);
            strcat(gWindowTitle, "CENG477 - HW3 [");
            strcat(gWindowTitle, stream.str().c_str());
            strcat(gWindowTitle, " FPS]");
        }
            glfwSetWindowTitle(win, gWindowTitle);
            glfwWaitEvents();
    }

    glfwDestroyWindow(win);
    glfwTerminate();

    exit(EXIT_SUCCESS);

    return 0;
}

