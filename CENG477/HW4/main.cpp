#include "helper.h"
#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"
#include "glm/gtx/rotate_vector.hpp"
#include "glm/gtc/type_ptr.hpp"

static GLFWwindow* win = NULL;
int widthWindow = 1000, heightWindow = 1000;

// Shaders
GLuint idProgramShader;
GLuint idFragmentShader;
GLuint idVertexShader;
GLuint idJpegTexture;
GLuint idHeightTexture;
GLuint idMVPMatrix;

// Buffers
GLuint idVertexBuffer;
GLuint idAttributeBuffer;

int textureWidth, textureHeight;
float heightFactor;
glm::vec3 lightPos;
GLuint depthMapFBO;
GLuint depthCubemap;
bool lightPosFlag = false;

// VARIABLE WE HAVE ADDED
std::string vs = "shader.vert";
std::string fs = "shader.frag";
glm::vec3* vertices;
GLuint vertices_size;
glm::vec3 cameraPos;
glm::vec3 up;
glm::vec3 gaze;
glm::vec3 LEFT;
glm::mat4 P;
glm::mat4 MV;
glm::mat4 MVP;
glm::vec3 speed;
bool firstTime = true;
GLuint MOVE;
bool* keystates = new bool[18];
bool fullscreen = false;

void toggleFullscreen(){
  if (fullscreen) {
    glfwSetWindowMonitor(win, nullptr, widthWindow, heightWindow, widthWindow, heightWindow, 0);
    glViewport(0, 0, widthWindow, heightWindow);
  }
  else {
    glfwGetWindowPos(win, &widthWindow, &heightWindow);
    glfwGetWindowSize(win, &widthWindow, &heightWindow);
    int newW = glfwGetVideoMode(glfwGetPrimaryMonitor())->width;
    int newH = glfwGetVideoMode(glfwGetPrimaryMonitor())->height;
    glfwSetWindowMonitor(win, glfwGetPrimaryMonitor(), 0, 0, newW, newH, 0);
    glViewport(0, 0, newW, newH);
  }
  fullscreen = !fullscreen;
}

void creatTheFlatEarth(){
  vertices_size = 6 * textureHeight * textureWidth;
  vertices = new glm::vec3[vertices_size];
  unsigned k = 0;
  for(int i=0; i<textureHeight; i++){
    for(int j=0; j<textureWidth; j++){
      vertices[k++] = glm::vec3(j  ,0.0,i  );
      vertices[k++] = glm::vec3(j  ,0.0,i+1);
      vertices[k++] = glm::vec3(j+1,0.0,i  );
      vertices[k++] = glm::vec3(j+1,0.0,i  );
      vertices[k++] = glm::vec3(j  ,0.0,i+1);
      vertices[k++] = glm::vec3(j+1,0.0,i+1);
    }
  }
  glGenBuffers(1, &idVertexBuffer);
  glGenVertexArrays(1, &idAttributeBuffer);
  glBindBuffer(GL_ARRAY_BUFFER, idVertexBuffer);
  glBindVertexArray(idAttributeBuffer);
  glBufferData(GL_ARRAY_BUFFER, vertices_size * sizeof(glm::vec3), vertices, GL_STATIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
  glEnableVertexAttribArray(0);
  delete [] vertices;
}

void setScene(){

  if(firstTime){
    for(int i=0; i<18; i++)
      keystates[i] = false;
    firstTime = false;
    up = glm::vec3(0.0, 1.0, 0.0);
    gaze = glm::vec3(0.0, 0.0, 1.0);
    cameraPos = glm::vec3(textureWidth/2.0, textureWidth/10.0, -textureWidth/4.0);
    lightPos = glm::vec3(textureWidth/2.0, 100.0, textureHeight/2.0);
    P = glm::perspective(45.0f, 1.0f, 0.1f, 1000.0f);
    MV = glm::lookAt(cameraPos, cameraPos + gaze, up);
    MVP = P * MV * glm::mat4(1.0);
    heightFactor = 10;
    speed = glm::vec3(0.0,0.0,0.0);
    MOVE = 0;
  }
  else{
    cameraPos += speed*gaze;
    MV = glm::mat4(glm::lookAt(cameraPos, cameraPos + gaze, up));
    MVP = P * MV * glm::mat4(1.0);
  }
  
  glUniformMatrix4fv(glGetUniformLocation(idProgramShader, "MVP"), 1, GL_FALSE, glm::value_ptr(MVP));
  glUniformMatrix4fv(glGetUniformLocation(idProgramShader, "MV"), 1, GL_FALSE, glm::value_ptr(MV));
  glUniform1f(glGetUniformLocation(idProgramShader, "heightFactor"), heightFactor);
  glUniform1i(glGetUniformLocation(idProgramShader, "widthTexture"), textureWidth);
  glUniform1i(glGetUniformLocation(idProgramShader, "heightTexture"), textureHeight);
  glUniform3fv(glGetUniformLocation(idProgramShader, "cameraPosition"), 1, glm::value_ptr(cameraPos));
  glUniform3fv(glGetUniformLocation(idProgramShader, "lightPosition"), 1, glm::value_ptr(lightPos));
  glUniform1i(glGetUniformLocation(idProgramShader, "rgbTexture"), 0);
  glUniform1i(glGetUniformLocation(idProgramShader, "heightMap"), 1);
  glUniform1i(glGetUniformLocation(idProgramShader, "move"), MOVE);

}

static void errorCallback(int error, const char* description)
{
  fprintf(stderr, "Error: %s\n", description);
}


static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) glfwSetWindowShouldClose(window, GLFW_TRUE);
    else if(key == GLFW_KEY_R && action == GLFW_PRESS) keystates[0] = true;
    else if(key == GLFW_KEY_R && action == GLFW_RELEASE) keystates[0] = false;
    else if(key == GLFW_KEY_F && action == GLFW_PRESS) keystates[1] = true;
    else if(key == GLFW_KEY_F && action == GLFW_RELEASE) keystates[1] = false;
    else if(key == GLFW_KEY_T && action == GLFW_PRESS) keystates[2] = true;
    else if(key == GLFW_KEY_T && action == GLFW_RELEASE) keystates[2] = false;
    else if(key == GLFW_KEY_G && action == GLFW_PRESS) keystates[3] = true;
    else if(key == GLFW_KEY_G && action == GLFW_RELEASE) keystates[3] = false;
    else if(key == GLFW_KEY_LEFT && action == GLFW_PRESS) keystates[4] = true;
    else if(key == GLFW_KEY_LEFT && action == GLFW_RELEASE) keystates[4] = false;
    else if(key == GLFW_KEY_RIGHT && action == GLFW_PRESS) keystates[5] = true;
    else if(key == GLFW_KEY_RIGHT && action == GLFW_RELEASE) keystates[5] = false;
    else if(key == GLFW_KEY_DOWN && action == GLFW_PRESS) keystates[6] = true;
    else if(key == GLFW_KEY_DOWN && action == GLFW_RELEASE) keystates[6] = false;
    else if(key == GLFW_KEY_UP && action == GLFW_PRESS) keystates[7] = true;
    else if(key == GLFW_KEY_UP && action == GLFW_RELEASE) keystates[7] = false;
    else if(key == GLFW_KEY_I && action == GLFW_PRESS) keystates[8] = true;
    else if(key == GLFW_KEY_I && action == GLFW_RELEASE) keystates[8] = false;
    else if(key == GLFW_KEY_Q && action == GLFW_PRESS) keystates[9] = true;
    else if(key == GLFW_KEY_Q && action == GLFW_RELEASE) keystates[9] = false;
    else if(key == GLFW_KEY_E && action == GLFW_PRESS) keystates[10] = true;
    else if(key == GLFW_KEY_E && action == GLFW_RELEASE) keystates[10] = false;
    else if(key == GLFW_KEY_Y && action == GLFW_PRESS) keystates[11] = true;
    else if(key == GLFW_KEY_Y && action == GLFW_RELEASE) keystates[11] = false;
    else if(key == GLFW_KEY_H && action == GLFW_PRESS) keystates[12] = true;
    else if(key == GLFW_KEY_H && action == GLFW_RELEASE) keystates[12] = false;
    else if(key == GLFW_KEY_X && action == GLFW_PRESS) keystates[13] = true;
    else if(key == GLFW_KEY_X && action == GLFW_RELEASE) keystates[13] = false;
    else if(key == GLFW_KEY_S && action == GLFW_PRESS) keystates[14] = true;
    else if(key == GLFW_KEY_S && action == GLFW_RELEASE) keystates[14] = false;
    else if(key == GLFW_KEY_A && action == GLFW_PRESS) keystates[15] = true;
    else if(key == GLFW_KEY_A && action == GLFW_RELEASE) keystates[15] = false;
    else if(key == GLFW_KEY_W && action == GLFW_PRESS) keystates[16] = true;
    else if(key == GLFW_KEY_W && action == GLFW_RELEASE) keystates[16] = false;
    else if(key == GLFW_KEY_D && action == GLFW_PRESS) keystates[17] = true;
    else if(key == GLFW_KEY_D && action == GLFW_RELEASE) keystates[17] = false;
    else if (key == GLFW_KEY_P && action == GLFW_PRESS) toggleFullscreen();
    else ;

    if(keystates[0])
      heightFactor += 0.5;
    else if(keystates[1])
      heightFactor -= 0.5;
    else if(keystates[2])
      lightPos.y += 5;
    else if(keystates[3])
      lightPos.y -= 5;
    else if(keystates[4])
      lightPos.x += 5;
    else if(keystates[5])
      lightPos.x -= 5;
    else if(keystates[6])
      lightPos.z -= 5;
    else if(keystates[7])
      lightPos.z += 5;
    else if(keystates[8])
	    firstTime = true;
    else if(keystates[9])
      MOVE--;
    else if(keystates[10])
      MOVE++;
    else if(keystates[11])
      speed += 0.25;
    else if(keystates[12])
      speed -= 0.25;
    else if(keystates[13])
      speed *= 0.0; 
    else if(keystates[14]){
      float angle = 0.05;
      LEFT = normalize(cross(up,gaze));
	    gaze = normalize(glm::rotate(gaze,angle,LEFT));
      up = normalize(glm::rotate(up,angle,LEFT));
    }
    else if(keystates[15]){
      float angle = 0.05;
	    gaze = normalize(glm::rotate(gaze,angle,up));
    }
    else if(keystates[16]){
      float angle = -0.05;
      LEFT = normalize(cross(up,gaze));
	    gaze = normalize(glm::rotate(gaze,angle,LEFT));
      up = normalize(glm::rotate(up,angle,LEFT));
    }
    else if(keystates[17]){
      float angle = -0.05;
	    gaze = normalize(glm::rotate(gaze,angle,up));
    }
    else ;
}

int main(int argc, char *argv[]) {

  if (argc != 3) {
    printf("Please provide height and texture image files!\n");
    exit(-1);
  }
  glfwSetErrorCallback(errorCallback);
  if (!glfwInit()) {
    exit(-1);
  }

  //glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  //glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_ANY_PROFILE);
  
  win = glfwCreateWindow(widthWindow, heightWindow, "CENG477 - HW4", NULL, NULL);
  if (!win) {
      glfwTerminate();
      exit(-1);
  }
  glfwMakeContextCurrent(win);
  GLenum err = glewInit();
  if (err != GLEW_OK) {
      fprintf(stderr, "Error: %s\n", glewGetErrorString(err));

      glfwTerminate();
      exit(-1);
  }
  
  initShaders(idProgramShader,vs,fs);
  glUseProgram(idProgramShader);
  glfwSetKeyCallback(win, keyCallback);
  initTexture(argv[1], argv[2], &textureWidth, &textureHeight);
  creatTheFlatEarth();
  setScene();
  glEnable(GL_DEPTH_TEST);
  while(!glfwWindowShouldClose(win)) {
    glClearColor(0,0,0,1);
    glClearDepth(1.0f);
    glClearStencil(0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    setScene();
    glDrawArrays(GL_TRIANGLES,0,18*textureWidth*textureHeight);
    glfwSwapBuffers(win);
    glfwPollEvents();
  }

  glfwDestroyWindow(win);
  glfwTerminate();
  delete [] keystates;
  return 0;
}
