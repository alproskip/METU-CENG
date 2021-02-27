#version 330

layout(location = 0) in vec3 position;

// Data from CPU 
uniform mat4 MVP; // ModelViewProjection Matrix
uniform mat4 MV; // ModelView idMVPMatrix
uniform vec3 cameraPosition;
uniform vec3 lightPosition;
uniform float heightFactor;

// Texture-related data
uniform sampler2D heightMap;
uniform int widthTexture;
uniform int heightTexture;
uniform int move;


// Output to Fragment Shader
out vec2 textureCoordinate; // For texture-color
out vec3 vertexNormal; // For Lighting computation
out vec3 ToLightVector; // Vector from Vertex to Light;
out vec3 ToCameraVector; // Vector from Vertex to Camera;


void main()
{
    textureCoordinate = vec2(float(-position.x)/float(widthTexture+1), float(-position.z)/float(heightTexture+1));
    vec4 textureColor = texture(heightMap, textureCoordinate);

    vec4 h1 = texture(heightMap, vec2(float(-(position.x+1))/float(widthTexture) ,float(-(position.z-1))/float(heightTexture)));
    vec4 h2 = texture(heightMap, vec2(float(-(position.x  ))/float(widthTexture) ,float(-(position.z-1))/float(heightTexture)));
    vec4 h3 = texture(heightMap, vec2(float(-(position.x-1))/float(widthTexture) ,float(-(position.z  ))/float(heightTexture)));
    vec4 h4 = texture(heightMap, vec2(float(-(position.x-1))/float(widthTexture) ,float(-(position.z+1))/float(heightTexture)));
    vec4 h5 = texture(heightMap, vec2(float(-(position.x  ))/float(widthTexture) ,float(-(position.z+1))/float(heightTexture)));
    vec4 h6 = texture(heightMap, vec2(float(-(position.x+1))/float(widthTexture) ,float(-(position.z  ))/float(heightTexture)));
    
    vec3 ray1 = vec3(1 , (h1.r - textureColor.r) * heightFactor, -1);
    vec3 ray2 = vec3(0 , (h2.r - textureColor.r) * heightFactor, -1);
    vec3 ray3 = vec3(-1, (h3.r - textureColor.r) * heightFactor, 0 );
    vec3 ray4 = vec3(-1, (h4.r - textureColor.r) * heightFactor, 1 );
    vec3 ray5 = vec3(0 , (h5.r - textureColor.r) * heightFactor, 1 );
    vec3 ray6 = vec3(1 , (h6.r - textureColor.r) * heightFactor, 0 );
    vec3 thenormal = vec3(0.0,0.0,0.0);

    thenormal += cross(ray1, ray2);
    thenormal += cross(ray2, ray3);
    thenormal += cross(ray3, ray4);
    thenormal += cross(ray4, ray5);
    thenormal += cross(ray5, ray6);
    thenormal += cross(ray6, ray1);

    vertexNormal = normalize(thenormal);
    ToCameraVector = normalize(cameraPosition - vec3(position.x, textureColor.r * heightFactor, position.z));
    ToLightVector = normalize(lightPosition - vec3(position.x, textureColor.r * heightFactor, position.z));
    gl_Position = MVP * vec4(position.x-move, textureColor.r * heightFactor, position.z, 1);
}
