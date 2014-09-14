uniform mat4 RSMMatrix;

layout(location = 0) in vec3 Position;
layout(location = 1) in vec3 Normal;
layout(location = 2) in vec4 Color;
layout(location = 3) in vec2 Texcoord;
layout(location = 4) in vec2 SecondTexcoord;
layout(location = 7) in vec3 Origin;
layout(location = 8) in vec3 Orientation;
layout(location = 9) in vec3 Scale;
#ifdef GL_ARB_bindless_texture
layout(location = 10) in uvec2 Handle;
#endif

out vec3 nor;
out vec2 uv;
out vec2 uv_bis;
out vec4 color;
#ifdef GL_ARB_bindless_texture
flat out uvec2 handle;
#endif


mat4 getWorldMatrix(vec3 translation, vec3 rotation, vec3 scale);
mat4 getInverseWorldMatrix(vec3 translation, vec3 rotation, vec3 scale);

void main(void)
{
    mat4 ModelMatrix = getWorldMatrix(Origin, Orientation, Scale);
    mat4 ModelViewProjectionMatrix = RSMMatrix * ModelMatrix;
    mat4 TransposeInverseModel = transpose(inverse(ModelMatrix));
    gl_Position = ModelViewProjectionMatrix * vec4(Position, 1.);
    nor = (TransposeInverseModel * vec4(Normal, 0.)).xyz;
    uv = Texcoord;
    color = Color.zyxw;
#ifdef GL_ARB_bindless_texture
    handle = Handle;
#endif
}
