layout (std140) uniform MatrixesData
{
    mat4 ViewMatrix;
    mat4 ProjectionMatrix;
    mat4 InverseViewMatrix;
    mat4 InverseProjectionMatrix;
    mat4 ShadowViewProjMatrixes[4];
    vec2 screen;
};

uniform mat4 ModelMatrix;
uniform mat4 InverseModelMatrix;

uniform mat4 TextureMatrix =
    mat4(1., 0., 0., 0.,
         0., 1., 0., 0.,
         0., 0., 1., 0.,
         0., 0., 0., 1.);


in vec3 Position;
in vec2 Texcoord;
in vec3 Normal;
out vec3 nor;
out vec2 uv;


void main(void)
{
    mat4 ModelViewProjectionMatrix = ShadowViewProjMatrixes[2] * ModelMatrix;
    mat4 TransposeInverseModelView = transpose(InverseModelMatrix * inverse(ShadowViewProjMatrixes[2]));
    gl_Position = ModelViewProjectionMatrix * vec4(Position, 1.);
    nor = Normal;
    uv = (TextureMatrix * vec4(Texcoord, 1., 1.)).xy;
}
