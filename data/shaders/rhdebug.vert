uniform vec3 extents;

layout (std140) uniform MatrixesData
{
    mat4 ViewMatrix;
    mat4 ProjectionMatrix;
    mat4 InverseViewMatrix;
    mat4 InverseProjectionMatrix;
    mat4 ShadowViewProjMatrixes[4];
    vec2 screen;
};

#define DIM_X 32
#define DIM_Y 32
#define DIM_Z 32

out vec3 uvw;

void main(void)
{
    // Determine the RH center
    float   gx = gl_VertexID & 31;
    float   gy = (gl_VertexID >> 5) & 31;
    float   gz = (gl_VertexID >> 10) & 31;
    uvw = vec3(gx, gy, gz) / 32.;
    gx -= 16.;
    gy -= 16.;
    gz -= 16.;
    vec3  stratum = extents / vec3(DIM_X - 1., DIM_Y - 1., DIM_Z - 1.);
    gl_Position = ProjectionMatrix * ViewMatrix * vec4(vec3(gx, gy, gz) * stratum, 1.);
    gl_PointSize = 500. / gl_Position.w;

}