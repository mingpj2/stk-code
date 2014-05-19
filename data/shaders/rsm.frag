uniform sampler2D tex;

in vec2 uv;
in vec3 nor;
out vec3 RSMColor;
out vec3 RSMNormals;

void main()
{
    RSMColor = texture(tex, uv).xyz;
    RSMNormals = nor;
}
